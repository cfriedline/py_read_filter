from IPython.parallel import Client
import sys
import os
import socket
import numpy
import gzip
import tempfile
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import izip
import shutil
import stopwatch
import traceback
from collections import defaultdict, deque
from multiprocessing import Pool, Manager
import multiprocessing
import argparse
import fabric
import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
log.addHandler(ch)

def format_fastq_tuple(title, seq, qual):
    assert len(seq) == len(qual)
    return "@%s\n%s\n+\n%s\n" % (title, seq, qual)


def process_single_file(f, args):
    tmp = tempfile.NamedTemporaryFile(delete=False)
    basename = os.path.basename(f)
    count = 0
    n = 0
    trimmed = 0
    for title, seq, qual in FastqGeneralIterator(open(f)):
        if seq.startswith("N"):
            seq = seq[1:]
            qual = qual[1:]
            
        if "N" not in seq:
            scores = eval_quality(qual, args)
            if scores:
                if len(scores) != len(seq):
                    seq = seq[0:len(scores)]
                    qual = qual[0:len(scores)]
                    trimmed += 1
                tmp.write(format_fastq_tuple(title, seq, qual))
        else:
            n += 1
            
        count += 1
        
        if count % 10000 == 0:
            print "%s, %s, %d, %d, %d" % (socket.gethostname(),
                                          basename,
                                          count,
                                          n,
                                          trimmed)
    tmp.close()
    return tmp.name


def collapse_results(source, results):
    out = source.replace(".gz", "_processed.fastq")
    temp = tempfile.NamedTemporaryFile(delete=False)
    for r in results:
        for line in open(r):
            temp.write(line)
    temp.close()
    shutil.copy(temp.name, out)
    os.remove(temp.name)
    x = [os.remove(x) for x in results]
    return out


def process_single(seqs, args):
    timer = stopwatch.Timer()
    pool = Pool()
    hostname = socket.gethostname()
    splits = split_file(seqs)
    results = []
    source = None

    for k, temp_files in splits.items():
        source = k
        for f in temp_files:
            try:
                p = pool.apply_async(process_single_file, (f,args))
                results.append(p)
#                 print f
#                 process_single_file(f)
            except:
                traceback.print_exc()

    pool.close()
    pool.join()
    
    # collapse processed temp files
    res = collapse_results(source, [x.get() for x in results])
    
    # remove temp split source files
    for k, v in splits.items():
        x = [os.remove(x) for x in v]
    timer.stop()
    return hostname, source, res, timer.elapsed


def split_file(seqs):
    d = defaultdict(list)
    num_cpu = multiprocessing.cpu_count()
    for seq in seqs:
        print "seq=", seq
        f, num = seq
        print f
        reads_per_file = float(num)//num_cpu
        read_idx = 0
        file_num = 0
        for title, seq, qual in FastqGeneralIterator(gzip.open(f)):
            if read_idx == 0:
                t = tempfile.NamedTemporaryFile(delete=False)
                print socket.gethostname(), t.name, file_num + 1, "/", num_cpu
                d[f].append(t)
            t.write(format_fastq_tuple(title, seq, qual))
            read_idx += 1
            
            if read_idx == reads_per_file:
                read_idx = 0
                file_num += 1
    for k, l in d.items():
        [x.close() for x in l]
        d[k] = [x.name for x in l]
    return d


def convert_qual(q):
    return ord(q)-33


def get_qual_scores(q):
    qual = [ord(x)-33 for x in q]  # list comps seems to be fastest here
    return numpy.array([qual, numpy.mean(qual)])


def eval_quality(q, args):
    qual = get_qual_scores(q)
    scores = qual[0]

    if qual[1] < args.qual_cutoff:
        return False
    
    below_cutoff = 0.0
    window = deque(maxlen=args.win_size)
    win_end = win_size
    last_good = None
    for s in scores:
        window.append(s)
        if s < args.qual_cutoff:
            below_cutoff += 1  # keep track of scores below the quality cutoff
        if len(window) == win_size:
            if numpy.mean(window) < qual_cutoff:
                if last_good is None:
                    last_good = win_end
                    if float(last_good)/len(scores) < len_cutoff:
                        return False  # then it's too short
            win_end += 1
    perc_below = below_cutoff/len(scores)

    if last_good:
        # trim the scores if it will be long enough
        scores = scores[0:(last_good-1)]
        
    # perc_len = float(len(scores))/len(qual[0])

    if perc_below > args.qual_perc_cutoff:
        # drop reads if overall bases have quality values < cutoff,
        # even if average is ok
        return False
    return scores


def process_paired_files(file1, file2, queue, args):
    f1 = FastqGeneralIterator(open(file1))
    f2 = FastqGeneralIterator(open(file2))

    tmp1 = tempfile.NamedTemporaryFile(delete=False)
    tmp2 = open(tmp1.name + ".1", "w")

    basename = [os.path.basename(x) for x in [file1, file2]]
    queue.put(basename)
    count = 0
    n = 0
    trimmed = 0
    for (h1, s1, q1), (h2, s2, q2) in izip(f1, f2):
        for pair in [[s1, q1], [s2, q2]]:
            if pair[0].startswith("N"):
                pair[0] = pair[0][1:]
                pair[1] = pair[1][1:]

        if "N" not in s1 and "N" not in s2:
            scores1 = eval_quality(q1, args)
            scores2 = eval_quality(q2, args)
            if scores1 and scores2:
                if len(scores1) != len(s1):
                    s1 = s1[0:len(scores1)]
                    q1 = q1[0:len(scores1)]
                    trimmed += 1
                if len(scores2) != len(s2):
                    s2 = s2[0:len(scores2)]
                    q2 = q2[0:len(scores2)]
                    trimmed += 1
                tmp1.write(format_fastq_tuple(h1, s1, q1))
                tmp2.write(format_fastq_tuple(h2, s2, q2))
        else:
            n += 1

        count += 1
        if count % 10000 == 0:
            queue.put("%s, %s, %d, %d, %d" % (socket.gethostname(),
                                              basename,
                                              count,
                                              n,
                                              trimmed))
    [x.close() for x in [tmp1, tmp2]]
    queue.put("DONE")
    return tmp1.name, tmp2.name


def collapse_paired_results(sources, results):
    outs = []
    for i in xrange(len(sources)):
        out = sources[i].replace(".gz", "_processed.fastq")
        outs.append(out)
        temp = tempfile.NamedTemporaryFile(delete=False)
        for r in [x[i] for x in results]:
            for line in open(r):
                temp.write(line)
        temp.close()
        shutil.copy(temp.name, out)
        os.remove(temp.name)
    for pair in results:
        for p in pair:
            os.remove(p)
    return outs


def process_paired(seqs, args):
    timer = stopwatch.Timer()
    splits = split_file([seqs[0], seqs[1]])
    sources = []
    tmpfiles = []
    pool = Pool()
    manager = Manager()
    queue = manager.Queue()

    for k, v in splits.items():
        sources.append(k)
        tmpfiles.append(v)
        results = []

        pairs = 0

    for temp1, temp2 in izip(tmpfiles[0], tmpfiles[1]):
        p = pool.apply_async(process_paired_files, (temp1, temp2, queue, args))
        pairs += 1
        results.append(p)
    pool.close()
    completed = 0
    while True:
        item = queue.get()
        print item, completed
        if item == "DONE":
            completed += 1
        if completed == pairs:
            break
    pool.join()

    res = collapse_paired_results(sources, [x.get() for x in results])

    for i in xrange(len(tmpfiles)):
        for j in tmpfiles[i]:
            os.remove(j)
    timer.stop()
    return socket.gethostname(), sources, res, timer.elapsed


def setup_cluster_nodes(dview):
    dview['process_paired'] = process_paired
    dview['process_paired_files'] = process_paired_files
    dview['collapse_paired_results'] = collapse_paired_results
    dview['eval_quality'] = eval_quality
    dview['convert_qual'] = convert_qual
    dview['get_qual_scores'] = get_qual_scores
    dview['split_file'] = split_file
    dview['process_single'] = process_single
    dview['collapse_results'] = collapse_results
    dview['process_single_file'] = process_single_file
    dview['format_fastq_tuple'] = format_fastq_tuple

    
def get_args():
    log.info("getting args")
    p = argparse.ArgumentParser(description="Python Read Filterer")
    p.add_argument("--read1")
    p.add_argument("--read2")
    p.add_argument("--cluster_nodes", default=1)
    p.add_argument("--cluster_profile", default="default")
    p.add_argument("--cluster_delay", default=15)
    p.add_argument("--sge", type=bool, default=False)
    p.add_argument("--win_size", default=5)
    p.add_arguymetn("--qual_cutoff", default=30)
    p.add_argument("--len_cutoff", default=0.5)
    p.add_argument("--qual_perc_cutoff", default=0.20)
    
    if len(sys.argv) < 2:
        p.print_help()
        sys.exit()

    args = p.parse_args()
    return args

def start_cluster(args):
    pass

def main():
    args = get_args()
    start_cluster(args)

if __name__ == '__main__':
    main()
