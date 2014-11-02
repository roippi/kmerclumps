#!/usr/bin/python

from collections import defaultdict, deque, Counter
from itertools import islice

def sliding_window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   

       this function courtesy http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def get_clumps(genome, k, L, t):
    """
    Given a genome, returns all (L,t)-clumps of k-mers.

    See http://stackoverflow.com/a/26695030/2581969 for some explanation if that doesn't mean anything to you.
    """
    kmers = KmerSequence(L-k, t)

    for kmer in sliding_window(genome, k):
        kmers.add(kmer)
        
    return kmers.clumps

class KmerSequence(object):
    __slots__ = ['order', 'counts', 'limit', 'clumps', 't']

    def __init__(self, limit, threshold):
        self.order = deque()
        self.counts = Counter()
        self.limit = limit
        self.clumps = set()
        self.t = threshold

    def add(self, kmer):
        if len(self.order) > self.limit:
            self._remove_oldest()
        self._add_one(kmer)

    def _add_one(self,kmer):
        self.order.append(kmer)
        old_count = self.counts[kmer]
        new_count = old_count + 1
        self.counts[kmer] = new_count

        if new_count == self.t:
            self.clumps.add(kmer)

    def _remove_oldest(self):
        toremove = self.order.popleft()
        old_count = self.counts[toremove]
        self.counts[toremove] -= 1


if __name__ == '__main__':
    import sys

    with open(sys.argv[1]) as f:
        genome = f.read()

    if len(sys.argv) > 2:
        #user specified k,L,t; use those
        k,L,t = map(int, sys.argv[2:])
    else:
        #defaults
        k = 9
        L = 500
        t = 3

    clumps = get_clumps(genome, k,L,t)

    print '({}, {})-clumps of {}-mers found in that file: {}'.format(
        L,
        t,
        k,
        len(clumps)
        )
