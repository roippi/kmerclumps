kmerclumps
==========

A tool for finding (L,t) clumps of k-mers in a genome.

python usage:

    from kmers import get_clumps

    clumps = get_clumps(genome, k, L, t)

command line usage:

    python kmers.py <path_to_genome_file> [k L t]
