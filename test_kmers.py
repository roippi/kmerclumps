#!/usr/bin/python

import unittest
from kmers import sliding_window
from kmers import get_clumps

genome = 'ACTAGACATGAACATGAACATACTCACCAGACATACTACT'

class TestKmerSequence(unittest.TestCase):
    def test_15_3_clump(self):
        k = 3
        L = 15
        t = 3

        clumps = get_clumps(genome, k,L,t)

        self.assertEqual(clumps, {tuple('ACA'), tuple('CAT')})

    def test_t(self):
        k = 3
        L = 15

        clumps_2 = {('A', 'C', 'A'),('T', 'G', 'A'),('A', 'A', 'C'),
                    ('A', 'C', 'T'),('C', 'A', 'T'),('G', 'A', 'A'),
                    ('A', 'T', 'G'),('T', 'A', 'C')}

        self.assertEqual(get_clumps(genome, k,L,2), clumps_2)

        clumps_3 = {tuple('ACA'), tuple('CAT')}
        self.assertEqual(get_clumps(genome, k,L,3), clumps_3)

        clumps_4 = set()
        self.assertEqual(get_clumps(genome, k,L,4), clumps_4)

    def test_k(self):
        L = 15
        t = 3

        clumps_2 = {tuple('CA'), tuple('AT'), tuple('AC'), tuple('GA')}
        self.assertEqual(get_clumps(genome, 2,L,t), clumps_2)

        clumps_3 = {tuple('ACA'), tuple('CAT')}
        self.assertEqual(get_clumps(genome, 3,L,t), clumps_3)

        clumps_4 = set()
        self.assertEqual(get_clumps(genome, 4,L,t), clumps_4)

    def test_L(self):
        k=3
        t=3

        clumps_10 = set()
        self.assertEqual(get_clumps(genome, k,10,t), clumps_10)

        clumps_15 = {tuple('ACA'), tuple('CAT')}
        self.assertEqual(get_clumps(genome, k,15,t), clumps_15)

        clumps_25 = {('A', 'C', 'A'),('C', 'A', 'T'),('T', 'A', 'C'),
                    ('A', 'C', 'T')}

        self.assertEqual(get_clumps(genome, k,25,t), clumps_25)
        

if __name__ == "__main__":
    unittest.main()
