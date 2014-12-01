#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 27 09:52:49 2014

@author: tgibbons
"""

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC

# Sequence statistics
ss = dict()

fOut1 = open(sys.argv[3], 'w')
fOut2 = open(sys.argv[4], 'w')

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    if GC(rec.seq) < float(sys.argv[2]):
        fOut1.write('>{0}\n{1}\n'.format(rec.description, rec.seq))
    else:
        fOut2.write('>{0}\n{1}\n'.format(rec.description, rec.seq))

fOut1.close()
fOut2.close()

sys.exit()
