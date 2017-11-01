# Curated Hepatitis B Virus Alignments from GenBank Data
# Copyright (C) 2017 University of the Witwatersrand, Johannesburg, South Africa
# Author: Dr Trevor G. Bell, TrevorGrahamBell@gmail.com

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# This file is part of the pipeline described in the following publication:

# Bell, T. G., Yousif, M. and Kramvis, A. 2016
# Bioinformatic curation and alignment of genotyped hepatitis B virus (HBV)
# sequence data from the GenBank public database.
# SpringerPlus 5: 1896.
# https://springerplus.springeropen.com/articles/10.1186/s40064-016-3312-0

#!/usr/bin/python

import climb
import sys
import re

S = climb.Sequence()
S.load(sys.argv[1])

threshold = int(sys.argv[2])

length = int(sys.argv[3])

cons = ''
for i in S.baseDistribution('1-%i' % (S.seqLength()[0][1])):
        highNum = 0
        cc = 0
        for j in i:
                if (j != climb.GAP) and (i[j] > highNum):
                        highNum = i[j]
                        highBase = j
        cons += highBase

print 'ID,Mismatches,Bases,Percentage'

mismatchList = []
overLength = 0

for i in S.seq:
        pos = 0
        bases = 0
        mismatches = 0

        for j in i['seq']:
                if j != climb.GAP:
                        bases += 1
                        if pos <= len(cons)-1:                    # some sequences are longer than the consensus
                                if j != cons[pos]:
                                        mismatches += 1
                pos += 1

        if mismatches/float(bases)*100 > threshold:
                mismatchList.append(i['id'])
                f = open('_%s.fasta' % (i['id'][:8]), 'w')

                # f.write('>Cons\n')                    # to include consensus with query in the file (for muscle)
                # f.write('%s\n' % cons)

                f.write('>%s\n' % (i['id'][:8]))
                f.write('%s\n' % (i['seq']))
                f.close()

        print '%s,%04i,%04i,%3.2f' % (i['id'][:8], mismatches, bases, mismatches/float(bases)*100)

f = open('%s_cons.fasta' % (sys.argv[1][0]), 'w')          # first character is the genotype
f.write('>Cons\n')
f.write('%s\n' % cons)
f.close()

f = open('%s_log.txt' % (sys.argv[1][0]), 'a')          # first character is the genotype
f.write('splitcons:\n')
f.write('Threshold: %i\n' % threshold)
f.write('\tMismatch List:\n')
for i in mismatchList:
        f.write('\t%s\n' % i)

f.write('Number excluded: %i\n' % len(mismatchList))
f.write('of which overlengths: %i\n' % overLength)
f.write('Sequence count before: %i\n' % S.seqCount())
S.seqRemoveByID(mismatchList, escape=True)
f.write('Sequence count after: %i\n' % S.seqCount())
f.close()
S.save(S.fileName, overWrite=True)
S.unload()

# Ends
