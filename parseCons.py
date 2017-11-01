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
import glob

S = climb.Sequence()

geno = sys.argv[1]
threshold = int(sys.argv[2])

files = glob.glob('n_*.fasta')

files.sort()

withGaps = 0
withMismatches = 0
for fn in files:
        S.load(fn)
        cons = S.seq[1]['seq']
        i = S.seq[0]

        m = re.search(climb.GAP, cons)
        if m is not None:
                # print 'Gaps found in consensus aligned with %s' % (i['id'])
                withGaps += 1
                S.unload()
                continue

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

        # print '%s*,%04i,%04i,%3.2f' % (i['id'][:8], mismatches, bases, mismatches/float(bases)*100)

        if mismatches/float(bases)*100 <= threshold:
                withMismatches += 1
                print '>%s' % (i['id'])
                print '%s' % (i['seq'])
        else:
                withGaps += 1

        S.unload()

f = open('%s_log.txt' % (geno), 'a')
f.write('parsecons:\n')
f.write('sequences with gaps (excluded): %i\n' % withGaps)
f.write('sequence with mismatches (included): %i\n' % withMismatches)
f.close()

# Ends

