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

S = climb.Sequence()
S.load(sys.argv[1])

order = []
p = 0
for i in S.seq:
        c = 0
        firstBase = None
        pos = 0
        for j in i['seq']:
                if j != '-':
                        c += 1
                if (firstBase is None) and (j != '-'):
                        firstBase = pos
                pos += 1
        order.append([p,c,firstBase,i['id']])
        p += 1

# http://stackoverflow.com/a/4233482/580010
order.sort(key=lambda keys: (keys[2], keys[1]))

for i in order:
        print '>%s' % (i[3])
        print '%s' % (S.seq[i[0]]['seq'])

S.unload()

# Ends
