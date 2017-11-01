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

import re
import sys

result_handle = open(sys.argv[1])

from Bio.Blast import NCBIXML
blast_records = NCBIXML.parse(result_handle)

# if alignment length is less than query length, check second hit

def placeFrag(ss, pp, bb):	# accepts the ACTUAL position (1-indexed) as this is what is specified by the user and passed in for non-wraps positions
	return bb[:pp-1] + ss + bb[pp+len(ss)-1:]

# print "Fragment,QueryLength,FullLength,AlignLength"
for blast_record in blast_records:

        if len(blast_record.alignments) == 0:
                f = open('%s_log.txt' % (sys.argv[1][0]), 'a')          # first letter of filename is genotype
                f.write('No match for: %s\n' % blast_record.query)
                f.close()
                continue
        out = blast_record.alignments[0].hsps[0].align_length
        tag = '='
        if blast_record.query_length > blast_record.alignments[0].hsps[0].align_length:
                if len(blast_record.alignments[0].hsps) > 1:
                        out += blast_record.alignments[0].hsps[1].align_length
                        tag = '+'
                else:
                        tag = '-'

        # print '%-10s,%04i,%-10s,%04i,%-01s' % (blast_record.query, blast_record.query_length, blast_record.alignments[0].title.split()[1], out, tag)

        backbone = '-' * int(sys.argv[2])
        for i in blast_record.alignments[0].hsps:
                backbone = placeFrag(i.query, i.sbjct_start, backbone)

        if backbone.find('N') >= 0:
                flagN = 'N'
        else:
                flagN = '.'

        if re.search('[^ACGTN-]', backbone) == None:
                flagW = '.'
        else:
                flagW = 'W'

        c = 0
        for i in backbone:
                if i == '-':
                        c += 1
        c = len(backbone) - c

        if c == int(sys.argv[2]):
                flagL = 'C'
        else:
                flagL = 'S'

        flagUn = blast_record.query[-1] # saved from GBGeno2
        if flagUn == 'V':
                flagUn = '.'    # period dropped after BLASTing, so "V" used instead
        flags = '[%04i]_[%s%s%s%s]' % (c, flagL, flagN, flagW, flagUn)
        out = blast_record.query[:-1] + flags

        # print ">%s_%s_[%s]" % (blast_record.query, tag, blast_record.alignments[0].title.split()[1])
        print ">%s" % (out)
        print backbone + '\n'

# Ends

