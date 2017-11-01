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

if len(sys.argv) < 2:
	sys.exit('Specify raw GenBank input file.')

fullLengths = {'A': 3221, 'B': 3215, 'C': 3215, 'D': 3182, 'E': 3212, 'F': 3215, 'G': 3248, 'H': 3215, 'I': 3215, 'J': 3215}

aF = open('A.fasta','w')
bF = open('B.fasta','w')
cF = open('C.fasta','w')
dF = open('D.fasta','w')
eF = open('E.fasta','w')
fF = open('F.fasta','w')
gF = open('G.fasta','w')
hF = open('H.fasta','w')
iF = open('I.fasta','w')
jF = open('J.fasta','w')

aFull = open('A_full.fasta','w')
bFull = open('B_full.fasta','w')
cFull = open('C_full.fasta','w')
dFull = open('D_full.fasta','w')
eFull = open('E_full.fasta','w')
fFull = open('F_full.fasta','w')
gFull = open('G_full.fasta','w')
hFull = open('H_full.fasta','w')
iFull = open('I_full.fasta','w')
jFull = open('J_full.fasta','w')

sumF = open('summary.csv', 'w')
exportF = open('export.csv', 'w')

from Bio import SeqIO
for seq_record in SeqIO.parse(sys.argv[1], 'genbank'):
        out = seq_record.annotations['accessions']
        source = False
        genotype = '?'
        for f in seq_record.features:
                if f.type == 'source':
                        source = True
                        if f.qualifiers.has_key('note'):
                                if (re.search('type', f.qualifiers['note'][0], flags=re.IGNORECASE)) and not (re.search('recomb', f.qualifiers['note'][0], flags=re.IGNORECASE)):
                                        nn = [tt for tt in re.split('[,:; ]', f.qualifiers['note'][0]) if tt]
                                        nflag = False
                                        o = []
                                        for i in nn:
                                                if nflag:       # if 'type' found in previous entry
                                                        o.append(i)
                                                if re.search('type', i, flags=re.IGNORECASE):
                                                        nflag = True
                                                else:
                                                        nflag = False
                                        out.append(o)

                                        # single letter, and then other data
                                        genonotes = []
                                        for i in o:
                                                if len(i) == 1 and i.upper() in ['A','B','C','D','E','F','G','H','I','J'] and genotype == '?':
                                                        genotype = i.upper()
                                                elif genotype == '?' and re.match('^[ABCDEFGHIJ][0-9]', i.upper()):
                                                        genotype = i[0].upper()
                                else:
                                        out.append('NO TYPE')
                        else:
                                out.append('NO NOTE')
        if not source:
                out.append('NO SOURCE')

	exportF.write('%s,%s\n' % (out[0], genotype))	# move into condition above to export only genotyped accession numbers

        if genotype in ['A','B','C','D','E','F','G','H','I','J']:

                if re.search('UNVERIFIED', seq_record.description, flags=re.IGNORECASE) == None:
                        flagUn = 'V'    # '.' is dropped after BLASTing
                else:
                        flagUn = 'U'

                flags = '_[%s]_%s' % (genotype, flagUn) # passing flagUn to parseXML2b later

        if genotype == 'A':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        aF.write('>%s%s\n' % (out[0], flags))
                        aF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        aFull.write('>%s%s\n' % (out[0], flags))
                        aFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'B':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        bF.write('>%s%s\n' % (out[0], flags))
                        bF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        bFull.write('>%s%s\n' % (out[0], flags))
                        bFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'C':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        cF.write('>%s%s\n' % (out[0], flags))
                        cF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        cFull.write('>%s%s\n' % (out[0], flags))
                        cFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'D':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        dF.write('>%s%s\n' % (out[0], flags))
                        dF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        dFull.write('>%s%s\n' % (out[0], flags))
                        dFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'E':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        eF.write('>%s%s\n' % (out[0], flags))
                        eF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        eFull.write('>%s%s\n' % (out[0], flags))
                        eFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'F':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        fF.write('>%s%s\n' % (out[0], flags))
                        fF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        fFull.write('>%s%s\n' % (out[0], flags))
                        fFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'G':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        gF.write('>%s%s\n' % (out[0], flags))
                        gF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        gFull.write('>%s%s\n' % (out[0], flags))
                        gFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'H':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        hF.write('>%s%s\n' % (out[0], flags))
                        hF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        hFull.write('>%s%s\n' % (out[0], flags))
                        hFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'I':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        iF.write('>%s%s\n' % (out[0], flags))
                        iF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        iFull.write('>%s%s\n' % (out[0], flags))
                        iFull.write('%s\n' % seq_record.seq.tostring())

        if genotype == 'J':
                if len(seq_record.seq.tostring()) <= fullLengths[genotype]:
                        jF.write('>%s%s\n' % (out[0], flags))
                        jF.write('%s\n' % seq_record.seq.tostring())
                else:
                        print '%s\t%s\t%i' % (genotype, out[0], len(seq_record.seq.tostring()))
                if len(seq_record.seq.tostring()) == fullLengths[genotype]:
                        jFull.write('>%s%s\n' % (out[0], flags))
                        jFull.write('%s\n' % seq_record.seq.tostring())

        sumF.write('%s,%i\n' % (genotype,len(seq_record.seq.tostring())))

aF.close()
bF.close()
cF.close()
dF.close()
eF.close()
fF.close()
gF.close()
hF.close()
iF.close()
jF.close()

aFull.close()
bFull.close()
cFull.close()
dFull.close()
eFull.close()
fFull.close()
gFull.close()
hFull.close()
iFull.close()
jFull.close()

sumF.close()
exportF.close()

# Ends
