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

#!/bin/bash

THEFILE=$1
THEDATE=$2

if [ -z "$THEFILE" ] || [ -z "$THEDATE" ] ; then

    echo "Usage: runAll.sh <input file> <dated folder>"

else

    python GBGeno2.py "$THEFILE" > OverLengths.txt

    touch CheckFinalAlignment.txt

    ./RunBlastN.sh A 8 8 3221
    ./RunBlastN.sh B 8 8 3215
    ./RunBlastN.sh C 8 8 3215
    ./RunBlastN.sh D 8 8 3182
    ./RunBlastN.sh E 8 8 3212
    ./RunBlastN.sh F 8 8 3215
    ./RunBlastN.sh G 8 8 3248
    ./RunBlastN.sh H 8 8 3215
    ./RunBlastN.sh I 8 8 3215
    ./RunBlastN.sh J 10 10 3215

    mkdir $THEDATE

    for i in {A..J} ;
    do
	cp "$i"8cleanSorted.fasta "$THEDATE"/Genotype"$i".fasta
	zip "$THEDATE"/Genotype"$i".zip "$THEDATE"/Genotype"$i".fasta
    done

    for i in {A..J} ;
    do
	grep -c ">" "$THEDATE"/Genotype"$i".fasta
	rm "$THEDATE"/Genotype"$i".fasta
    done

    # Tidy up
    zip $THEDATE.zip $THEFILE ?.xml *cleanSorted* ?_log.txt *8.csv OverLengths.txt CheckFinalAlignment.txt formatdb.log export.csv summary.csv

    for i in {A..J} ;
    do
	rm "$i".xml
	rm "$i"_log.txt
	rm "$i"8*.fasta
	rm "$i".fasta
	rm "$i"_full.fasta
	rm "$i"placed.fasta
	rm *8.csv
    done

    rm OverLengths.txt CheckFinalAlignment.txt HBV_* formatdb.log export.csv summary.csv $THEFILE

    mv $THEDATE.zip $THEDATE/

fi

# Ends
