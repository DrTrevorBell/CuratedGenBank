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

GENO=$1
THRESHOLDSTART=$2
THRESHOLDEND=$3
LENGTH=$4

if [ -z "$GENO" ] || [ -z "$THRESHOLDSTART" ] || [ -z "$THRESHOLDEND" ] || [ -z "$LENGTH" ] ; then

	echo "Specify genotype, cut-off start, cut-off end and length."

else

	rm -v "$GENO"_log.txt
	touch "$GENO"_log.txt
	rm -v HBV_"$GENO"*
	rm -v "$GENO"*clean*
	rm -v "$GENO"*placed*
	rm -v "$GENO"*merge*
	rm -v "$GENO"*.csv
	rm -v "$GENO".xml

	echo "Extracting reference sequence..."
	python getRef.py "$GENO"_full.fasta "$LENGTH" > HBV_"$GENO"
	ret=$?
	if [ $ret -ne 0 ] ; then
		echo "Error with reference sequence length."
		exit $?
	fi

	echo "Generating BLAST reference library..."
	formatdb -i HBV_"$GENO" -p F -t HBV_"$GENO"

	echo "BLASTing..."
	blastn -out "$GENO".xml -outfmt 5 -query "$GENO".fasta -db HBV_"$GENO" -evalue 0.001 -max_target_seqs 5 -gapopen 5 -gapextend 5

	echo "Parsing..."
        python parseXML2b.py "$GENO".xml "$LENGTH" > "$GENO"placed.fasta

	# Loop from here; replace the "sorted" file each time, as it is modified
	for l in `seq $THRESHOLDSTART 1 $THRESHOLDEND` ;
	do

		TEMPFILE=$GENO$l

		echo "Threshold $l..."

		echo "Sorting..."
		python FASTAsort.py "$GENO"placed.fasta > "$TEMPFILE"placedSorted.fasta
		cp "$TEMPFILE"placedSorted.fasta "$TEMPFILE"placedSortedBackup.fasta

		echo "Analyzing consensus..."
		python splitCons.py "$TEMPFILE"placedSorted.fasta "$l" "$LENGTH" > "$TEMPFILE".csv

		echo "Aligning..."
		for i in _*.fasta ; do needle -asequence $i -bsequence "$GENO"_cons.fasta -aformat fasta -outfile n$i -gapopen 10 -gapextend 0.5 ; done

		echo "Parsing consensus..."
	        python parseCons.py "$GENO" "$l" > "$TEMPFILE"mergeBack.fasta

		echo "Building alignment..."
		cat "$TEMPFILE"placedSorted.fasta "$TEMPFILE"mergeBack.fasta > "$TEMPFILE"clean.fasta
		python FASTAsort.py "$TEMPFILE"clean.fasta > "$TEMPFILE"cleanSorted.fasta

		echo "Checking alignment..."
		python CheckFinalAlignment.py "$TEMPFILE"cleanSorted.fasta "$LENGTH" >> CheckFinalAlignment.txt

		echo "Done."

		rm -v _*.fasta
		rm -v n_*.fasta
		rm -v "$GENO"_cons.fasta

	done

fi

# Ends
