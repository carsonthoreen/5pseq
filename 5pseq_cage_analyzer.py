# 5pseq_cage_analyzer.py
# Author: Carson Thoreen
# Date: Thu Oct 27 17:21:02 2022
# Script for quantifying kmers in a BAM file of aligned CAGE reads.
# Args: BAM file of CAGE reads.
# Intended to be run on a BAM of CAGE reads that is subsetted for expected promoter regions.
# Produces a .kmers file that lists counts of all 7mers at the 5' ends of CAGE reads.

import sys
import pysam


def main() :

	filename = sys.argv[1]
	bamreader = pysam.AlignmentFile( filename, "rb" )

	seqs = dict()
	seen = set()
	stats = {'total':0, 'clipped':0, 'unique':0}

	for aln in bamreader:
		qname = aln.query_name
		stats['total'] += 1

		#remove soft clipped bases from 5p end of read (4 is soft_clip)
		#Note: pysam get_reference_sequence() removes soft-clipped bases automatically.
		cigar = aln.cigartuples
		softclip = 0
		if aln.is_reverse : 
			if cigar[-1][0] == 4 : softclip = cigar[-1][1]
		else : 
			if cigar[0][0] == 4 : softclip = cigar[0][1] 

		if softclip > 0 : 
			stats['clipped'] += 1
			continue
			
		#only count each read once -- ignore if seen already
		if not aln.is_unmapped and (qname not in seen) :

			seen.add(qname)
			stats['unique'] += 1

			seq = aln.get_reference_sequence().upper()
			if aln.is_reverse : 
				seq = reverse_complement(seq)

			subseq = seq[0:7]
			if subseq in seqs : seqs[subseq] += 1
			else : seqs[subseq] = 1
			
	output_file = filename + ".kmers"
	with open(output_file, "w") as ofh : 
		ofh.write("seq\tnumreads\n")
		for k,v in seqs.items() : 
			ofh.write( "%s\t%d\n" % (k,v) )

	print("Total: %d" % (stats['total']))
	print("Soft-clipped: {} ({})".format(stats['clipped'], stats['clipped']/stats['total']))
	print("Unique: {} ({})".format(stats['unique'], stats['unique']/stats['total']))

def reverse_complement(seq) : 
	tr = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
	return "".join( [ tr[nt] for nt in seq[::-1] ] )

if __name__=='__main__' : 
	main()
