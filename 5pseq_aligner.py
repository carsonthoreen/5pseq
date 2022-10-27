# 5pseq_aligner
# Author: Carson Thoreen
# Date: Thu Oct  7 10:56:16 2021
# Script for extracting TSS frequencies from the 5pseq library reads.
# Args: filename of FASTQ file containing 5pseq reads.
# Script will iterate through fastq reads, align them to the promoter sequence defined in main(), 
# and report counts of each sequence, along with other alignment info (pos, identity, etc.).

import sys
import gzip
import re
from operator import itemgetter
from collections import namedtuple

Alignment = namedtuple('Alignment', 'query alignment_str position identity')

PLASMID_REFERENCE_SEQ = 'AGAGCTGGTTTAGTGAACCGTCNNNNNNNGCCGCAGCCGCCGCCATCGTCGACGCGCGCTTCCCTGTTCACCTCTG'
SPIKE_IN = 'GGGGCTCTTCCCATGGCCGCAGCCGGCCGCCATCGTCGACGCGCGCTTCCCTGTTCACCT'

class AlignerSimple : 
	"""
	A simple sequence aligner that searches for a seed sequence in the query rather than
	analyzing identity at all positions. 
	Args: reference sequence, seed region coordinate start, seed region coordinate end
	Returns: A named tuple describing the alignment
	"""

	def __init__(self, reference, seed_start, seed_end) : 
		self.reference = reference
		self.seed_start = seed_start 
		self.seed_end = seed_end 
		self.seed_seq = self.reference[ self.seed_start:self.seed_end ]

	def align( self, query ) :
			
		pos = query.find( self.seed_seq )
		if pos < 0 : return None

		# determine alignment position as coordinates for reference sequence.
		qpos = self.seed_start - pos
		alignment_str  = self.get_alignment_string( query, qpos )
		identity = alignment_str.count('+')

		return Alignment( position=qpos, alignment_str=alignment_str, query=query, identity=identity )
		
	def get_alignment_string(self, query, pos) : 
	
		alignment = [ '+' if (query[i] == self.reference[i+pos] or self.reference[i+pos] == 'N')
								else '-' for i in range(0, len(query)) ]

		return ''.join(alignment)

def main() :

	filename = sys.argv[1]

	maxreads = 100000000
	query_size = 30
	nreads = 0

	# initalize both aligners to search for a 10 nt seed sequence.
	aligner = AlignerSimple( PLASMID_REFERENCE_SEQ, 34, 44 )
	spike_aligner = AlignerSimple( SPIKE_IN, 5,15 )
	alignments = dict()

	with gzip.open( filename, "rt" ) as fh :
	
		for line in fh : 

			if line.startswith('@') : line = fh.readline().strip()
			else: continue

			query = line[0:query_size]

			# ignore reads with ambiguous bases at the start
			if 'N' in query[0:10] : continue

			# test for alignment to spike-in
			a = spike_aligner.align( query )
			if a : 
				if 'SPIKE_IN' in alignments : 
					alignments['SPIKE_IN'][1] += 1
				else : 
					alignments['SPIKE_IN'] = [a,1]
			else : 
				
				a = aligner.align( query )
				if not a : continue

				# use only the 'random' segment of the read plus any 5p extensions
				# for the key so that mismatches later in the read don't affect grouping.
				key = query[0:(7+(22-a.position))] 

				if key in alignments : 
					alignments[key][1] += 1
				else : 
					alignments[key] = [a, 1]

			nreads += 1
			if nreads > maxreads : break

	output_file = filename + '.5pseqs'
	write_output( alignments, output_file )

def write_output(alignments, output_file ) :

	with open(output_file, "w") as ofh : 
		
		header = ('Seq', 'Match', 'Pos', 'Count')
		ofh.write( "\t".join(header) + "\n")

		for key, val in alignments.items() : 
			a, count = val

			query_str = a.query
			if key=='SPIKE_IN' : query_str = 'SPIKE_IN'

			ofh.write( "%s\t%s\t%d\t%d\n" % (query_str, a.alignment_str, a.position, count) )

if __name__=='__main__' : 
	main()

