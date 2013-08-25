'''
Oct 18, 2012

This is the smith-waterman local alignment algorithm 
Fixed gap penalties
For BENG182
'''

import sys
from datetime import datetime

class ScoreParam:
	def __init__(self, gap, match, mismatch):
		self.gap = gap
		self.match = match
		self.mismatch = mismatch

def main():
	
	seq1 = ''
	seq2 = ''
	s = []
	t = []
	
	##get concatenated sequences as strings seq1, seq2
	file1 = open('human.txt', 'r')
	pulled1 = file1.readlines()
	for line in pulled1:
		current = line.strip()
		if not current.startswith('>'):
			seq1 = seq1 + current
	file2 = open('mouse.txt', 'r')
	pulled2 = file2.readlines()
	for line in pulled2:
		current = line.strip()
		if not current.startswith('>'):
			seq2 = seq2 + current
	##convert to list
	s = list(seq1)
	t = list(seq2)
	
	print datetime.now()
	##print best score and location in matrix
	print local_align(s, t)
	print datetime.now()


##function to calculate best score and its location, aligning 2 sequences, given score parameters
def local_align(seqA = [], seqB = [], score=ScoreParam(-1, 1, -1)):
	##create score matrix
	S = make_matrix(len(seqA)+1, len(seqB)+1)
	best = 0
	optloc = (0,0)
	for j in xrange(1, len(seqB)):
		for i in xrange(1, len(seqA)):
			S[i][j] = max(
				S[i][j-1] + score.gap,
				S[i-1][j] + score.gap,
				S[i-1][j-1] + (score.match if seqB[j] == seqA[i] else score.mismatch),
			)
			if S[i][j] >= best:
				best = S[i][j]
				optloc = (i+1, j+1)
	return best, optloc

##function to create size x by size y matrix filled with 0s
def make_matrix(sizex, sizey):
	return [[0]*sizey for i in xrange(sizex)]
	

if __name__ == '__main__':
	main()
