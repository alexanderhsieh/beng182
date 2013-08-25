'''
Oct 18,2012

This is the local alignment algorithm using affine gap penalties
For BENG182

'''
import sys
from datetime import datetime

class ScoreParam:
	def __init__(self, g_o, g_e, match, mismatch):
		self.g_o = g_o
		self.g_e = g_e
		self.match = match
		self.mismatch = mismatch
	def matchchar(self, a,b):
		assert len(a) == len(b) == 1
		if a==b:
			return self.match
		else:
			return self.mismatch
	def __str__(self):
		return "match = %d; mismatch = %d; gap_o = %d, g_e = %d" % (
			self.match, self.mismatch, self.g_o, self.g_e
		)

def main():
	Infinity = float('inf')
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
	affine_align(s, t)
	print datetime.now()


##function to calculate best score and its location, aligning 2 sequences, given score parameters
def affine_align(seqA = [], seqB = [], score=ScoreParam(-2, -2, 1, -3)):
	Infinity = float('inf')
	##create score matrix, deletion matrix, insertion matrix
	S = make_matrix(len(seqA)+1, len(seqB)+1)
	D = make_matrix(len(seqA)+1, len(seqB)+1)
	I = make_matrix(len(seqA)+1, len(seqB)+1)
	choice = make_matrix(len(seqA)+1, len(seqB)+1)
	best = 0
	optloc = (0,0)
	for i in xrange(1, len(seqA)+1):
		S[i][0] = -Infinity
		I[i][0] = -Infinity
		D[i][0] = score.g_o + i*score.g_e
	for i in xrange(1, len(seqB)+1):
		S[0][i] = -Infinity
		I[0][i] = score.g_o + i*score.g_e
		D[0][i] = -Infinity
	for j in xrange(1, len(seqB)+1):
		for i in xrange(1, len(seqA)+1):
			#tmp variable to hold the cost of opening+extending, just so I dont recalculate 100000x
			tmp = score.g_o + score.g_e
			D[i][j] = max(
				D[i-1][j] + score.g_e,
				S[i-1][j] + tmp
			)
			I[i][j] = max(
				I[i][j-1] + score.g_e,
				S[i][j-1] + tmp
			)	
			S[i][j] =max(
				0,
				D[i][j],
				I[i][j],
				S[i-1][j-1] + score.matchchar(seqA[i-1], seqB[j-1])
			)
			if S[i][j] >= best:
				best = S[i][j]
				optloc = (i,j)
	print best, optloc
	return best, optloc
	
##function to create size x by size y matrix filled with 0s
def make_matrix(sizex, sizey):
	return [[0]*sizey for i in xrange(sizex)]
if __name__ == '__main__':
	main()
