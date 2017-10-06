#!/usr/bin/python2

# A script that generates Arrowbot genomes with random morphology and a hand-designed controller with a guaranteed performance

import numpy as np

id = 0

def randomMorphology(n):
	return np.random.randint(n+1, size=n)

def jMatrix(morph):
	n = morph.shape[0]
	j = np.zeros((n, n), dtype=np.int)
	for i,sp in enumerate(morph):
		if sp > 0:
			j[i][sp-1] = 1
	return j

def goodController(morph):
	n = morph.shape[0]
	j = jMatrix(morph)
	w = np.eye(n, dtype=np.int)
	k = np.array([ [ 1 if i<=l else 0 for i in range(n) ] for l in range(n) ], dtype=np.int)
	y = (j-w).dot(k)
	return w,y

def goodGenomeString(morph):
	global id
	s = '{} '.format(id)
	id += 1
	s += ' '.join(map(str, morph)) + ' '
	w,y = goodController(morph)
	n = morph.shape[0]
	s += ' '.join(map(str, w.reshape((n*n,)))) + ' '
	s += ' '.join(map(str, y.reshape((n*n,))))
	return s

n = 10
sampleSize = 100000
for _ in range(sampleSize):
	m = randomMorphology(n)
	print goodGenomeString(m)
