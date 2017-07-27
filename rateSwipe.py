from os.path import join, exists, expanduser
from time import sleep
import subprocess
import sys
sys.path.append(join(expanduser('~'), 'morphMod'))

import pbsGridWalker.grid as gr
import pbsGridWalker.tools.algorithms as tal
import pbsGridWalker.tools.iniWriter as tiniw
import pbsGridWalker.tools.fsutils as tfs
import pbsGridWalker.tools.fileIO as tfio # for writeColumns()


import morphModRoutes as mmr
import classifiers
import gctools

#dryRun = False

# Tunable hyperparameters
numTrials = 20
segments = 3
computationName = 'rateSwipe_N' + str(segments)

# Constant hyperparameters
evsAdditionalParams = {'individual': 'compositeFixedProbabilities', 'evolver': 'cluneSimplifiedMorphologyControlIndividuals', 'communicator': 'chunkedUnixPipe',
                       'compositeClass0': 'integerVectorSymmetricRangeMutations',
                       'lengthClass0': segments, 'initLowerLimitClass0': 0, 'initUpperLimitClass0': segments, 'lowerCapClass0': 0, 'upperCapClass0': segments,
                       'mutationAmplitudeClass0': 1,
                       'compositeClass1': 'integerWeightsSwitchableConnections',
                       'lengthClass1': 2*(segments**2), 'initLowerLimitClass1': -1, 'initUpperLimitClass1': 1, 'lowerCapClass1': -1, 'upperCapClass1': 1,
                       'mutExplorationClass1': 0.8, 'mutInsDelRatioClass1': 1, 'mutationAmplitudeClass1': 1,
                       'genStopAfter': 200, 'populationSize': 50, 'secondObjectiveProbability': 1.,
                       'logParetoFront': 'yes', 'logBestIndividual': 'yes', 'logParetoFrontKeepAllGenerations': 'yes', 'logParetoFrontPeriod': 1,
                       'backup': 'no', 'trackAncestry': 'no'}
arrowbotsAdditionalParams = {'segments': segments, 'sensorAttachmentType': 'variable',
                             'simulationTime': 3., 'timeStep': 0.05, 'integrateError': 'false', 'writeTrajectories': 'false'}
arrowbotInitialConditions = [[0]*segments]*segments # segmentsXsegments null matrix
arrowbotTargetOrientations = [ [1 if i==j else 0 for i in range(segments)] for j in range(segments) ] # segmentsXsegments identity matrix

# Optional definitions for pbsGridWalker that depend on the number of segments
pointsPerJob = 1
queue = 'shortq'
expectedWallClockTime = '01:00:00'

# Optional definitions for pbsGridWalker that are constant
maxJobs = 2
involvedGitRepositories = mmr.involvedGitRepositories

# Required pbsGridWalker definitions
_morphologicalMutationProbabilityGrid = gr.LinGrid('probabilityOfMutatingClass0', 0.1, 0.1, 0, 8)
_initialPopulationTypeGrid = gr.Grid1d('initialPopulationType', ['sparse', 'random'])
parametricGrid = _morphologicalMutationProbabilityGrid * \
                 _initialPopulationTypeGrid * \
                 gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=numTrials)

def prepareEnvironment(experiment):
	if not exists(mmr.arrowbotsExecutable):
		raise RuntimeError('Arrowbots executable not found at ' + mmr.arrowbotsExecutable)
	if not exists(mmr.evsExecutable):
		raise RuntimeError('EVS executable not found at ' + mmr.evsExecutable)

def runComputationAtPoint(worker, params):
	print('Running evs-arrowbots pair with the following parameters: ' + str(params))
	parsedParams = tal.classifyDictWithRegexps(params, classifiers.serverClientClassifier)
	serverParams = tal.sumOfDicts(parsedParams['server'], evsAdditionalParams)
	print('Server params: ' + str(serverParams))
	clientParams = tal.sumOfDicts(parsedParams['client'], arrowbotsAdditionalParams)
	print('Client params: ' + str(clientParams))
	tiniw.write(serverParams, classifiers.evsClassifier, 'evs.ini')
	tiniw.write(clientParams, classifiers.arrowbotsClassifier, 'arrowbot.ini')
	tfio.writeColumns(arrowbotInitialConditions, 'initialConditions.dat')
	tfio.writeColumns(arrowbotTargetOrientations, 'targetOrientations.dat')

	geneFifo = tfs.makeUniqueFifo('.', 'genes')
	evalFifo = tfs.makeUniqueFifo('.', 'evals')

	clientProc = worker.spawnProcess([mmr.arrowbotsExecutable, geneFifo, evalFifo])
	if not worker.runCommand([mmr.evsExecutable, evalFifo, geneFifo, str(serverParams['randomSeed']), 'evs.ini']):
		return False
	worker.killProcess(clientProc, label='client')
	# TODO: Validation of the obtained files here
	return True

def processResults(experiment):
	#import matplotlib
	#matplotlib.use('Agg')
	#import matplotlib.pyplot as plt
	import os
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)

	##### Extracting and plotting the distance to the maximally modular morphology (MMM) for various values relative mutation rate #####
	# mmmmdist and similar abbreviations stand for "minimal distance to the maximally modular morphology" (across the Pareto front)
	def plotMinMMMDistVSMorphologicalMutationRate():
		stagesToConsider = 5
		stages = tal.splitIntegerRangeIntoStages(1, evsAdditionalParams['genStopAfter'], stagesToConsider)
		def minMMMDistFileName(gridPoint, generation):
			return '../results/IP{}_gen{}_MMP{}_minMMMDist'.format(gridPoint['initialPopulationType'], generation, gridPoint['probabilityOfMutatingClass0'])
		def generateMinMMMDistTimeSlices(gridPoint):
			for gen in stages:
				with open(minMMMDistFileName(gridPoint, gen), 'a') as file:
					file.write('{}\n'.format(gctools.minParetoFrontHammingDistanceToMMM(gen)))
		experiment.executeAtEveryGridPointDir(generateMinMMMDistTimeSlices)

		os.chdir('results')

		mmprobslist = sorted([ gp['probabilityOfMutatingClass0'] for gp in _morphologicalMutationProbabilityGrid ])
		iptypeslist = [ gp['initialPopulationType'] for gp in _initialPopulationTypeGrid ]

		mmmmddata = {}
		for ipt in iptypeslist:
			mmmmddata[ipt] = {}
			for gen in stages:
				mmmmddata[ipt][str(gen)] = []
				for mmprob in mmprobslist:
					fn = minMMMDistFileName({'initialPopulationType': ipt, 'probabilityOfMutatingClass0': mmprob}, gen)
					mmmmddata[ipt][str(gen)].append(np.loadtxt(fn, dtype=np.int))
					os.remove(fn)
				mmmmddata[ipt][str(gen)] = np.stack(mmmmddata[ipt][str(gen)]).T

		def plotMMMMDistVSMorphMod(initPopType):
			filename = 'mmmmdist_vs_mmrate_IP{}.png'.format(initPopType)
			title = 'Minimal distance to maximally modular morphology\nvs morphological mutation rate ({} initial population)'.format(initPopType)
			margins = 0.5

			tplt.plotAverageTimeSeries(mmmmddata[initPopType], 'mmmdist', filename, title=title, legendLocation=1, xlabel='mmrate', xlimit=1, ylimit=None, margins=margins)

		map(plotMMMMDistVSMorphMod, iptypeslist)

		os.chdir('..')

	plotMinMMMDistVSMorphologicalMutationRate()
