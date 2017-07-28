from os.path import join, expanduser
from time import sleep
import subprocess
import sys
sys.path.append(join(expanduser('~'), 'morphMod'))

import pbsGridWalker.grid as gr
import pbsGridWalker.tools.algorithms as tal
import pbsGridWalker.tools.fsutils as tfs

import morphModRoutes as mmr
import classifiers
import gctools
import gccommons

#dryRun = False

# Tunable hyperparameters
numTrials = 1
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
                       'genStopAfter': 200, 'populationSize': 25, 'secondObjectiveProbability': 1.,
                       'logParetoFront': 'yes', 'logBestIndividual': 'yes', 'logParetoFrontKeepAllGenerations': 'yes', 'logParetoFrontPeriod': 1,
                       'backup': 'no', 'trackAncestry': 'no'}
arrowbotsAdditionalParams = {'segments': segments, 'sensorAttachmentType': 'variable',
                             'simulationTime': 3., 'timeStep': 0.05, 'integrateError': 'false', 'writeTrajectories': 'false'}
arrowbotInitialConditions = [[0]*segments]*segments # segmentsXsegments null matrix
arrowbotTargetOrientations = [ [1 if i==j else 0 for i in range(segments)] for j in range(segments) ] # segmentsXsegments identity matrix

# Optional definitions for pbsGridWalker that depend on the number of segments
pointsPerJob = 20
queue = 'shortq'
expectedWallClockTime = '01:00:00'

# Optional definitions for pbsGridWalker that are constant
maxJobs = 2
involvedGitRepositories = mmr.involvedGitRepositories

# Required pbsGridWalker definitions
_morphologicalMutationProbabilityGrid = gr.LinGrid('probabilityOfMutatingClass0', 0.0, 0.1, 0, 10)
_initialPopulationTypeGrid = gr.Grid1d('initialPopulationType', ['sparse', 'random'])
parametricGrid = _morphologicalMutationProbabilityGrid * \
                 _initialPopulationTypeGrid * \
                 gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=numTrials)

def prepareEnvironment(experiment):
	gccommons.prepareEnvironment(experiment)

def runComputationAtPoint(worker, params):
	return gccommons.runComputationAtPoint(worker, params,
		evsAdditionalParams, arrowbotsAdditionalParams,
		arrowbotInitialConditions,
		arrowbotTargetOrientations)

def processResults(experiment):
	#import matplotlib
	#matplotlib.use('Agg')
	#import matplotlib.pyplot as plt
	import os
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)

	# We'll take a look at some parameters vs relative mutation rate at several stages (generation counts) along the evolutionary process

	# Linear stages
#	stagesToConsider = 5
#	stages = tal.splitIntegerRangeIntoStages(0, evsAdditionalParams['genStopAfter'], stagesToConsider)

	# Exponential stages
	stages = [0, 1]
	mult = 5
	while stages[-1] < evsAdditionalParams['genStopAfter']:
		stages.append(stages[-1]*mult)
	stages.pop()

	##### Extracting and plotting the distance to the maximally modular morphology (MMM) for various values relative mutation rate #####
	# mmmmdist and similar abbreviations stand for "minimal distance to the maximally modular morphology" (across the Pareto front)

	def plotMinMMMDistVSMorphologicalMutationRate():
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
					#os.remove(fn)
				mmmmddata[ipt][str(gen)] = np.stack(mmmmddata[ipt][str(gen)]).T

		def plotMMMMDistVSMorphMod(initPopType):
			filename = 'mmmmdist_vs_mmrate_IP{}.png'.format(initPopType)
			title = 'Minimal distance to maximally modular morphology\nvs morphological mutation rate ({} initial population)'.format(initPopType)
			margins = 0.5

			tplt.plotAverageTimeSeries(mmmmddata[initPopType], 'mmmdist', filename, title=title, legendLocation=1, xlabel='mmrate', xlimit=1, ylimit=None, margins=margins, timeRange=mmprobslist)

		map(plotMMMMDistVSMorphMod, iptypeslist)
		#print(str(mmmmddata))

		os.chdir('..')

	plotMinMMMDistVSMorphologicalMutationRate()

	def plotFitnessVSMorphologicalMutationRate():
		import glob
		def fitnessFileName(gridPoint, generation):
			return '../results/IP{}_gen{}_MMP{}_fitness'.format(gridPoint['initialPopulationType'], generation, gridPoint['probabilityOfMutatingClass0'])
		def generateFitnessTimeSlices(gridPoint):
			bestIndividualFileName = glob.glob('./bestIndividual*.log')[0]
			bestIndividualData = np.loadtxt(bestIndividualFileName)
			for genRec in range(bestIndividualData.shape[0]):
				gen = int(bestIndividualData[genRec,0])
				if gen in stages:
					with open(fitnessFileName(gridPoint, gen), 'a') as file:
						file.write('{}\n'.format(bestIndividualData[genRec,1]))
		experiment.executeAtEveryGridPointDir(generateFitnessTimeSlices)

		os.chdir('results')

		mmprobslist = sorted([ gp['probabilityOfMutatingClass0'] for gp in _morphologicalMutationProbabilityGrid ])
		iptypeslist = [ gp['initialPopulationType'] for gp in _initialPopulationTypeGrid ]

		mmmmddata = {}
		for ipt in iptypeslist:
			mmmmddata[ipt] = {}
			for gen in stages:
				mmmmddata[ipt][str(gen)] = []
				for mmprob in mmprobslist:
					fn = fitnessFileName({'initialPopulationType': ipt, 'probabilityOfMutatingClass0': mmprob}, gen)
					mmmmddata[ipt][str(gen)].append(np.loadtxt(fn))
					#os.remove(fn)
				mmmmddata[ipt][str(gen)] = -1.*np.stack(mmmmddata[ipt][str(gen)]).T

		def plotFitnessVSMorphMod(initPopType):
			filename = 'fitness_vs_mmrate_IP{}.png'.format(initPopType)
			title = 'Fitness vs morphological mutation rate ({} initial population)'.format(initPopType)
			margins = 0.5

			tplt.plotAverageTimeSeries(mmmmddata[initPopType], 'error', filename, title=title, legendLocation=1, xlabel='mmrate', xlimit=1, ylimit=None, margins=margins, timeRange=mmprobslist, yscale='log')

		map(plotFitnessVSMorphMod, iptypeslist)
		#print(str(mmmmddata))

		os.chdir('..')

	plotFitnessVSMorphologicalMutationRate()
