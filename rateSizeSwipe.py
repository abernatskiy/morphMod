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

# Tunable hyperparameters
numTrials = 100
segments = 3
# Optional definitions for pbsGridWalker that depend on the number of segments
pointsPerJob = 20
queue = 'shortq'
expectedWallClockTime = '01:00:00'
maxJobs = 8

# Constant hyperparameters
evsDefaults = {'individual': 'compositeFixedProbabilities', 'evolver': 'cluneSimplifiedMorphologyControlIndividuals', 'communicator': 'chunkedUnixPipe',
               'compositeClass0': 'integerVectorSymmetricRangeMutations', 'probabilityOfMutatingClass0': 0.2,
               'lengthClass0': segments, 'initLowerLimitClass0': 0, 'initUpperLimitClass0': segments, 'lowerCapClass0': 0, 'upperCapClass0': segments,
               'mutationAmplitudeClass0': 1,
               'compositeClass1': 'integerWeightsSwitchableConnections',
               'lengthClass1': 2*(segments**2), 'initLowerLimitClass1': -1, 'initUpperLimitClass1': 1, 'lowerCapClass1': -1, 'upperCapClass1': 1,
               'mutExplorationClass1': 0.8, 'mutInsDelRatioClass1': 1, 'mutationAmplitudeClass1': 1,
               'genStopAfter': 125, 'populationSize': 25,
               'initialPopulationType': 'sparse', 'secondObjectiveProbability': 1.,
               'logParetoFront': 'yes', 'logBestIndividual': 'yes', 'logParetoFrontKeepAllGenerations': 'yes', 'logParetoFrontPeriod': 1, 'logParetoSize': 'yes',
               'backup': 'no', 'trackAncestry': 'no',
               'randomSeed': 0}
arrowbotsDefaults = {'segments': segments, 'sensorAttachmentType': 'variable',
                     'simulationTime': 3., 'timeStep': 0.05,
                     'integrateError': 'false', 'writeTrajectories': 'false'}
arrowbotInitialConditions = [[0]*segments]*segments # segmentsXsegments null matrix
arrowbotTargetOrientations = [ [1 if i==j else 0 for i in range(segments)] for j in range(segments) ] # segmentsXsegments identity matrix
# Optional definitions for pbsGridWalker that are constant
involvedGitRepositories = mmr.involvedGitRepositories
# dryRun = False

### Required pbsGridWalker definitions
computationName = 'rateSizeSwipe_N' + str(segments)

evsDefaults.pop('probabilityOfMutatingClass0')
evsDefaults.pop('initialPopulationType')
evsDefaults.pop('populationSize')
evsDefaults.pop('randomSeed')
nonRSGrid = gr.LinGrid('probabilityOfMutatingClass0', 0.0, 0.2, 0, 5) * \
            gr.Grid1d('initialPopulationType', ['sparse', 'random']) * \
            gr.Grid1d('populationSize', [15, 25, 40, 60])
parametricGrid = nonRSGrid*numTrials + gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=len(nonRSGrid)*numTrials)

def prepareEnvironment(experiment):
	gccommons.prepareEnvironment(experiment)

def runComputationAtPoint(worker, params):
	return gccommons.runComputationAtPoint(worker, params,
		evsDefaults, arrowbotsDefaults,
		arrowbotInitialConditions,
		arrowbotTargetOrientations)

def processResults(experiment):
	import os
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)

	# We'll take a look at some parameters vs relative mutation rate at several stages (generation counts) along the evolutionary process

	# Linear stages
#	stagesToConsider = 5
#	stages = tal.splitIntegerRangeIntoStages(0, evsAdditionalParams['genStopAfter'], stagesToConsider)
	# Exponential stages
	stages = [0, 5]
	mult = 5
	while stages[-1] <= evsDefaults['genStopAfter']:
		stages.append(stages[-1]*mult)
	stages.pop()

	##### Extracting and plotting the distance to the maximally modular morphology (MMM) for various values relative mutation rate #####
	# mmmmdist and similar abbreviations stand for "minimal distance to the maximally modular morphology" (across the Pareto front)

	xlabel = r'$P_{mm}$'
	fieldNames = [ 'gen {}'.format(st) for st in stages ]

	def generateMinMMMDistTimeSlices(gridPoint):
		return [ gctools.minParetoFrontHammingDistanceToMMM(gen) for gen in stages ]
	tplt.plotComputationVariableAgainstParameter(experiment, 'mmmmd', generateMinMMMDistTimeSlices, 'probabilityOfMutatingClass0',
	                                     fieldNames=fieldNames, xlabel=r'$P_{mm}$', ylabel=r'$\mu$')

	def generateFitnessTimeSlices(gridPoint):
		bestIndividualData = np.loadtxt('bestIndividual{}.log'.format(gridPoint['randomSeed']))
		fitnessData = []
		for genRec in range(bestIndividualData.shape[0]):
			gen = int(bestIndividualData[genRec,0])
			if gen in stages:
				fitnessData.append(bestIndividualData[genRec,1]) # WILL break if the best individual records are not in the order of increasing generation
		return fitnessData
	tplt.plotComputationVariableAgainstParameter(experiment, 'error', generateFitnessTimeSlices, 'probabilityOfMutatingClass0',
	                                     fieldNames=fieldNames, transform=lambda x: -1.*x, yscale='log', xlabel=r'$P_{mm}$', ylabel=r'$E$')
