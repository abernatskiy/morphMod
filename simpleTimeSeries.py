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
numTrials = 4
segments = 3
# Optional definitions for pbsGridWalker that depend on the number of segments
pointsPerJob = 5
maxJobs = 2
queue = 'shortq'
expectedWallClockTime = '03:00:00'

# Constant hyperparameters
evsDefaults = {'individual': 'compositeFixedProbabilities', 'evolver': 'cluneSimplifiedMorphologyControlIndividuals', 'communicator': 'chunkedUnixPipe',
               'compositeClass0': 'integerVectorSymmetricRangeMutations', 'probabilityOfMutatingClass0': 0.2,
               'lengthClass0': segments, 'initLowerLimitClass0': 0, 'initUpperLimitClass0': segments, 'lowerCapClass0': 0, 'upperCapClass0': segments,
               'mutationAmplitudeClass0': 1,
               'compositeClass1': 'integerWeightsSwitchableConnections',
               'lengthClass1': 2*(segments**2), 'initLowerLimitClass1': -1, 'initUpperLimitClass1': 1, 'lowerCapClass1': -1, 'upperCapClass1': 1,
               'mutExplorationClass1': 0.8, 'mutInsDelRatioClass1': 1, 'mutationAmplitudeClass1': 1,
               'genStopAfter': 600, 'populationSize': 50,
               'initialPopulationType': 'random', 'secondObjectiveProbability': 1.,
               'logParetoFront': 'yes', 'logBestIndividual': 'yes', 'logParetoFrontKeepAllGenerations': 'yes', 'logParetoFrontPeriod': 5, 'logParetoSize': 'yes',
               'backup': 'no', 'trackAncestry': 'no',
               'randomSeed': 0}
evsDefaults['logParetoFrontPeriod'] = 1 # overriding because the field in defaults gets changed during the scripts' autogeneration
arrowbotsDefaults = {'segments': segments, 'sensorAttachmentType': 'variable',
                     'simulationTime': 10., 'timeStep': 0.1,
                     'integrateError': 'false', 'writeTrajectories': 'false'}
arrowbotInitialConditions = [[0]*segments]*segments # segmentsXsegments null matrix
arrowbotTargetOrientations = [ [1 if i==j else 0 for i in range(segments)] for j in range(segments) ] # segmentsXsegments identity matrix
# Optional definitions for pbsGridWalker that are constant
involvedGitRepositories = mmr.involvedGitRepositories
# dryRun = False

### Required pbsGridWalker definitions
computationName = 'simpleTimeSeries_N' + str(segments)

gcvsrandGrid = gr.Grid1d('compositeClass0', ['integerVectorSymmetricRangeMutations', 'integerVectorRandomJumps'])*gr.Grid1d('probabilityOfMutatingClass0', [0.2])
constmorphGrid = gr.Grid1d('compositeClass0', ['integerVectorSymmetricRangeMutations'])*gr.Grid1d('probabilityOfMutatingClass0', [0.0])
nonRSGrid = gcvsrandGrid.concatenate(constmorphGrid)
parametricGrid = nonRSGrid*numTrials + gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=len(nonRSGrid)*numTrials)

for par in parametricGrid.paramNames():
	evsDefaults.pop(par)

def prepareEnvironment(experiment):
	gccommons.prepareEnvironment(experiment)

def runComputationAtPoint(worker, params):
	return gccommons.runComputationAtPoint(worker, params,
		evsDefaults, arrowbotsDefaults,
		arrowbotInitialConditions,
		arrowbotTargetOrientations)

def processResults(experiment):
	'''
	import os
	import shutil
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)

	def gridFileNamePrefix(gridPoint):
		return 'IP{}'.format(gridPoint['initialPopulationType'])

	##### Extracting and plotting fitness time series #####

	xlabel = r'$T$'
	figureDims = None
	xlimit = evsDefaults['genStopAfter']
	margins = 0.5
	strips = 'conf95'

	title = None
	legendLocation = None

	def plotFitnessTSs():
		def fitnessFileName(gp):
			return gridFileNamePrefix(gp) + '_fitness'
		def columnExtractor(gp):
			outFile = fitnessFileName(gp)
			subprocess.call('cut -d \' \' -f 2 bestIndividual*.log | tail -n +4 | tr \'\n\' \' \' >> ../results/' + outFile, shell=True)
			subprocess.call('echo >> ../results/' + outFile, shell=True)
		experiment.executeAtEveryGridPointDir(columnExtractor)

		os.chdir('results')

		ylabel = r'$E$'
		ylimit = None

		title = None
		dataDict = {x: -1.*np.loadtxt(fitnessFileName({'initialPopulationType': x})) for x in ['random']}

		# Plotting averages in logarithmic scale on y
		yscale = 'log'

		xscale = 'lin'
		tplt.plotAverageTimeSeries(dataDict, ylabel, 'errorComparisonLinLin.png', title=title, legendLocation=legendLocation, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, strips=strips, figureDims=figureDims)
		xscale = 'log'
		tplt.plotAverageTimeSeries(dataDict, ylabel, 'errorComparisonLogLin.png', title=title, legendLocation=legendLocation, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, strips=strips, figureDims=figureDims)

		# Plotting the trajectory scatter in logarithmic scale on y
		alpha = 0.3
		yscale = 'log'

		xscale = 'lin'
		tplt.plotAllTimeSeries(dataDict, ylabel, 'errorAllTrajectoriesLinLog.png', title=title, legendLocation=legendLocation, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha, figureDims=figureDims)
		xscale = 'log'
		tplt.plotAllTimeSeries(dataDict, ylabel, 'errorAllTrajectoriesLogLog.png', title=title, legendLocation=legendLocation, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha, figureDims=figureDims)

		os.chdir('..')

	plotFitnessTSs()

	##### Extracting and plotting time series for distance to the maximally modular morphology (MMM) #####
	def plotMinMMMDistTSs():
		def minMMMDistFileName(gridPoint):
			return '../results/' + gridFileNamePrefix(gridPoint) + '_minMMMDist'
		def generateMinMMMDistTimeSeries(gridPoint):
			minMMMDistTS = [ gctools.minParetoFrontHammingDistanceToMMM(gen) for gen in range(1, evsDefaults['genStopAfter']+1) ]
			filename = minMMMDistFileName(gridPoint)
			with open(filename, 'a') as file:
				file.write(' '.join(map(str, minMMMDistTS)) + '\n')
		experiment.executeAtEveryGridPointDir(generateMinMMMDistTimeSeries)

		os.chdir('results')
		ylabel = r'$\mu$'
		ylimit = None

		title = None
		dataDict = {x: np.loadtxt(minMMMDistFileName({'initialPopulationType': x})) for x in ['random']}

		# Plotting averages in linear time scales on y
		yscale = 'lin'
		xscale = 'lin'
		tplt.plotAverageTimeSeries(dataDict, ylabel, 'minMMMDistTS.png', title=title, legendLocation=legendLocation, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, strips=strips, figureDims=figureDims)

		os.chdir('..')

	plotMinMMMDistTSs()
	'''
