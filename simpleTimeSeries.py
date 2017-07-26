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
computationName = 'simpleTimeSeries_N' + str(segments)

# Constant hyperparameters
initialPopulationTypes = ['sparse', 'random']
evsAdditionalParams = {'individual': 'compositeFixedProbabilities', 'evolver': 'cluneSimplifiedMorphologyControlIndividuals', 'communicator': 'chunkedUnixPipe',
                       'compositeClass0': 'integerVectorSymmetricRangeMutations', 'probabilityOfMutatingClass0': 0.2,
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
parametricGrid = gr.Grid1d('initialPopulationType', initialPopulationTypes)*gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=numTrials)

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
	import os
	import shutil
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)

	def gridFileNamePrefix(gridPoint):
		return 'IP{}'.format(gridPoint['initialPopulationType'])

	##### Extracting and plotting fitness time series #####

	def plotFitnessTSs():
		def fitnessFileName(gp):
			return gridFileNamePrefix(gp) + '_fitness'
		def columnExtractor(gp):
			outFile = fitnessFileName(gp)
			subprocess.call('cut -d \' \' -f 2 bestIndividual*.log | tail -n +4 | tr \'\n\' \' \' >> ../results/' + outFile, shell=True)
			subprocess.call('echo >> ../results/' + outFile, shell=True)
		experiment.executeAtEveryGridPointDir(columnExtractor)

		os.chdir('results')

		xlabel = 'Generations'
		xlimit = evsAdditionalParams['genStopAfter']
		ylimit = None
		margins = 0.5

		title = 'Fitness time series for the two types of initial populations'
		dataDict = {x: -1.*np.loadtxt(fitnessFileName({'initialPopulationType': x})) for x in ['random', 'sparse']}

		# Plotting averages in linear scale on y
		yscale = 'lin'

		xscale = 'lin'
		tplt.plotAverageTimeSeries(dataDict, 'Error', 'errorComparisonLinLin.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins)
		xscale = 'log'
		tplt.plotAverageTimeSeries(dataDict, 'Error', 'errorComparisonLogLin.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins)

		# Plotting the trajectory scatter in logarithmic scale on y
		alpha = 0.3
		yscale = 'log'

		xscale = 'lin'
		tplt.plotAllTimeSeries(dataDict, 'Error', 'errorAllTrajectoriesLinLog.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha)
		xscale = 'log'
		tplt.plotAllTimeSeries(dataDict, 'Error', 'errorAllTrajectoriesLogLog.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha)

		os.chdir('..')

	plotFitnessTSs()

	##### Extracting and plotting time series for distance to the maximally modular morphology (MMM) #####
	def plotMinMMMDistTSs():
		def minMMMDistFileName(gridPoint):
			return '../results/' + gridFileNamePrefix(gridPoint) + '_minMMMDist'
		def generateMinMMMDistTimeSeries(gridPoint):
			minMMMDistTS = [ gctools.minParetoFrontHammingDistanceToMMM(gen) for gen in range(1, evsAdditionalParams['genStopAfter']+1) ]
			filename = minMMMDistFileName(gridPoint)
			with open(filename, 'a') as file:
				file.write(' '.join(map(str, minMMMDistTS)) + '\n')
		experiment.executeAtEveryGridPointDir(generateMinMMMDistTimeSeries)

		os.chdir('results')
		xlabel = 'Generations'
		ylimit = None
		xlimit = evsAdditionalParams['genStopAfter']
		margins = 0.5

		title = 'Hamming distance to maximally modular morphology'
		dataDict = {x: np.loadtxt(minMMMDistFileName({'initialPopulationType': x})) for x in ['random', 'sparse']}

		# Plotting averages in linear time scales on y
		yscale = 'lin'
		xscale = 'lin'
		tplt.plotAverageTimeSeries(dataDict, 'Mutations to MMM', 'minMMMDistTS.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins)

		os.chdir('..')

	plotMinMMMDistTSs()
