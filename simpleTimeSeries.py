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

#dryRun = False

# Tunable hyperparameters
numTrials = 20
segments = 3
computationName = 'basisScript_N' + str(segments)

# Constant hyperparameters
# Regarding convergence times: about 100 gens for 3 segments at probability of morphological mutation of 0.2, 400 for 5, 10 gets stuck indefinitely

# attachments = {'identity': 'J=I', 'null': 'J=0'}
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
# parametricGrid = gr.Grid1d('sensorAttachmentType', attachments.keys())*gr.Grid1d('initialPopulationType', initialPopulationTypes)*gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=numTrials)
parametricGrid = gr.Grid1d('initialPopulationType', initialPopulationTypes)*gr.Grid1dFromFile('randomSeed', mmr.randSeedFile, size=numTrials)

def prepareEnvironment(experiment):
	if not exists(mmr.arrowbotsExecutable):
		raise RuntimeError('Arrowbots executable not found at ' + mmr.arrowbotsExecutable)
	if not exists(mmr.evsExecutable):
		raise RuntimeError('EVS executable not found at ' + mmr.evsExecutable)

def runComputationAtPoint(worker, params):
	print('Running evs-arrowbots pair with the following parameters: ' + str(params))
	parsedParams = tal.classifyDict(params, classifiers.serverClientClassifier)
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
	def fitnessFileName(initPopType):
		return 'IP' + initPopType + '_fitness'
	def columnExtractor(gp):
		outFile = fitnessFileName(gp['initialPopulationType'])
		subprocess.call('cut -d \' \' -f 2 bestIndividual*.log | tail -n +4 | tr \'\n\' \' \' >> ../results/' + outFile, shell=True)
		subprocess.call('echo >> ../results/' + outFile, shell=True)
	experiment.executeAtEveryGridPointDir(columnExtractor)
	os.chdir('results')
	xlabel = 'Generations'
	ylimit = None
	margins = 0.5
	xlimit = evsAdditionalParams['genStopAfter']
	alpha = 0.3
	yscale = 'log'

	def plotTSForTheTwoInitalPopulationTypes():
		title = 'Fitness time series for the two types of initial populations'
		dataDict = {x: -1.*np.loadtxt(fitnessFileName(x)) for x in ['random', 'sparse']}

		# Plotting logscale average and trajectory scatter
		xscale = 'log'
		tplt.plotAverageTimeSeries(dataDict, 'Error', 'errorComparisonLog.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins)
		tplt.plotAllTimeSeries(dataDict, 'Error', 'errorAllTrajectoriesLog.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha)

		# Plotting linscale average and trajectory scatter
		xscale = 'lin'
		tplt.plotAverageTimeSeries(dataDict, 'Error', 'errorComparisonLin.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins)
		tplt.plotAllTimeSeries(dataDict, 'Error', 'errorAllTrajectoriesLin.png', title=title, legendLocation=1, xlabel=xlabel, xlimit=xlimit, ylimit=ylimit, xscale=xscale, yscale=yscale, margins=margins, alpha=alpha)

	plotTSForTheTwoInitalPopulationTypes()
	os.chdir('..')
