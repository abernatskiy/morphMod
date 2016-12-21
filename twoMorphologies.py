from os.path import join, exists, expanduser
from time import sleep
import subprocess
import sys
sys.path.append(join(expanduser('~'), 'morphMod'))

import pbsGridWalker.grid as gr
import pbsGridWalker.routes as rt
import pbsGridWalker.tools.algorithms as tal
import pbsGridWalker.tools.iniWriter as tiniw
import pbsGridWalker.tools.fsutils as tfs

# constant auxiliary definitions
randSeedFile = join(rt.pbsGridWalker, 'seedFiles', 'randints1416551751.dat')
serverClientClassifier = {'server': ['randomSeed'], 'client': ['sensorAttachmentType']}
evsClassifier = {'classes': ['individual', 'communicator', 'evolver'],
             'indivParams': ['length', 'mutProbability', 'mutInsDelRatio', 'mutExploration',
                             'initLowerLimit', 'initUpperLimit', 'lowerCap', 'upperCap',
                             'initProbabilityOfConnection', 'mutationAmplitude'],
              'evolParams': ['populationSize', 'genStopAfter', 'initialPopulationType', 'trackAncestry',
                             'secondObjectiveProbability', 'logPopulation', 'logPopulationPeriod',
                             'logBestIndividual', 'logBestIndividualPeriod', 'printBestIndividual',
                             'printBestIndividualPeriod', 'printParetoFront', 'printParetoFrontPeriod',
                             'printPopulation', 'printPopulationPeriod', 'printGeneration',
                             'printGenerationPeriod', 'backup', 'backupPeriod', 'logParetoFront',
                             'logParetoFrontPeriod', 'logParetoFrontKeepAllGenerations',
                             'randomSeed']
             }
arrowbotsClassifier = {'arrowbot parameters': ['segments', 'sensorAttachmentType'],
                    'simulation parameters': ['simulationTime', 'timeStep', 'integrateError', 'writeTrajectories']
                   }
evsExecutable = join(rt.home, 'morphMod', 'evs', 'evsServer.py')
arrowbotsExecutable = join(rt.home, 'morphMod', 'arrowbots', 'arrowbotEvaluator')

def writeSSV(listOfLists, filename):
	with open(filename, 'w') as file:
		for list in listOfLists:
			file.write(' '.join(map(str, list)) + '\n')

# definition of hyperparameters
evsAdditionalParams = {'individual': 'integerWeightsSwitchableConnections', 'communicator': 'unixPipe', 'evolver': 'cluneSimplified',
                       'length': 18, 'initLowerLimit': -1, 'initUpperLimit': 1, 'lowerCap': -1, 'upperCap': 1, 'mutExploration': 0.5, 'mutInsDelRatio': 1, 'mutationAmplitude': 1,
                       'populationSize': 30, 'initialPopulationType': 'sparse', 'genStopAfter': 1000, 'secondObjectiveProbability': 1, 'logParetoFront': 'yes', 'logBestIndividual': 'yes',
                         'logParetoFrontKeepAllGenerations': 'yes', 'logParetoFrontPeriod': 100
                      } # length is 2*segments^2!!!
arrowbotsAdditionalParams = {'segments': 3,
                             'simulationTime': 10., 'timeStep': 0.1, 'integrateError': 'no', 'writeTrajectories': 'no'
                            }
numTrials = 50
arrowbotInitialConditions = [[1-0.57735,-1+2*0.57735, -2*0.57735], [-0.57735,1+2*0.57735,-2*0.57735], [-0.57735,2*0.57735,1-2*0.57735]]
arrowbotTargetOrientations = [[1,0,0], [0,1,0], [0,0,1]]

# definitions required for pbsGridWalker
computationName = 'twoMorphologiesInitSqrt3'
attTypes = ['identity', 'null'] # auxiliary
parametricGrid = gr.Grid1d('sensorAttachmentType', attTypes)*gr.Grid1dFromFile('randomSeed', randSeedFile, size=numTrials)

def prepareEnvironment(experiment):
	if not exists(arrowbotsExecutable):
		raise RuntimeError('Arrowbots executable not found at ' + arrowbotsExecutable)
	if not exists(evsExecutable):
		raise RuntimeError('EVS executable not found at ' + evsExecutable)

def processResults(experiment):
	import os
	import shutil
	import numpy as np
	import pbsGridWalker.tools.plotutils as tplt
	tfs.makeDirCarefully('results', maxBackups=100)
	def columnExtractor(gp):
		outFile = gp['sensorAttachmentType'] + 'Fitness'
		subprocess.call('cut -d \' \' -f 2 bestIndividual*.log | tail -n +4 | tr \'\n\' \' \' >> ../results/' + outFile, shell=True)
		subprocess.call('echo >> ../results/' + outFile, shell=True)
	experiment.executeAtEveryGridPointDir(columnExtractor)
	os.chdir('results')
	tplt.plotAverageTimeSeries({x: np.loadtxt(x + 'Fitness') for x in attTypes}, '-Error', 'fitnessComparison50.png', title='Fitness time series for the two types of sensors attachment', xlimit=50)
	tplt.plotAllTimeSeries({x: np.loadtxt(x + 'Fitness') for x in attTypes}, '-Error', 'fitnessAllTrajectories50.png', title='Fitness time series for the two types of sensors attachment', xlimit=50)
	tplt.plotAverageTimeSeries({x: np.loadtxt(x + 'Fitness') for x in attTypes}, '-Error', 'fitnessComparison.png', title='Fitness time series for the two types of sensors attachment')
	tplt.plotAllTimeSeries({x: np.loadtxt(x + 'Fitness') for x in attTypes}, '-Error', 'fitnessAllTrajectories.png', title='Fitness time series for the two types of sensors attachment')
	os.chdir('..')

def runComputationAtPoint(worker, params):
	print('Running evs-arrowbots pair with the following parameters: ' + str(params))
	parsedParams = tal.classifyDict(params, serverClientClassifier)
	serverParams = tal.sumOfDicts(parsedParams['server'], evsAdditionalParams)
	clientParams = tal.sumOfDicts(parsedParams['client'], arrowbotsAdditionalParams)
	tiniw.write(serverParams, evsClassifier, 'evs.ini')
	tiniw.write(clientParams, arrowbotsClassifier, 'arrowbot.ini')
	writeSSV(arrowbotInitialConditions, 'initialConditions.dat')
	writeSSV(arrowbotTargetOrientations, 'targetOrientations.dat')

	geneFifo = tfs.makeUniqueFifo('.', 'genes')
	evalFifo = tfs.makeUniqueFifo('.', 'evals')

	clientProc = worker.spawnProcess([arrowbotsExecutable, geneFifo, evalFifo])
	if not worker.runCommand([evsExecutable, evalFifo, geneFifo, str(serverParams['randomSeed']), 'evs.ini']):
		return False
	worker.killProcess(clientProc)

	# TODO: Validation of the obtained files here

	return True

# auxiliary definitions for pbsGridWalker
pointsPerJob = 20
#passes = 1
queue = 'shortq'
maxJobs = 50
expectedWallClockTime = '00:30:00'
involvedGitRepositories = {'evs': join(rt.home, 'morphMod', 'evs'), 'arrowbots': join(rt.home, 'morphMod', 'arrowbots'), 'morphMod': join(rt.home, 'morphMod')}
#dryRun = False
