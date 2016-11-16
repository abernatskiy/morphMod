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
serverClientClassifier = {'server': ['randomSeed'], 'client': ['sensorAttachment']}
evsClassifier = {'classes': ['individual', 'communicator', 'evolver'],
             'indivParams': ['length', 'mutProbability', 'mutInsDelRatio', 'mutExploration',
                             'initLowerLimit', 'initUpperLimit', 'lowerCap', 'upperCap',
                             'initProbabilityOfConnection', 'mutationAmplitude'],
              'evolParams': ['populationSize', 'genStopAfter', 'initialPopulationType', 'trackAncestry'
                             'secondObjectiveProbability', 'logPopulation', 'logPopulationPeriod',
                             'logBestIndividual', 'logBestIndividualPeriod', 'printBestIndividual',
                             'printBestIndividualPeriod', 'printParetoFront', 'printParetoFrontPeriod',
                             'printPopulation', 'printPopulationPeriod', 'printGeneration',
                             'printGenerationPeriod', 'backup', 'backupPeriod', 'logParetoFront',
                             'logParetoFrontPeriod', 'logParetoFrontKeepAllGenerations']
             }
arrowbotsClassifier = {'arrowbot parameters': ['segments', 'sensorAttachment'],
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
arrowbotInitialConditions = [[0,0,0], [0,0,0], [0,0,0]]
arrowbotTargetOrientations = [[1,0,0], [0,1,0], [0,0,1]]

# definitions required for pbsGridWalker
computationName = 'twoMorphologies'
parametricGrid = gr.Grid1d('sensorAttachment', ['identity', 'null'])*gr.Grid1dFromFile('randomSeed', randSeedFile, size=numTrials)

def prepareEnvironment(experiment):
	if not exists(arrowbotsExecutable):
		raise RuntimeError('Arrowbots executable not found at ' + arrowbotsExecutable)
	if not exists(evsExecutable):
		raise RuntimeError('EVS executable not found at ' + evsExecutable)

def processResults(experiment):
	pass

def runComputationAtPoint(worker, params):
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
involvedGitRepositories = {'evs': join(rt.home, 'morphMod', 'evs'), 'arrowbots': join(rt.home, 'morphMod', 'arrowbots')}
#dryRun = False
