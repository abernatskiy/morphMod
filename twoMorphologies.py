from os.path import join, exists
from time import sleep
import subprocess

import pbsGridWalker.grid as gr
import pbsGridWalker.routes as rt
import pbsGridWalker.tools.algorithms as tal
import pbsGridWalker.tools.iniWriter as tiniw
import pbsGridWalker.tools.fsutils as tfs

# constant auxiliary definitions
randSeedFile = join(rt.pbsGridWalker, 'seedFiles', 'randints1416551751.dat')
serverClientClassifier = {'server': ['randSeed'], 'client': ['sensorAttachment']}
evsClassifier = {'classes': ['individual', 'communicator', 'evolver'],
             'indivParams': ['length', 'mutProbability', 'mutInsDelRatio', 'mutExploration',
                             'initLowerLimit', 'initUpperLimit', 'lowerCap', 'upperCap',
                             'initProbabilityOfConnection', 'mutationAmplitude'],
              'evolParams': ['populationSize', 'genStopAfter', 'initialPopulationType',
                             'secondObjectiveProbability']
             }
arrowbotsClassifier = {'arrowbot parameters': ['segments', 'sensorAttachmentType'],
                    'simulation parameters': ['simulationTime', 'timeStep', 'integrateError', 'writeTrajectories']
                   }
evsExecutable = join(rt.home, 'morphMod', 'evs', 'evsServer.py')
arrowbotsExecutable = join(rt.home, 'morphmod', 'arrowbots', 'arrowbotEvaluator')

# definition of hyperparameters
evsAdditionalParams = {'individual': 'integerWeightsSwitchableConnections', 'communicator': 'unixPipe', 'evolver': 'cluneSimplified',
                       'length': 3, 'initLowerLimit': -1, 'initUpperLimit': 1, 'lowerCap': -1, 'upperCap': 1, 'mutExploration': 0.5, 'mutInsDelRatio': 1, 'mutationAmplitude': 1,
                       'populationSize': 30, 'initialPopulationType': 'sparse', 'genStopAfter': 10, 'secondObjectiveProbability': 1, 'logParetoFront': 'yes', 'logBestIndividual': 'yes'
                      }
arrowbotsAdditionalParams = {'segments': 3,
                             'simulationTime': 10., 'timeStep': 0.1, 'integrateError': 'no', 'writeTrajectories': 'yes'
                            }

# definitions required for pbsGridWalker
computationName = 'twoMorphologies'
parametricGrid = gr.Grid1d('sensorAttachment', ['identity', 'null'])*gr.Grid1dFromFile('randomSeed', randSeedFile, size=50)

def prepareEnvironment(experiment):
	if not exists(arrowbotsExecutable):
		raise RuntimeError('Arrowbots executable not found at ' + arrowbotsExecutable)
	if not exists(evsExecutable):
		raise RuntimeError('EVS executable not found at ' + evsExecutable)

def processResults(experiment):
	pass

def runComputationAtPoint(worker, params):
	serverParams, clientParams = tal.classifyDict(params, serverClientClassifier)
	serverParams = tal.sumOfDicts(serverParams, evsAdditionalParams)
	clientParams = tal.sumOfDicts(clientParams, arrowbotsAdditionalParams)
	tiniw.write(serverParams, evsClassifier, 'evs.ini')
	tiniw.write(clientParams, arrowbotsClassifier, 'arrowbot.ini')

	geneFifo = tfs.makeUniqueFifo('.', 'genes')
	evalFifo = tfs.makeUniqueFifo('.', 'evals')

	clientProc = worker.spawnProcess([arrowbotsExecutable, geneFifo, evalFifo])
	if not worker.runCommand([evsExecutable, evalFifo, geneFifo, str(serverParams['randomSeed']), 'evs.ini']):
		return False
	worker.killProcess(clientProc)

	# TODO: Validation of the obtained files here

	return True

# auxiliary definitions for pbsGridWalker
#pointsPerJob = 1
#passes = 1
queue = 'shortq'
#maxJobs = 1
expectedWallClockTime = '00:05:00'
involvedGitRepositories = {'evs': join(rt.home, 'morphMod', 'evs'), 'arrowbots': join(rt.home, 'morphMod', 'arrowbots')}
#dryRun = False
