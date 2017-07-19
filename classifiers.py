serverClientClassifier = {'server': ['randomSeed', 'initialPopulationType'], 'client': ['sensorAttachmentType']}

# Since EVS has supported composite Individuals, it now requires a regexp classifier
_indivClassParams = ['length', 'mutProbability', 'mutInsDelRatio', 'mutExploration',
                     'initLowerLimit', 'initUpperLimit', 'lowerCap', 'upperCap',
                     'initProbabilityOfConnection', 'mutationAmplitude']
_indivClassParamsRegexp = '(' + '|'.join(_indivClassParams) + ')'
_classSuffixRegexp = 'Class[0-9]+'
evsClassifier = {'classes': ['individual', 'communicator', 'evolver'],
             'indivParams': [ _indivClassParamsRegexp,
                              _indivClassParamsRegexp + _classSuffixRegexp,
                              '(composite|probabilityOfMutating)' + _classSuffixRegexp ],
              'evolParams': ['populationSize', 'genStopAfter', 'initialPopulationType', 'trackAncestry',
                             'secondObjectiveProbability', 'logPopulation', 'logPopulationPeriod',
                             'logBestIndividual', 'logBestIndividualPeriod', 'printBestIndividual',
                             'printBestIndividualPeriod', 'printParetoFront', 'printParetoFrontPeriod',
                             'printPopulation', 'printPopulationPeriod', 'printGeneration',
                             'printGenerationPeriod', 'backup', 'backupPeriod', 'logParetoFront',
                             'logParetoFrontPeriod', 'logParetoFrontKeepAllGenerations',
                             'randomSeed', 'morphologyControlIndivs']}

arrowbotsClassifier = {'arrowbot parameters': ['segments', 'sensorAttachmentType'],
                    'simulation parameters': ['simulationTime', 'timeStep', 'integrateError', 'writeTrajectories']}
