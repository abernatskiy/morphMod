serverClientClassifier = {'server': ['randomSeed', 'initialPopulationType'], 'client': ['sensorAttachmentType']}

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
                             'randomSeed']}

arrowbotsClassifier = {'arrowbot parameters': ['segments', 'sensorAttachmentType'],
                    'simulation parameters': ['simulationTime', 'timeStep', 'integrateError', 'writeTrajectories']}
