import numpy as np
phists = { 'n={}'.format(n) : { t : np.loadtxt('N{}/paretoHistogram{}'.format(n, t)) for t in ['random', 'sparse' ] } for n in [3,5,10] }
pmaxsize = { ns : { t : int(max(phists[ns][t][:,1])) for t in ['random', 'sparse' ] } for ns in ['n=3','n=5', 'n=10'] }
pavg = { ns : { t : sum(np.array(phists[ns][t][:,1])*phists[ns][t][:,0])/sum(phists[ns][t][:,0]) for t in ['random', 'sparse' ] } for ns in ['n=3','n=5', 'n=10'] }
pmut = { ns : { t : 50-pavg[ns][t] for t in ['random', 'sparse' ] } for ns in ['n=3','n=5', 'n=10'] }
