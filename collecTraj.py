
from func.eq_helpers.FI import *
import sys
from func.helpers import *
from func.network import *
from func.model_setup_experimental import *
from func.ABChelpers import *
import matplotlib
# matplotlib.use('PDF')
import matplotlib.pylab as plt
from matplotlib import rc
plt.rcParams['ps.useafm'] = True
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
from functools import partial
from multiprocessing import Pool
from func.eq_helpers.MFcomp import *

params= {'J': 1.4,#10.0,
 'g': 4.0,
 'N': 10000,#40000,
 'epsilon': 0.2,
 'eta': 0.0,
 'p_rate': 950.,#2000.,#2884.,#500.,#np.random.uniform(2070,2080,1)[0],#2077.4792278063533,
 'J_ext':1.0,#.5,#3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'constantI':0.0,
 'Ks': [80,20],
 'V_m': 0.0,
 'b': 0.05,#0.05
 'a': 0.0,
 'tau_w':3000.,#5000 - original from 26.04.20. #10000.0,#17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 100,
 'voltage':True,
 'chunk': False,
 'chunk_size': 1000.0,
 'directory': 'sim/',
 'simulation': 'hash',
 'simtime':500500.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 20} #40
#

params= {'J': 1.4,#10.0,
 'g': 4.0,
 'N': 10000,#40000,
 'epsilon': 0.2,
 'eta': 0.0,
 'p_rate': 850.,#2000.,#2884.,#500.,#np.random.uniform(2070,2080,1)[0],#2077.4792278063533,
 'J_ext':1.0,#.5,#3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'constantI':0.0,
 'Ks': [80,20],
 'V_m': 0.0,
 'b': 0.05,#0.006,#1.0,
 'a': 0.0,
 'tau_w':5000.,#10000.0,#17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 100,
 'voltage':True,
 'chunk': False,
 'chunk_size': 1000.0,
 'directory': 'sim/',
 'simulation': 'hash',
 'simtime':500500.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 40}


def getTrajectories(params):
    """ 
    Sample burst trajectories from random seeds 
    """
    seeds = np.arange(1000,1010,1)
    #seeds = [1000]
    #seeds = [1000]
    Res = []
    for seed in seeds:
        params['master_seed']=int(seed)
    #    time,Vs, Ws = collectVoltage(params, force = False)
    #    meanW = np.mean(Ws,0)
     #   stdW = np.std(Ws,0)
        name = get_hash(params)
        print(name)
        #st,gid = read_gdf(params['directory'],name,(0,params['simtime']),threads = params['threads'])
        #meanW_burst,meanSC_burst = meanTraj(meanW,st,gid,params,bin_size =10, primer=(20,50),
        #                            interp =False,smooth=False,smooth_par=(9,5), interp_dt=0.05)
        #Res.append([meanW_burst,meanSC_burst])
    return Res

ResBic = getTrajectories(params)
np.save('ResTraj_b0.01t8000',ResBic)

