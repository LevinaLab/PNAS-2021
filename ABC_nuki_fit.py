# Fit ext_nu and K^I for the adLIF networks
import numpy as np
from scipy import stats
import sys
from func.helpers import *
from func.network import *
from func.model_setup_experimental import *
from func.ABChelpers import *
na = np.array
import pandas as pd

#Inputs
epsilon = np.float(sys.argv[1])
print(epsilon)
data_index= np.int(sys.argv[2])

N =1000
NI = np.int(N*epsilon)
NE = N-NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
params= {'J': 2.0,#10.0,
 'g': 4.0,
 'N': 1000,
 'epsilon': epsilon,
 'eta': 0.0,
 'p_rate': 640.,#1000.,#np.random.uniform(2070,2080,1)[0],#2077.4792278063533,
 'J_ext': 1.,#3.820498723458609,
 'tauMem': 20.0,
 'CMem': 1.0,
 'theta': 20.0,
 'V_res': 10.0,
 'Ks': [int(KE),int(KE)/4],
 'V_m': 0.0,
 'b': 0.05,#1.0,
 'a': 0.0,
 'tau_w':8000.,#10000.0,#17000.0,
 'p': 0.1,
 't_ref': 2.0,
 'd': 3.5,
 'N_rec': 200,
 'voltage': False,
 'chunk': True,
 'chunk_size':100100.0,
 'directory': 'sim/nuki_fit/%s/'%epsilon,
 'simulation': 'hash',
 'simtime': 500500.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 18}

# Original data
origEps2 = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.8, 0.95]
IBIs = pd.read_excel('data/meanFeatures.xlsx','meanIBI')
sdIBI = pd.read_excel('data/meanFeatures.xlsx','SD')
CVs = pd.read_excel('data/meanFeatures.xlsx','meanCV')
DatamIBI =np.nanmean(IBIs,0)*1000
DatamCV =np.nanmean(CVs,0)

nest.set_verbosity('M_FATAL')
Eps = [epsilon]
# Fit into real data 
#Set a random seed 
np.random.seed(914)
#We need some real values so we know it is working. 
#Test: 
theta_0 = [params['p_rate'], params['Ks'][1]]#[2080.,3.8,8]
model =nuki_fit(params,orig_data=[DatamIBI[data_index],DatamCV[data_index]],
        alpha=1,beta =1)

data = model.generate_data(theta_0)
print(data)
#Now we need to set the prior distributions. We were clever and set up our draw theta method to call 
#frozen scipy.stats objects, so we jest need to give the model our prior distributions 
#model.set_prior([stats.uniform(1,10000),stats.uniform(0.1,2),stats.uniform(0,80)])
#oprg priors 
#100,5000 
#5, 100
model.set_prior([stats.uniform(100,2000),stats.uniform(2,50)])
#And while we are at it, we will give the model our observed data as well
model.set_data([DatamIBI.copy()[data_index],DatamCV.copy()[data_index]])
print('first simulation is done')
#print(data)
name = 'nuki_fit_%s'%epsilon
# Save paramters to the folder with simulations
save_obj(params,params['directory']+'/params')

if np.int(sys.argv[3])==1:
    post =np.load('ABC_results/'+ name+'.npy', allow_pickle = True)
    post = post[na([p[4] for p in post])>0]
    print('continuing')
    print('Fitting eps %s , %s'%(epsilon,origEps2[data_index]))
    posterior2 = pmc_abc(model,[DatamIBI.copy()[data_index],DatamCV.copy()[data_index]] ,
                        epsilon_0=1.5, min_samples=50, steps=15,
                        file=name,
                        minError=0.05,resume = post )
else:
    posterior2 = pmc_abc(model,[DatamIBI.copy()[data_index],DatamCV.copy()[data_index]] ,
                    epsilon_0=1.5, min_samples=50, steps=15,
                    file=name,
                    minError=0.05)
