
# Visualize BIC grid results
import numpy as np
from scipy import stats
import sys
from func_.model_setup_experimental import *
import func_.helpers as helpers
na = np.array
import pandas as pd
na = np.array
import seaborn as sns
import pandas as pd

#Generate a color palette
## Load Experimental Data
origEps2 = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.8, 0.95]
IBIs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanIBI')
sdIBI = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','SD')
CVs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanCV')
DatamIBI =np.nanmean(IBIs,0)*1000
DatamCV =np.nanmean(CVs,0)
#Load bic data

Eps = [0.05,0.2,0.3,0.5,0.8]#[0.05,0.1,0.2,0.3,0.5,0.7,0.8,0.92]#[0.05,0.1,0.2,0.5,0.7,0.8,0.92]
#data_index = 0
#eps_i = 0 #when only one simulation
meanBic = []
mean_bic=[ ]
std_bic = []
for i,epsilon in enumerate(Eps):
    if epsilon==0.92:
        epsilon =0.95
    bic_ = na(pd.read_excel('data/bic_IBI_raw_corr.xlsx',str(epsilon)))[:7,1:]
    data_len = bic_.shape[1]
    bic_dat = bic_/bic_[0,:]
    mean_bic.append(np.nanmean(bic_dat,1))
    std_bic.append(np.nanstd(na(bic_dat,dtype = 'float'),1))
popMean_bic = np.nanmean(mean_bic,0)
# set BIC concentrations
c = [0,0.5,1,1.5,3.,10.,40.]
Kd =3
bic_conc =(1/(1+(na(c)/Kd)))
G_mod = np.round(bic_conc,2) # decrease in synaptic strength
# Get the model
# Load Posterioes
Post = []
for epsilon in Eps:
    name =  'ABC_results/nuki_fit_%s.npy'%epsilon
    posterior_= np.load(name, allow_pickle = True)
    posterior_ = [posterior_[na([p[4] for p in posterior_])>0]]
    Post.append(posterior_[0])
# Get posteriors 
X=[]
for post in Post:
    i =-2
    nu = post[i][0][0]
    ki= post[i][0][1]
    x1 = pd.Series(nu, name="$nu^{ext}$")
    x2 = pd.Series(ki, name="$K^I$")
    X.append([x1,x2])

# Load the results 
Res = []
for eps in Eps:
    mypath = 'BICgrid_results/'
    files = [mypath+i  for i in listdir(mypath) if str(eps) in i]
    results = np.load(files[0],allow_pickle =True)
    Res.append(results[0])

N =1000
NI = np.floor(N*na(Eps))
NE = N-NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
balance = KE/4
#Set degrees for sampling
#deg_inh = np.arange(int(balance)-5,int(balance)+10,1)

deg_inh = np.arange(int(balance[0])-4,int(balance[0])+6,1)
deg_inh= deg_inh[deg_inh>=0]
nu_grid = np.arange(450.,700.,20.)
colors, colors_mapping = helpers.return_colours(balance[0],deg_inh)

# Plot the IBI on the grid
eps_i = 0
data_index=0
plt.figure(figsize = (10,5))
for i,KI in enumerate(deg_inh):
    plt.subplot(121)
    plt.plot(nu_grid,(Res[eps_i][0,0,i,:]/1000),'o',color = colors[i],label =
             int(KI-balance[eps_i]))
    plt.yscale('symlog')
    plt.subplot(122)
    plt.plot(nu_grid,(Res[eps_i][1,0,i,:]),'-',color = colors[i])
plt.subplot(121)
sns.despine()
plt.axhline(DatamIBI[data_index]/1000,color = 'k')
x0=np.min(nu_grid)
x1 =np.max(nu_grid)
plt.fill_between([x0,x1],np.nanmin(na(IBIs)[:,data_index]),np.nanmax(na(IBIs)[:,data_index]),color=
                 'gray',alpha =.2)
plt.ylabel('IBI(s)')
plt.xlabel('$\\nu_{ext}$')
plt.legend(title= '$K^I-K^E/g$')
sns.despine()
plt.subplot(122)
plt.fill_between([x0,x1],np.nanmin(na(CVs)[:,data_index]),np.nanmax(na(CVs)[:,data_index]),color='gray',alpha=.2)
plt.axhline(DatamCV[data_index],color = 'k')
plt.ylabel('CV')
plt.xlabel('$\\nu_{ext}$')
sns.despine()
plt.tight_layout()
plt.show(block = False)

n_samples = 1
MSE = np.zeros([len(Eps),len(deg_inh),len(nu_grid)])
for eps_i,eps in enumerate(Eps):
    #plt.figure(figsize = (10,5))
    data = na(mean_bic[eps_i],dtype ='float')#/(mean_bic[0][0])
    all_data = Res[eps_i]
    for deg_i,KI in enumerate(deg_inh):
        for nu_i,nu in enumerate(nu_grid):
            #plt.plot(nu_grid,(DatamIBI[2]- Res[0][0,0,i,:])**2,'o',color = colors[i],label =
            #         KI-balance)
            bic_data = all_data[0,:,deg_i,nu_i]/all_data[0,0,deg_i,nu_i]
            if bic_data[2]==0:
                MSE[eps_i,deg_i,nu_i] = np.nan
            else:
                MSE[eps_i,deg_i,nu_i]= np.nanmean((data-bic_data)**2)
        #plt.yscale('symlog')

color_nu = sns.color_palette("bone",n_colors = len(nu_grid))
import matplotlib.colors as mcolors
import matplotlib.cm as cm
plt.figure(figsize= (10,3))
for i,eps in enumerate(Eps):
    plt.subplot(1,len(Eps),i+1)
    deg_inh = np.arange(int(balance[i])-4,int(balance[i])+6,1)
    deg_inh= deg_inh[deg_inh>=0]
    for nu_i,nu in enumerate(nu_grid):
        plt.plot(deg_inh-balance[i],MSE[i,:,nu_i],'-o',color = color_nu[nu_i],label = nu)
        plt.ylim([0,5])
        plt.xlabel('$K^I-K^E/g$')
        plt.axvline(0,color = 'k')
    if i==0:
        plt.ylabel('MSE')
    if i==len(Eps)-1:
        normalize = mcolors.Normalize(vmin=nu_grid.min(), vmax=nu_grid.max())
        colormap = cm.bone
        # setup the colorbar
        scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
        scalarmappaple.set_array(nu_grid)
        cbar = plt.colorbar(scalarmappaple)
        cbar.set_label('$\\nu_{ext}$', labelpad=3)
plt.tight_layout()
sns.despine()
plt.show(block = False)

plt.figure(figsize= (10,3))
for i,eps in enumerate(Eps):
    deg_inh = np.arange(int(balance[i])-4,int(balance[i])+6,1)
    deg_inh= deg_inh[deg_inh>=0]
    plt.subplot(1,len(Eps),i+1)
    plt.title(eps)
    plt.imshow(np.log(MSE[i,:,:]),cmap
               ='Blues_r',origin='lower',aspect='auto',extent=[np.min(nu_grid),np.max(nu_grid),-4,6],vmin=-4,vmax=4)#np.min(deg_inh-balance[i]),np.max(deg_inh-balance[i])])
    plt.axhline(0)
    plt.xlabel('$\\nu_{ext}$')

    if i==0:
        plt.ylabel('$K^I-K^E/g$')
    if i==len(Eps)-1:
        cbar = plt.colorbar()
        cbar.set_label('ln(MSE)', labelpad=3)
#plt.colorbar()
plt.tight_layout()
print(np.nanmin(MSE))
plt.show(block = False)


#Visualize the simulations
eps_i=4
deg_inh = np.arange(int(balance[eps_i])-4,int(balance[eps_i])+6,1)
deg_inh= deg_inh[deg_inh>=0]
colors, colors_mapping = return_colours(balance[eps_i],deg_inh)
BIC3Dgrid(mean_bic[eps_i],std_bic[eps_i],deg_inh,na(MSE[eps_i]),Res[eps_i],colors_mapping,title=Eps[eps_i])
plt.show()

from mpl_toolkits.mplot3d import Axes3D
def BIC3Dgrid(mean_bic,std_bic,unique_deg,MSE_matrix,all_data,color_mapping,title='0'):
    """ 
    foo 
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """

    mean_lines = {}#[[0]]*(len(unique_deg)+1)
    divisor = np.zeros([len(unique_deg)+1,7])
    stds = np.zeros([len(unique_deg)+1,7])
    data = na(mean_bic[0])#/mean_bic[0][0]
    ii=0
    c[0]=0.1

    in_degs = unique_deg
    data = na(mean_bic[0],dtype ='float')#/(mean_bic[0][0])
    for deg_i,deg in enumerate(in_degs):
        for sample_i in np.arange(13):
            #bic_data = all_data[0,:,deg_i,sample_i]/all_data[0,sample_i,deg_i,0]
            if all_data[0,1,deg_i,sample_i]>0:
                bic_data = all_data[0,:,deg_i,sample_i]/all_data[0,0,deg_i,sample_i]
                color =colors[deg_i]
                try:
                    mean_lines[deg].append(bic_data/bic_data[0])
                except:
                    mean_lines[deg] = [bic_data/bic_data[0]]
    
    # degrees_vector
    # min_ = np.nanargmin(MSE)
    min_ = np.nanargmin(np.nanmean(MSE_matrix,1))#np.nanargmin(MSE)

    best_deg = in_degs[min_]#
    # best_deg = balance[0]#degrees_vector[min_]
    plt.figure(figsize=(10,6))
    ax = plt.subplot(projection='3d')
    ax.set_title(title)
    xs = mean_bic
    ys = np.ones(len(bic_data))*best_deg
    zs = std_bic/np.sqrt(4)
    ax.plot(np.log(c),ys,xs,'--o',color = 'C1',alpha = 1.8)
    for i, j in enumerate(xs):
        ax.plot([np.log(c)[i], np.log(c)[i]], [ys[i], ys[i]], [j+zs[i], j-zs[i]], marker="_", color='C1')

    for deg in mean_lines.keys():
    #     if ~np.isnan(mean_lines[deg][0][-1]):
        print(mean_lines[deg][0][-1])
    #     print(mean_lines[deg]/mean_lines[deg][0])
        mean_data = np.nanmean(na(mean_lines[deg]),0)
        std_data = np.nanstd(na(mean_lines[deg]),0)/np.sqrt(len(mean_lines[deg]))
        x = np.ones(len(bic_data))*deg
        col_ind = color_mapping[deg]
        ax.plot(np.log(c),x,mean_data,'o-', color = colors[col_ind],alpha =0.8)    
        ax.add_collection3d(plt.fill_between(np.log(c), mean_data-std_data,mean_data+std_data, color=colors[col_ind], alpha=0.3), zs=deg, zdir='y')
    ax.view_init(30,220)

    ax.set_xlabel('ln [BIC]')
    ax.set_ylabel('$K^I$')
    # a = ax.zaxis.label.get_rotation() # put the actual angle in the z-label
    # ax.set_zlabel(u'z-rot = {:.1f}°'.format(a))
    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel('normalized IBI',rotation=90)


