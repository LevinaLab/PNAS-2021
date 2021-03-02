
# Check posteriors of ABC fitting
import numpy as np
from scipy import stats
import sys
from func_.model_setup_experimental import *
na = np.array
import pandas as pd
na = np.array
import seaborn as sns
import pandas as pd

#Generate a color palette
def return_colours(balance,in_degs, plot = False):
    """ 
    get colour palette for BIC plot
    Args:
            balance (int): number of degrees to center around
            in_degs (array): array of in-degrees sampled

    Returns:
            colors(list): colors
            color_mapping(dict): degrees - color

    """

    colors_example =sns.diverging_palette(17, 230,sep=1,s = 99.9, l=60, n=21, center="dark",as_cmap=0)
    real_deg =in_degs#np.unique(np.floor(Post[eps_i][-1][0][2,:]))
    unique_deg = -1*(real_deg-int(balance))
    indis = np.arange(0,len(unique_deg))
    # unique_deg = np.arange(unique_deg.min(),unique_deg.max())
    colors = np.zeros([len(unique_deg)+1,3])
    i_n = 0
    for i,ud in enumerate(unique_deg):
        if ud==0:
            print('a')
            colors[i,:]=(0,0,0)
        elif ud>0:
            colors[i,:] = colors_example[indis[i]]#int(ud+balance[eps_i])]
        else:
            n_ind = -int(ud+np.abs(unique_deg[-1])+1)
            print(n_ind)
            colors[i,:] = colors_example[n_ind]#int(ud-unique_deg[-1]-2)]
    colors[-1,:]=colors_example[-1]
    #color mapping
    color_mapping = {}
    color_mapping = {int(key):i for i,key in enumerate(real_deg)}

    if plot:
        sns.palplot(colors)
    return colors,color_mapping


## Load Experimental Data
origEps2 = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.8, 0.95]
IBIs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanIBI')
sdIBI = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','SD')
CVs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanCV')
DatamIBI =np.nanmean(IBIs,0)*1000
DatamCV =np.nanmean(CVs,0)
#Load bic data

Eps = [0.05,0.1,0.2,0.3,0.5,0.7,0.8,0.92]#[0.05,0.1,0.2,0.5,0.7,0.8,0.92]
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
    mypath = 'ABC_results/BIC20/'
    files = [mypath+i  for i in listdir(mypath) if str(eps) in i]
    results = np.load(files[0],allow_pickle =True)
    Res.append(results[0])
N =1000
NI = N*na(Eps)
NE = N-NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
N= 1000
# KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
balance = KE/4
n_samples = 20
pl_d=10
mn_d=5
n_deg = len(np.arange(int(balance[0])-mn_d,int(balance[0])+pl_d,1))
MSE = np.zeros([len(Eps),n_deg,n_samples])

plt.figure(figsize = (5,5))
for eps_i,eps in enumerate(Eps):
    plt.subplot(4,2,eps_i+1)
    all_data = Res[eps_i]
    in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
    in_degs= in_degs[in_degs>=0]
    print(balance[eps_i]-in_degs)
    data = na(mean_bic[eps_i],dtype ='float')#/(mean_bic[0][0])
    colors, colors_mapping = return_colours(balance[eps_i],in_degs)
    for deg_i,deg in enumerate(in_degs):
        for sample_i in np.arange(n_samples):
            bic_data = all_data[0,:,deg_i,sample_i]/all_data[0,0,deg_i,sample_i]
            MSE[eps_i,deg_i,sample_i]= np.nanmean((data-bic_data)**2)
            color =colors[deg_i]
            plt.plot(c,bic_data,'-', color = color,alpha =0.2)
    plt.errorbar(c,data,na(std_bic[eps_i]),fmt='--k',linewidth = 2)
    plt.xscale('symlog')
    plt.xlabel('ln([BIC])')
    plt.ylabel('realtive increase in IBI')
    plt.legend()
# plt.errorbar(c,np.mean(bic_*1000,1),np.std(na(bic_,dtype = 'float')*1000,1), fmt= 'ok')
plt.show(block=False)

MSE[~np.isfinite(MSE)] = np.nan
sns.set_context('paper')
#plt.figure(figsize=(4,10))
fig, ax = plt.subplots(len(Eps),1,figsize = (4,10))
from mpl_toolkits.axes_grid1 import make_axes_locatable
for eps_i,eps in enumerate(Eps):
    #plt.subplot(len(Eps),1,eps_i+1)
    ax2 = ax[eps_i].twinx()
    all_data = Res[eps_i]
    ax[eps_i].text(2,0.5,eps)
    in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
    in_degs= in_degs[in_degs>=0]
    MSE_matrix = na(MSE[eps_i,:,:])
    #plt.plot(np.repeat(in_degs,n_samples)-balance[eps_i],np.hstack(MSE[eps_i,:]),'.')
    # plt.plot(in_degs,np.nanmean(MSE_matrix,1),'-')
    if eps ==0.92:
        #add a axis break
        #see https://stackoverflow.com/questions/44731152/matplotlib-create-broken-axis-in-subplot
        divider = make_axes_locatable(ax[eps_i])
        ax3 = divider.new_vertical(size="20%", pad=0.1)
        fig.add_axes(ax3)
        ax3.errorbar(in_degs-balance[eps_i],np.nanmean(MSE_matrix,1),np.nanstd(MSE_matrix,1)/np.sqrt(20),fmt='r-')
        ax3.tick_params(bottom="off", labelbottom='off')
        ax3.get_xaxis().set_visible(False)
        ax3.spines['bottom'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.set_yticks([100])
        ax3.set_ylim([80,120])
        ax3.set_xlim([-10,10])

        d = .015  # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax3.transAxes, color='k', clip_on=False)
        ax3.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax3.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

        kwargs.update(transform=ax[eps_i].transAxes)  # switch to the bottom axes
        ax[eps_i].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax[eps_i].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    ax[eps_i].errorbar(in_degs-balance[eps_i],np.nanmean(MSE_matrix,1),np.nanstd(MSE_matrix,1)/np.sqrt(20),fmt='r-')
    ax[eps_i].axvline(0,color ='k')
    all_kde = stats.gaussian_kde(X[eps_i][1])
    x = np.arange(in_degs[0]-10,in_degs[-1]+10,1)
    y_ = all_kde.evaluate(x)
    ax2.plot(x-balance[eps_i],y_,'o',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.xlabel('$K^I$')
    #plt.plot(x-balance[eps_i],np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),'o',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.fill_between(x-balance[eps_i],0,np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),color= 'gray',alpha =0.2)
    ax2.fill_between(x-balance[eps_i],0,y_,color= 'gray',alpha =0.2)

    ax2.set_xlim([-10,10])
    ax[eps_i].set_xlim([-10,10])
    ax[eps_i].set_ylim([0,8])
    ax[eps_i].set_ylabel('MSE')
    ax2.set_ylabel('$P(K^I)$',color = 'gray')
    ax[eps_i].spines['top'].set_visible(False)
    ax[eps_i].spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_edgecolor('gray')
    plt.tight_layout()
    #sns.despine()
plt.show(block=False)

from mpl_toolkits.mplot3d import Axes3D

eps_i=4
in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
colors, colors_mapping = return_colours(balance[eps_i],in_degs)
BIC3D(mean_bic[eps_i],std_bic[eps_i],in_degs,na(MSE[eps_i]),Res[eps_i],colors_mapping,title=Eps[eps_i])
plt.show()

def BIC3D(mean_bic,std_bic,unique_deg,MSE_matrix,all_data,color_mapping,title='0'):
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

    data = na(mean_bic[0],dtype ='float')#/(mean_bic[0][0])
    for deg_i,deg in enumerate(in_degs):
        for sample_i in np.arange(20):
            #bic_data = all_data[0,:,deg_i,sample_i]/all_data[0,sample_i,deg_i,0]
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
