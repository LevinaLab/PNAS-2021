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

Eps = [0.05,0.1,0.2,0.3,0.5,0.7,0.8,0.92]#,0.1,0.2,0.3,0.5,0.7,0.8,0.92]#[0.05,0.1,0.2,0.5,0.7,0.8,0.92]
meanBic = []
mean_bic=[ ]
std_bic = []
AllBIC = []
for i,epsilon in enumerate(Eps):
    if epsilon==0.92:
        epsilon =0.95
    bic_ = na(pd.read_excel('data/bic_IBI_raw_corr_full.xlsx',str(epsilon)))[:7,1:]
    data_len = bic_.shape[1]
    bic_dat = bic_/bic_[0,:]
    AllBIC.append(bic_dat)
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
    name =  'ABC_results/nuki_fit_LastStep%s.npy'%epsilon
    posterior_= np.load(name, allow_pickle = True)
    posterior_ = [posterior_[na([p[4] for p in posterior_])>0]]
    Post.append(posterior_[0])
# Get posteriors 
X=[]
for post in Post:
    i =-1
    nu = post[i][0][0]
    ki= post[i][0][1]
    x1 = pd.Series(nu, name="$nu^{ext}$")
    x2 = pd.Series(ki, name="$K^I$")
    X.append([x1,x2])

# Load the results 
Res = []
for eps in Eps:
    mypath = 'ABC_results/BICLastStep/'
    files = [mypath+i  for i in listdir(mypath) if str(eps) in i]
    results = np.load(files[0],allow_pickle =True)
    Res.append(results[0])

color_mapping_data = {}
col_data =na(pd.read_csv("nuki/Data_palette2.csv"))/255
col_model =na(pd.read_csv("nuki/Model_palette.csv"))/255
# col_data = sns.color_palette('Blues',n_colors=len(origEps2)*2)[len(origEps2)-4:]
# col_model = sns.color_palette('OrRd',n_colors=len(Eps)*2)[len(Eps):]

color_mapping_model = {}
color_mapping_data = {}
for i,eps in enumerate(origEps2):
    color_mapping_data[eps]=col_data[i]
    
for i,eps in enumerate(Eps):
    color_mapping_model[eps]=col_model[i] 

N =1000
NI = N*na(Eps)
NE = N-NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
N= 1000
# KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
balance = KE/4
n_samples = 10
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
    mse_norm = np.nanmean((data-np.ones(len(data)))**2)
    colors, colors_mapping = return_colours(balance[eps_i],in_degs)
    for deg_i,deg in enumerate(in_degs):
        for sample_i in np.arange(n_samples):
            bic_data = all_data[0,:,deg_i,sample_i]/all_data[0,0,deg_i,sample_i]
            MSE[eps_i,deg_i,sample_i]= np.nanmean((data-bic_data)**2)/mse_norm
            color =colors[deg_i]
            plt.plot(c,bic_data,'-', color = color,alpha =0.2)
    plt.errorbar(c,data,na(std_bic[eps_i]),fmt='--k',linewidth = 2)
    plt.xscale('symlog')
    plt.xlabel('ln([BIC])')
    plt.ylabel('realtive increase in IBI')
    plt.legend()
# plt.errorbar(c,np.mean(bic_*1000,1),np.std(na(bic_,dtype = 'float')*1000,1), fmt= 'ok')
plt.show(block=False)
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

MSE[~np.isfinite(MSE)] = np.nan
sns.set_context('paper')
#plt.figure(figsize=(4,10))
fig, ax = plt.subplots(len(Eps),1,figsize = (4,1.5*len(Eps)))
from mpl_toolkits.axes_grid1 import make_axes_locatable
for eps_i,eps in enumerate(Eps):
    #plt.subplot(len(Eps),1,eps_i+1)
    ax2 = ax[eps_i].twinx()
    all_data = Res[eps_i]
    in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
    in_degs= in_degs[in_degs>=0]
    MSE_matrix = na(MSE[eps_i,:,:])
    #plt.plot(np.repeat(in_degs,n_samples)-balance[eps_i],np.hstack(MSE[eps_i,:]),'.')
    # plt.plot(in_degs,np.nanmean(MSE_matrix,1),'-')
    if eps ==False:
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

    ax[eps_i].errorbar(in_degs-balance[eps_i],np.nanmean(MSE_matrix,1),np.nanstd(MSE_matrix,1),fmt='-o',color = color_mapping_model[eps])
    ax[eps_i].axvline(0,color ='k')
    all_kde = stats.gaussian_kde(X[eps_i][1])
    x = np.arange(in_degs[0]-10,in_degs[-1]+10,1)
    y_ = all_kde.evaluate(x)
    ax2.plot(x-balance[eps_i],y_,'-',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.xlabel('$K^I$')
    #plt.plot(x-balance[eps_i],np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),'o',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.fill_between(x-balance[eps_i],0,np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),color= 'gray',alpha =0.2)
    ax2.fill_between(x-balance[eps_i],0,y_,color= 'gray',alpha =0.2)
    ax[eps_i].text(-8,0.2,'%s%%'%np.floor(eps*100))
    ax2.set_xlim([-10,10])
    ax[eps_i].set_xlim([-10,10])
    #ax[eps_i].set_ylim([0,8])
    ax[eps_i].set_ylim([0,1.8])
    ax[eps_i].axhline(1,label = 'flat response MSE',color ='gray',linestyle = '--')
    ax[eps_i].set_ylabel('nMSE')
    ax2.set_ylabel('$P(K^I)$',color = 'gray')
    ax[eps_i].spines['top'].set_visible(False)
    ax[eps_i].spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_edgecolor('gray')
    ax2.set_yticks([])
    if eps_i==len(Eps)-1:
        ax[eps_i].set_xlabel('$K^I(\\epsilon) - K^E(\\epsilon)/g$')
    plt.tight_layout()
    #sns.despine()
plt.show(block=False)

# Cheche differences 
from scipy import stats
for eps_i,eps in enumerate(Eps):
    #plt.subplot(len(Eps),1,eps_i+1)
    all_data = Res[eps_i]
    in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
    in_degs= in_degs[in_degs>=0]
    MSE_matrix = na(MSE[eps_i,:,:])
    mse_mean  = np.nanmean(MSE_matrix,1)
    ind = np.nanargmin(mse_mean)
    p = [stats.ttest_rel(MSE_matrix[ind],MSE_matrix[ind+i]) for i in [-1,1]]
    print(eps,p[0][1],p[1][1])


from func_.report_helpers import BIC3D as BIC3D
# Plot in 3D 
eps_i = 2
colors, colors_mapping = return_colours(balance[eps_i],in_degs)
in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
in_degs= in_degs[in_degs>=0]
BIC3D(mean_bic[eps_i],std_bic[eps_i],in_degs,na(MSE[eps_i]),Res[eps_i],balance[eps_i],c,title=Eps[eps_i])
plt.savefig('figs/Figure7_supp.pdf',fmt = 'pdf')
plt.show()

# FIGURE 6
for i,epsilon in enumerate([0.25]):
    bic_ = na(pd.read_excel('data/bic_IBI_raw_corr.xlsx',str(epsilon)))[:7,1:]
    bic_dat = bic_/bic_[0,:]
    control_data =np.nanmean(bic_dat,1)
    control_std_bic = np.nanstd(na(bic_dat,dtype = 'float'),1)


# --------Panel A----------
# Prepare the data 
sns.set_context('talk')
n_samples = 10
fig = plt.figure(figsize=(12,6))
# grid = plt.GridSpec(1, 2, wspace=0.2, hspace=0.1)
StaticBIC = na([8.9066666666666663, 2.6400000000000001, 1.46875])
StaticBIC/=StaticBIC[0]
# plt.subplot(grid[0, 0])
plt.subplot(1,2,1)
indices= [2,3,4,5,6,7,8]
for i,eps in enumerate(Eps[2:-1]):
    if eps>0.:
#         print(i)
        i = indices[i]
        #Plot Data
        indi_c = i#control_indis[i]
        data = na(mean_bic[indi_c],dtype ='float')#/(mean_bic[0][0])
        dat_samples = AllBIC[indi_c].shape[1]
        sts =na(std_bic[indi_c],dtype ='float')/np.sqrt(dat_samples)#/(mean_bic[0][0])
        color = color_mapping_data[eps]
        print(origEps2[indi_c])
        plt.errorbar(c,data,sts,color= color,
                label='%s%%'%(np.int(Eps[indi_c]*100)),fmt=  'o',alpha =1.0,capsize=8,markersize = 8,zorder=i)
        
        #Plot Model
        # Get the right index 
        all_data= Res[i]

        in_degs = np.arange(int(balance[i])-mn_d,int(balance[i])+pl_d,1)
        in_degs= in_degs[in_degs>=0]
        MSE_matrix = na(MSE[i,:,:])
        MSE_matrix[MSE_matrix==np.inf] = np.nan
        MSE_matrix[MSE_matrix==0.] = np.nan
        cond=np.nanmin(MSE_matrix)
        best_ind = np.where(MSE_matrix==cond)
#         best_ind = (np.where(in_degs==int(balance[i]))[0],2)
        
        bic_data = all_data[0,:,best_ind[0],best_ind[1]][0]
        bic_data /=bic_data[0] 
        color = color_mapping_model[eps]
        if False:#eps==0.2:
            linewidth = 4
            markersize = 10
            zorder = 88
        else:
            linewidth = 2
            markersize = 8
            zorder = i+1
        plt.plot(c,bic_data,'-',color= color, #/mIBI[i,:][0]
            label='%s%%'%(np.int(Eps[i]*100)),alpha =1.0,markersize = markersize,linewidth = linewidth,zorder=zorder)
plt.xscale('symlog')
color = color_mapping_data[0.25]
plt.errorbar(c,control_data,control_std_bic,color= 'navy', mfc = 'white',label='control',fmt=  '--o',alpha =1.0,capsize=8,linewidth = 4,markersize = 8,zorder=99)


plt.ylim([-0.5,8.5])
plt.plot(c[0:len(StaticBIC)],StaticBIC,'-o',color = 'grey',markersize = 8,linewidth = 4, label ='no ad')
plt.plot(c[len(StaticBIC)-1],StaticBIC[-1],'*',color = 'grey',markersize = 18)
plt.xticks(c,c,fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(ncol=2,fontsize = 14,fancybox=False,frameon=False)
#plt.legend(bbox_to_anchor=(1,1))
sns.despine(trim = 1)
plt.ylabel('Change in IBI',fontsize = 18)
plt.xlabel("[BIC]",fontsize = 18)


plt.subplot(1,2,2)
indices = [0,1,7]
for i,eps in enumerate([0.05,0.1,0.92]):
    if eps>0.:
#         print(i)
        #Plot Data
        i = indices[i]
        indi_c = i#control_indis[i]
        data = na(mean_bic[indi_c],dtype ='float')#/(mean_bic[0][0])
        dat_samples = AllBIC[indi_c].shape[1]
        sts =na(std_bic[indi_c],dtype ='float')/np.sqrt(dat_samples)#/(mean_bic[0][0])
        if eps==0.92:
            eps_d = 0.95
        else:
            eps_d = eps
        color = color_mapping_data[eps_d]
        print(origEps2[indi_c])
        plt.errorbar(c,data,sts,color= color,
                label='%s%%'%(np.int(Eps[indi_c]*100)),fmt=  'o',alpha =1.0,capsize=8,markersize = 8,zorder=i)
        
        #Plot Model
        # Get the right index 
        all_data= Res[i]
        in_degs = np.arange(int(balance[i])-mn_d,int(balance[i])+pl_d,1)
        in_degs= in_degs[in_degs>=0]
        MSE_matrix = na(MSE[i,:,:])
        #MSE_matrix = np.reshape(MSE_all[i],[len(in_degs),n_samples])
        MSE_matrix[MSE_matrix==np.inf] = np.nan
        MSE_matrix[MSE_matrix==0.] = np.nan
        cond=np.nanmin(MSE_matrix)
        best_ind = np.where(MSE_matrix==cond)
#         best_ind = (np.where(in_degs==int(balance[i]))[0],2)
        
        bic_data = all_data[0,:,best_ind[0],best_ind[1]][0]
        bic_data /=bic_data[0] 
        color = color_mapping_model[eps]
        if eps==0.2:
            linewidth = 4
            markersize = 10
            zorder = 88
        else:
            linewidth = 2
            markersize = 8
            zorder = i+1
        plt.plot(c,bic_data,'-s',color= color, #/mIBI[i,:][0]
            label='%s%%'%(np.int(Eps[i]*100)),alpha =1.0,markersize = markersize,linewidth = linewidth,zorder=zorder)
        plt.xscale('symlog')
plt.ylim([-0.5,8.5])
plt.xticks(c,c,fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(ncol=2, fontsize = 14,fancybox=False,frameon=False)#bbox_to_anchor=(1,1))
sns.despine(trim = 1)
plt.tight_layout()
plt.ylabel('Change in IBI',fontsize = 18)
plt.xlabel("[BIC]",fontsize = 18)

plt.show(block =False)



## Load Static results
G_modStatic,staticBIC_results = np.load('static/BIC_results_0.5.npy',allow_pickle=True)
G_modStatic = na(G_modStatic)
Kd =3
bic_concStat =(Kd/G_modStatic) * (1-G_modStatic)

"""FIGURE 6 """

sns.set_context('paper')
# FIGURE 6 allrogether ogrid = plt.GridSpec(2, 3, wspace=0.4, hspace=0.3)
fig = plt.figure(figsize=(14, 7))
grid = plt.GridSpec(len(Eps), 8, wspace=2., hspace=0.0)
ax1 = fig.add_subplot(grid[2:, 0:3])
ax2 = fig.add_subplot(grid[2:, 3:6])
ax3 = [fig.add_subplot(grid[i, 6:]) for i in range(len(Eps))]
ymax = 10.0#set y-limi panel A and B

# --------Panel A----------
# Prepare the data 
n_samples = 10
# grid = plt.GridSpec(1, 2, wspace=0.2, hspace=0.1)
StaticBIC = na([8.9066666666666663, 2.6400000000000001, 1.46875])
StaticBIC/=StaticBIC[0]
# plt.subplot(grid[0, 0])
indices= [2,3,4,5,6,7,8]
for i,eps in enumerate(Eps[2:-1]):
    if eps>0.:
#         print(i)
        i = indices[i]
        #Plot Data
        indi_c = i#control_indis[i]
        data = na(mean_bic[indi_c],dtype ='float')#/(mean_bic[0][0])
        dat_samples = AllBIC[indi_c].shape[1]
        sts =na(std_bic[indi_c],dtype ='float')/np.sqrt(dat_samples)#/(mean_bic[0][0])
        color = color_mapping_data[eps]
        print(origEps2[indi_c])
        ax1.errorbar(c,data,sts,color= color,
                label='%s%%'%(np.int(Eps[indi_c]*100)),fmt=  'o',alpha =1.0,capsize=8,markersize = 8,zorder=i)
        
        #Plot Model
        # Get the right index 
        all_data= Res[i]

        in_degs = np.arange(int(balance[i])-mn_d,int(balance[i])+pl_d,1)
        in_degs= in_degs[in_degs>=0]
        MSE_matrix = na(MSE[i,:,:])
        MSE_matrix[MSE_matrix==np.inf] = np.nan
        MSE_matrix[MSE_matrix==0.] = np.nan
        cond=np.nanmin(MSE_matrix)
        best_ind = np.where(MSE_matrix==cond)
#         best_ind = (np.where(in_degs==int(balance[i]))[0],2)
        
        bic_data = all_data[0,:,best_ind[0],best_ind[1]][0]
        bic_data /=bic_data[0] 
        color = color_mapping_data[eps]
        if False:#eps==0.2:
            linewidth = 4
            markersize = 10
            zorder = 88
        else:
            linewidth = 2
            markersize = 8
            zorder = i+1
        ax1.plot(c,bic_data,'-',color= color, #/mIBI[i,:][0]
            label='%s%%'%(np.int(Eps[i]*100)),alpha =1.0,markersize = markersize,linewidth = linewidth,zorder=zorder)
ax1.set_xscale('symlog')
color = color_mapping_data[0.25]
ax1.errorbar(c,control_data,control_std_bic,color= 'navy', mfc = 'white',label='control',fmt=  '--o',alpha =1.0,capsize=8,linewidth = 4,markersize = 8,zorder=99)

ax1.set_ylim([-0.5,ymax])
ax1.plot(bic_concStat,staticBIC_results[0,:]/staticBIC_results[0,0],'-o',color = 'grey',label ='no ad')
#ax1.plot(c[len(StaticBIC)-1],StaticBIC[-1],'*',color = 'grey',markersize = 12)
ax1.set_xticks(c)
ax1.set_xticklabels(c)#,fontsize = 18)
#ax1.set_yticks(fontsize = 18)
ax1.legend(ncol=2,fancybox=False,frameon=False)
#plt.legend(bbox_to_anchor=(1,1))
ax1.set_ylabel('Change of IBI')
ax1.set_xlabel("$Bicuculline [\\mu M]$")
sns.despine(ax = ax1)

#-----PANEL B------------
indices = [0,1,7]
for i,eps in enumerate([0.05,0.1,0.92]):
    if eps>0.:
#         print(i)
        #Plot Data
        i = indices[i]
        indi_c = i#control_indis[i]
        data = na(mean_bic[indi_c],dtype ='float')#/(mean_bic[0][0])
        dat_samples = AllBIC[indi_c].shape[1]
        sts =na(std_bic[indi_c],dtype ='float')/np.sqrt(dat_samples)#/(mean_bic[0][0])
        if eps==0.92:
            eps_d = 0.95
        else:
            eps_d = eps
        color = color_mapping_data[eps_d]
        print(origEps2[indi_c])
        ax2.errorbar(c,data,sts,color= color,
                label='%s%%'%(np.int(Eps[indi_c]*100)),fmt=  'o',alpha =1.0,capsize=8,markersize = 8,zorder=i)
        #Plot Model
        # Get the right index 
        all_data= Res[i]
        in_degs = np.arange(int(balance[i])-mn_d,int(balance[i])+pl_d,1)
        in_degs= in_degs[in_degs>=0]
        MSE_matrix = na(MSE[i,:,:])
        #MSE_matrix = np.reshape(MSE_all[i],[len(in_degs),n_samples])
        MSE_matrix[MSE_matrix==np.inf] = np.nan
        MSE_matrix[MSE_matrix==0.] = np.nan
        cond=np.nanmin(MSE_matrix)
        best_ind = np.where(MSE_matrix==cond)
#         best_ind = (np.where(in_degs==int(balance[i]))[0],2)
        
        bic_data = all_data[0,:,best_ind[0],best_ind[1]][0]
        bic_data /=bic_data[0] 
        color = color_mapping_data[eps_d]
        if eps==0.2:
            linewidth = 4
            markersize = 10
            zorder = 88
        else:
            linewidth = 2
            markersize = 8
            zorder = i+1
        ax2.plot(c,bic_data,'-',color= color, #/mIBI[i,:][0]
            label='%s%%'%(np.int(Eps[i]*100)),alpha =1.0,markersize = markersize,linewidth = linewidth,zorder=zorder)
        ax2.set_xscale('symlog')
ax2.set_ylim([-0.5,ymax])
ax2.set_xticks(c)
ax2.set_xticklabels(c)#,fontsize = 18)
ax2.legend(ncol=2, fancybox=False,frameon=False)
sns.despine(ax=ax2)
ax2.set_ylabel('Change in IBI')
ax2.set_xlabel("$Bicuculline [\\mu M]$")

#------- Panel C-----------------
kde = True
for eps_i,eps in enumerate(Eps):
    if eps==0.92:
        eps_d =0.95
    else:
        eps_d = eps
    #plt.subplot(len(Eps),1,eps_i+1)
    if kde:
        ax4 = ax3[eps_i].twinx()
    all_data = Res[eps_i]
    in_degs = np.arange(int(balance[eps_i])-mn_d,int(balance[eps_i])+pl_d,1)
    in_degs= in_degs[in_degs>=0]
    MSE_matrix = na(MSE[eps_i,:,:])
    #plt.plot(np.repeat(in_degs,n_samples)-balance[eps_i],np.hstack(MSE[eps_i,:]),'.')
    # plt.plot(in_degs,np.nanmean(MSE_matrix,1),'-')
    if eps ==False:
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
        ax3[eps_i].plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        ax3[eps_i].plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    loc_min = np.nanargmin(np.nanmean(MSE_matrix,1))
    ax3[eps_i].plot((in_degs-balance[eps_i])[loc_min], np.nanmean(MSE_matrix,1)[loc_min]+0.2, 'v')
    ax3[eps_i].errorbar(in_degs-balance[eps_i],np.nanmean(MSE_matrix,1),np.nanstd(MSE_matrix,1),fmt='-o',color
                        = color_mapping_data[eps_d],markersize = 4)
    ax3[eps_i].axvline(0,color ='k')
    if kde:
        all_kde = stats.gaussian_kde(X[eps_i][1])
        x = np.arange(in_degs[0]-10,in_degs[-1]+10,0.5)
        y_ = all_kde.evaluate(x)
        ax4.plot(x-balance[eps_i],y_,'-',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.xlabel('$K^I$')
    #plt.plot(x-balance[eps_i],np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),'o',color = 'gray',label ='Posterior KDE ($K^I$) \n normalized to max MSE',alpha =0.3)
    #plt.fill_between(x-balance[eps_i],0,np.nanmax(MSE_matrix[np.isfinite(MSE_matrix)])*(y_/np.nanmax(y_)),color= 'gray',alpha =0.2)
        ax4.fill_between(x-balance[eps_i],0,y_,color= 'gray',alpha =0.2)
        ax4.set_xlim([-10,10])
        ax4.set_ylabel('$P(K^I)$',color = 'gray')
        ax4.spines['bottom'].set_visible(False)
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_edgecolor('gray')
        ax4.set_yticks([])

    ax3[eps_i].text(-8,0.2,'%s%%'%np.floor(eps*100))
    ax3[eps_i].set_xlim([-10,10])
    #ax[eps_i].set_ylim([0,8])
    ax3[eps_i].set_ylim([-0.2,1.2])
    #ax3[eps_i].axhline(1,label = 'flat response MSE',color ='gray',linestyle = '--')
    ax3[eps_i].axhline(0,color ='k')
    ax3[eps_i].set_ylabel('nMSE')
    ax3[eps_i].spines['top'].set_visible(False)
    ax3[eps_i].spines['right'].set_visible(False)
    ax3[eps_i].spines['bottom'].set_visible(False)
    if eps_i==len(Eps)-1:
        ax3[eps_i].set_xlabel('$K^I(\\epsilon) - K^E(\\epsilon)/g$')
        ax3[eps_i].spines['bottom'].set_visible(True)
    else:
        ax3[eps_i].set_xticks([])

    plt.tight_layout()
    #sns.despine()
plt.show(block=False)

plt.savefig('figs/Figure7v_fulldata.pdf',fmt='pdf')
