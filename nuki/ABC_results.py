# Check posteriors of ABC fitting
#
import numpy as np
from scipy import stats
import sys
sys.path.insert(0,'/home/ovinogradov/Projects/EI_project/adLIF_NEST/')
sys.path.insert(0,'/home/ovinogradov/Projects/EI_project/ABC_bursting/')
from func.helpers import *
from func.network import *
from func_.model_setup_experimental import *
na = np.array
import seaborn as sns
import pandas as pd

from func_.report_helpers import *

## Load Experimental Data
origEps2 = [0.05, 0.1, 0.2, 0.25, 0.3, 0.5, 0.7, 0.8, 0.95]
IBIs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanIBI')
sdIBI = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','SD')
CVs = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanCV')
Amps = pd.read_excel('../ABC_bursting/data/meanFeatures.xlsx','meanAmps')
DatamIBI =np.nanmean(IBIs,0)*1000
DatamCV =np.nanmean(CVs,0)


# Load Posterioes
Post = []
Eps = [0.05,0.1,0.2,0.3,0.7,0.5,0.7,0.8,0.92]#[0.05,0.2,0.5,0.7,0.8,0.92]#,0.5,0.8]#[0.05,0.1,0.2,0.3,0.5,0.65,0.8,0.95]
for epsilon in Eps:
    #name =  'ABC_results/nuki_fit_%s.npy'%epsilon
    #name =  'ABC_results/test_fit_%s.npy'%epsilon
    name =  'ABC_results/nuki_fit_LastStep%s.npy'%epsilon
    try:
        posterior_= np.load(name, allow_pickle = True)
    except:
        name =  'ABC_results/nuki_fit_%s.npy'%epsilon
        posterior_= np.load(name, allow_pickle = True)
    posterior_ = [posterior_[na([p[4] for p in posterior_])>0]]
    Post.append(posterior_[0])

# Plot convergence
plt.figure(figsize=(6,4))
sns.set_context('talk')
labels =Eps#[0.05,0.2,0.5,0.8,0.95]
ind__ = -1
for ind,pmc_posterior in enumerate(Post):
    eps = [i['epsilon'] for i in pmc_posterior]
    plt.plot(eps,'o-',label=' %s%% network, final $\\epsilon$ = %s' % (np.int(labels[ind]*100),np.round(eps[ind__],5)))
    plt.xlabel('Step', fontsize=24)
    plt.axhline(0,linestyle = '--')
    plt.ylabel(r'$\epsilon$', fontsize=30)
    plt.xlim(-1,30)
plt.legend(bbox_to_anchor=(.4, .6),fontsize =10)#loc='center right', bbox_to_anchor=(1.2, -0.6))
plt.show()
sns.despine(trim=1)
plt.tight_layout()
# plt.savefig('figs/fitting.pdf',fmt='pdf',bbox_inches = 'tight')

# Get posteriors 
X=[]
for post in Post:
    i =-1
    nu = post[i][0][0]
    ki= post[i][0][1]
    x1 = pd.Series(nu, name="$nu^{ext}$")
    x2 = pd.Series(ki, name="$K^I$")
    X.append([x1,x2])

# plot K^I distribution
p_dist = na(np.floor(X[2][1]))/200
N =1000
NI = N*na(Eps)
NE = N-NI
KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
N= 1000
# KE = 0.5490304776651601*((NE*0.14159415833146138)+32.43607416673515)
balance = KE/4
colors = sns.color_palette('Reds', len(balance))
NIs =1000*na(Eps)
plt.figure(figsize=(5,7))
Best_params  = []
for i,xs in enumerate(X):
    Best_params.append(find_MAP_2D(X[i]))
    all_kde = stats.gaussian_kde(xs[1])
    fake_x = NIs[i]*p_dist
    fake_kde = stats.gaussian_kde(fake_x)
    y_ = all_kde.evaluate(np.arange(-50,150,0.5))
    x  = np.arange(-50,150,0.5)-balance[i]
    fake_y = fake_kde.evaluate(np.arange(-50,150,0.5))
    plt.plot(x,0.5*(fake_y/fake_y.max())-i, '-', color = 'mistyrose',alpha =1)
    plt.fill_between(x,-i,0.5*(fake_y/fake_y.max())-i,color= 'mistyrose',alpha =1.)
    plt.plot(x,0.5*(y_/y_.max())-i,color =colors[i], label = '%s%%'%na(Eps[i]*100))
plt.ylabel('Networks')
plt.yticks(np.arange(-len(X)+1,1,1),['%s%%'%(i) for i in na(Eps)*100][::-1])
sns.despine()
plt.axvline(0,color = 'k',label = 'Equal Inputs', linestyle = '--',alpha =0.2)
plt.tight_layout()
plt.xlim([-30,140])
plt.xlabel('Inhibitory degrees \n relative to balance')
# plt.savefig('figs/postKI.pdf',fmt='pdf',bbox_inches = 'tight')
plt.show()



# iris = sns.load_dataset("iris")
def plot_2Dposterior(X):
    """ 
    plot 2D posteriors from ABC
    Args:
            X (array of pd series): values of accepted simualtions

    Returns:
            plot: grid plot of distributions
    """
    dat = pd.DataFrame(data= na(X).T, columns =[x.name for x in X])
    dat['$nu^{ext}$']=  dat['$nu^{ext}$']/1000
    #dat['$nu^{ext}$']=  np.log(dat['$nu^{ext}$'])
    g = sns.PairGrid(dat)
    g.map_upper(plt.scatter,s = .5)
    g.map_lower(sns.kdeplot,shade_lowest=False)#shade=True,
    g.map_diag(sns.kdeplot, lw=3, legend=False);
    g.fig.subplots_adjust(top=0.9)
    #g.fig.suptitle('%s%%'%(na(Eps[i])*100), fontsize = 16)
    plt.show()

plot_2Dposterior(X[0])


data = {}
data['DatamIBI']=DatamIBI
data['DatamCV'] =DatamCV
data['origEps2']=origEps2
data['Post'] = Post
data['X'] = X
data['Eps'] = Eps
data['i']= [-1]*len(Eps)

def get_activity(eps,sim_type = 'ML', plot=True,
                index=False):
    """ 
    Get MAP activity from the fitted model
    Args:
            eps (float): Percentage of inhibitory neurons
            type (string: 'ML' or 'MAP'): if ML try to retrieve the best model from the accepted paramters, MAP - simulate a model with MAP parameters
            plot(bool): plot raster
            index(int): which simulations out of accepted to plot, if 
    Returns:
        ST (tuple of 2 lists): spike times,gids
        sc (array): spike counts
        signal (array): Ca-activity
        bursts (list): burst times (see MI_burst)
    """
    # TODO: read paramters from a a simulation file
    #PAramters of the simualtions
    params = load_obj('sim/nuki_fit/%s/params'%eps)
    # get the data
    sim_index= np.where(na(Eps)==eps)[0][0]
    print(sim_index)
    #data_index= np.where(na(origEps2)==eps)[0][0]
    data_index = np.where(my_floor(na(origEps2),1)==my_floor(eps,1))[0][0]
    vs= True
    if sim_type=='MAP':
        print('estimating MAP...')
        best_params =find_MAP_2D(X[sim_index])#Best_params[i]#
        best_params = np.round(best_params,2)
        vs = False# need to run simualtions
    elif index:
        print('accepted simulation', index)
        best_params =Post[sim_index][-1][0][:,index]
    else:
        min_dist = np.argmin(Post[sim_index][-1][1])
        best_params =Post[sim_index][-1][0][:,min_dist]

    model =nuki_fit(params,orig_data=[DatamIBI[data_index],DatamCV[data_index]],
            alpha=1,beta =1,visual =vs)
    print('parameters:',best_params)
    params['epsilon']=eps
    # Generate data
    dat = model.generate_data(na([best_params[0],best_params[1]]))
    print('IBI,CV:',dat[0],dat[1])
    # Print the distance
    E = model.distance_function(model.summary_stats(model.orig_data),model.summary_stats(dat))
    print('Distance from Data',E)
    # Load the data
    name = get_hash(model.params)
    st,gid = read_gdf(params['directory'],name,(0,model.sim_time),threads = model.threads)
    st = st[gid<200]
    gid = gid[gid<200]
    bursts = MI_bursts(st,maxISIstart=4.5,maxISIb=4.5, minSburst= 50);

    bin_size  = 20
    sc,times = np.histogram(st,np.arange(0,model.sim_time,bin_size))
    times = times[1:]
    msc =np.median(sc)
    print('Median SC:',msc)
    tau =2000/bin_size
    signal  = np.convolve(sc,np.exp(-np.arange(10000/bin_size)/tau),'full')
    signal = signal[:-499]#cut the extra kernel width
    ST = (st,gid)
    burst_amps = get_amplitudes(signal,times,bursts)
    bursty = check_burstAmps(sc,times,bursts,thr = 3)
    print('amps are good',bursty)
    #  Amps.append(get_amplitudes(sc,time,bursts))
    if plot:
        plt.close()
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(st,gid,'|',markersize = 1.0)
        #plt.plot(time,0.2*params['N_rec']*signal/signal.max())
        plt.xlabel('time(ms)')
        plt.ylabel('neurons')
        if bursts and np.median(sc)<10.:
            bb = na(bursts)[:,0]
            b1 = na(bursts)[:,1]
            for i,b in enumerate(bb):
                plt.axvspan(b,b1[i],0,210,alpha =0.5)
            plt.plot(bb,np.ones(len(bb)),'*')
            plt.plot(b1,np.ones(len(b1))*201,'.k')
        plt.subplot(2,1,2)
        plt.plot(times,sc)
        #plt.xlim([200000,400000])
        #plt.ylim([0,params['N_rec']])
        plt.show()
    return ST,sc,signal,bursts,dat,burst_amps


# Visualize the fit
#Get all IBIs and CV from the posterior
def get_AllFeatures(eps,post,N_postSamples=50):
    """ 
    foo 
    Args:
            eps (float): inhibitiory percentage
            post (np object): ABC fit (output)
    Returns:
            allData (np array): FeaturesX N_postSamples 
    """
    params = load_obj('sim/nuki_LastStep/%s/params'%eps)
    # get the data
    sim_index= np.where(na(Eps)==eps)[0][0]
    print(eps,sim_index)
    data_index = np.where(my_floor(na(origEps2),1)==my_floor(eps,1))[0][0]
    vs= True
    model =nuki_fit(params,orig_data=[DatamIBI[data_index],DatamCV[data_index]],
            alpha=1,beta =1,visual =vs)
    # get data
    allData = np.zeros([2,N_postSamples])
    for index in range(N_postSamples):
        best_params =Post[sim_index][-1][0][:,index]
        dat = model.generate_data(na([best_params[0],best_params[1]]))
        allData[:,index] = dat
    return allData

AllData = np.zeros([2,50,len(Eps)])
for i,eps in enumerate(Eps):
    AllData[:,:,i] = get_AllFeatures(eps,i)

IBIs = na(IBIs)
CVs  = na(CVs)
fig, ax = plt.subplots(1,2)
for i,eps in enumerate(Eps):
    data_index = np.where(my_floor(na(origEps2),1)==my_floor(eps,1))[0][0]
    xs = np.random.normal(i,0.1,size=(50,))
    ax[0].plot(xs,AllData[0,:,i]/1000,'.',markersize = 2)
    xs_real = np.ones(len(IBIs[:,data_index]))*i
    ax[0].plot(xs_real,IBIs[:,data_index],'o',markersize= 3)
    #plt.xlim([0.5,1.5])
    ax[1].plot(xs,AllData[1,:,i],'.',markersize = 2)
    xs_real = np.ones(len(CVs[:,data_index]))*i
    ax[1].plot(xs_real,CVs[:,data_index],'o',markersize = 3)
    #plt.xlim([0.5,1.5])
plt.show(block=False)


# MAP vs real data
sim_type = 'ML'
mapData = np.zeros([2,len(Eps)])
mapAmp = np.zeros([len(Eps)])
par = np.zeros([len(Eps),4])
Signal = []#List of Convolved traces
ST = []
for i,eps in enumerate(Eps):
    if sim_type=='balanced':
        in_deg = balance[i]#+1
        #if eps==0.05:
        #    in_deg = balance[i]-2
    else:
        in_deg = None
    try:
        res= get_activity(eps,data,sim_type=sim_type, index=
                          0,network_type =
                          'nuki_LastStep',in_degs=in_deg,plot=False, verbosity =
                          False)
    except:
        res= get_activity(eps,data,sim_type=sim_type, index=
                          0,network_type =
                          'nuki_fit',in_degs=in_deg,plot=False, verbosity =
                          False)
    print('N of bursts:', len(res[4]))
    mapData[:,i] = res[-3]
    mapAmp[i] = np.mean(res[-2])
    par[i,1]  = res[-1][0]
    par[i,2]=int(res[-1][1])
    par[i,3] =np.round(int(res[-1][1])-balance[i],2)
    par[i,0]= eps
    Signal.append(res[2])
    ST.append(res[0])
times = res[3]

# MAP vs real data
mapData = np.zeros([2,len(Eps)])
mapAmp = np.zeros([len(Eps)])
for i,eps in enumerate(Eps):
    res= get_activity(eps,sim_type='ML', plot=False)
    mapData[:,i] = res[-2]
    mapAmp[i] = np.mean(res[-1])

# Plot MAP
n_samples= [np.sum(~np.isnan(IBIs[:,i])) for i in range(len(origEps2))]

fig, ax = plt.subplots(2,2)
for i,eps in enumerate(Eps):
    data_index = np.where(my_floor(na(origEps2),1)==my_floor(eps,1))[0][0]
    xs = np.random.normal(eps,0.01,size=(50,))
    ax[0][1].plot(xs,AllData[0,:,i]/1000,'.',markersize = 2)
    xs_real = np.ones(len(IBIs[:,data_index]))*eps
    ax[0][1].plot(xs_real,IBIs[:,data_index],'o',markersize= 3)
    #plt.xlim([0.5,1.5])
    ax[1][0].plot(xs,AllData[1,:,i],'.',markersize = 2)
    xs_real = np.ones(len(CVs[:,data_index]))*eps
    ax[1][0].plot(xs_real,CVs[:,data_index],'o',markersize = 3)
ax[0][1].errorbar(origEps2,np.nanmean(IBIs,0),np.nanstd(IBIs,0))
ax[0][1].plot(Eps,mapData[0,:]/1000,'o')
#ax[0].plot(Eps,AllData[0,:,:].T/1000,'.r',markersize = 2)
ax[1][0].errorbar(origEps2,np.nanmean(CVs,0),np.nanstd(CVs,0))
ax[1][0].plot(Eps,mapData[1,:],'o')
#ax[1].plot(Eps,AllData[1,:,:].T,'.r',markersize = 2)
ax[1][1].errorbar(origEps2,np.nanmean(Amps/Amps[0.2].mean(),0),np.nanstd(Amps/Amps[0.2].mean(),0))
ax[1][1].plot(Eps,mapAmp/mapAmp[1],'or')
fig.tight_layout()
plt.show(block=False)

# Plot Amplitudes

