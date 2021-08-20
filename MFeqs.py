# sctipt to solve the mean-field quations and plot the results

from func.eq_helpers.FI import *
import sys
from func.helpers import *
from func.network import *
from func.model_setup_experimental import *
from func.ABChelpers import *
# %
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


def getTrajectories(params):
    """ 
    Sample burst trajectories from random seeds 
    """
    seeds = np.arange(1000,1010,1)
    #seeds = [1000]
    Res = []
    for seed in seeds:
        params['master_seed']=int(seed)
        time,Vs, Ws = collectVoltage(params, force = False)
        meanW = np.mean(Ws,0)
        stdW = np.std(Ws,0)
        name = get_hash(params)
        st,gid = read_gdf(params['directory'],name,(0,params['simtime']),threads = params['threads'])
        meanW_burst,meanSC_burst = meanTraj(meanW,st,gid,params,bin_size =10, primer=(20,50),
                                    interp =False,smooth=True,smooth_par=(5,3), interp_dt=0.05)
        Res.append([meanW_burst,meanSC_burst])
    return Res


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
 'voltage':False,
 'chunk': False,
 'chunk_size': 1000.0,
 'directory': 'sim/',
 'simulation': 'hash',
 'simtime':100100.0,
 'master_seed': 1000,
 'dt': 0.5,
 'threads': 50}


params['g']=4.
params['voltage']= True
#fr_stat = firing_rate(params, primer = 5000)
time,Vs, Ws = collectVoltage(params, force = False)

params['g'] =0.8
#fr_stat = firing_rate(params, primer = 500)
#FPs = get_fp(params, parallel = True)
time,Vs, Ws = collectVoltage(params, force = False)

Gs = [1,0.2]
ws = np.arange(0.,500.,0.1)
nu_overG = SimplefixedW_blockInh(params.copy(),ws,Gs,proc = 50)
simple_trajPlot(ws,nu_overG[0],Ws, params)


params['b']= 0.0
params['g'] = 0.0
ws = np.arange(0.,20.,1.)
FRs = numericalnu_overW(params,ws,primer = 1000)


Gs = [1.0]
ws = np.arange(0,20,0.01)
nu_overG = SimplefixedW_blockInh(params,ws,Gs,proc = 60)

plt.figure()
ws = np.arange(0.,20.,1.)
plt.plot(ws/20,FRs[2,:],'o')
ws = np.arange(0,20,0.1)
plt.plot(ws/20,nu_overG[0]);
plt.show(block = False)

params['g']=4.
params['voltage']= True
time,Vs, Ws = collectVoltage(params, force = False)

params['voltage']= True
params['g'] =0.8
time,Vs, Ws = collectVoltage(params, force = False)


W = []
for tau in [600.,700.,1000.,3000.]:
    params['voltage']= True
    params['g'] =0.0
    params['tau_w']= tau
    time,Vs, Ws = collectVoltage(params, force = False)
    meanW = np.mean(Ws,0)
    stdW = np.std(Ws,0)
    W.append([meanW,stdW])
    #Standard fix points
[plt.plot(w[1000:16000],label =i) for i,w in enumerate(W)]
xs = np.arange(0,len(meanW))
for i in np.arange(4):
    plt.subplot(4,1,i+1)
    #plt.errorbar(xs,W[i][0],W[i][1])
    plt.plot(W[i][1])
plt.legend()
plt.show(block = False)
FPs = get_fp(params, parallel = True)
print(FPs)

Gs = [1.,0.2]#na([params['g']])#np.arange(0.1,1.2,0.1)[::-1]
ws = np.arange(0.,30.,0.05)
nu_overG = SimplefixedW_blockInh(params,ws,Gs,proc = 60)
simple_trajPlot(ws,nu_overG[0],Ws, params)

#Fixed points over different levels of adaptation 
try:
    nu_overW = np.load('nu_overW.npy')
except:
    nu_overW = fixedW_stat(params,ws, proc = 20)
    np.save('nu_overW',nu_overW)

time,Vs, Ws = collectVoltage(params, force = False)
fr_stat = firing_rate(params, primer = 50000)

meanW = np.mean(Ws,0)
stdW = np.std(Ws,0)

#Res = getTrajectories(params)

bin_size = 100
def plot_activity(params):
    name = get_hash(params)
    st,gid = read_gdf(params['directory'],name,(0,params['simtime']),threads = params['threads'])
    plt.figure()
    plt.plot(st,gid,'.',markersize = 1.0);
    plt.xlabel('sim time')
    plt.ylabel('neuron id')
    plt.show(block=False)
    #plt.xlim([50000,100000])

meanW_burst,meanSC_burst = [],[]
for res in Res:
    for i in np.arange(0,len(res[0])):
        meanW_burst.append(res[0][i])
        meanSC_burst.append(res[1][i])

name = get_hash(params)
st,gid = read_gdf(params['directory'],name,(0,params['simtime']),threads = params['threads'])
meanW_burst,meanSC_burst = meanTraj(meanW,st,gid,params,bin_size =100, primer=(20,30),
                            interp =False,smooth=False,smooth_par=(5,3), interp_dt=0.05)

#unNu_overW,unw = unfoldNu_overW(nu_overW,ws)
unNu_overG,unw = unfoldNu_overW(nu_overG[0],ws)
#meanW_burst,meanSC_burst = Res[1]
unNu_overW = unNu_overG

i0,i1 =0,2#len(meanW_burst)
[plt.plot(mw,mw/sc_,'-',alpha =0.1) for mw,sc_ in zip(meanW_burst[i0:i1],meanSC_burst[i0:i1])]
plt.show()


def simple_trajPlot(ws,nu_overW,Ws, params):
    """ 
    Plot simulated trajectoris and MF with fixed adaptation
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """
    #Preprocess
    unNu_overW,unw = unfoldNu_overW(nu_overW,ws)
    if len(np.shape(Ws))>1:
        meanW = np.mean(Ws,0)
        stdW = np.std(Ws,0)
    else:
        meanW = Ws
    name = get_hash(params)
    st,gid = read_gdf(params['directory'],name,(0,params['simtime']),threads = params['threads'])
    meanW_burst,meanSC_burst = meanTraj(meanW,st,gid,params,bin_size =20,
                                        primer=(0,50),
                                interp =False,smooth=False,smooth_par=(5,3), interp_dt=0.05)
    plt.figure()
    plt.plot(unw/(20),unNu_overW,color = 'k',label ='analytics',linewidth = 4)
    #plt.plot(ws/20,nu_overW[0,:,:],'.',color = 'r', label ='analytics')
    #plt.plot(ws,FRs[2,:],'.',label = 'numerics');
    #plt.plot(mw_stack,sc,'-',alpha = 0.3); plt.show(block =False)
    i0,i1 =0,len(meanW_burst)
    [plt.plot(mw,sc_,'-',alpha =0.1) for mw,sc_ in zip(meanW_burst[i0:i1],meanSC_burst[i0:i1])]
#    sns.kdeplot(np.hstack(meanW_burst),np.hstack(meanSC_burst),cmap="Blues",
#                shade=True, shade_lowest=False,bw=(0.02,2.))
    plt.xlabel('fixed w')
    plt.ylabel('Rate (Hz)')
    plt.legend()
#    plt.xlim([np.min(unw/20),np.max(meanW)])
#    plt.ylim([0,np.max(mean)])
    #plt.axvline(np.max(FPs)*params['b']/ 20)
#    sns.despine(trim =True)
    plt.show(block = False)


# nu ofer G and w
#nu_overG,_= np.load('nu_overWG.npy',allow_pickle=True)


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

Gs =np.arange(0.1,1.2,0.05)[::-1]
Gs =np.arange(0.1,1.2,0.01)[::-1]
ws = np.arange(0.,1000.,0.01)
if False:
    nu_overG = SimplefixedW_blockInh(params.copy(),ws,Gs,proc = 40)
    nu_overW = SimplefixedW_blockInh(params.copy(),ws,[1.0],proc = 40)

nu_overG = bif_WsGs(params.copy(),ws,Gs,proc = 40)
#np.save('nus',[nu_overG,nu_overW])
nu_overG,nu_overW = np.load('nus.npy',allow_pickle =True)

#nu_overW=nu_overW
# Prepare Inset
s1,s2,s3 = np.shape(nu_overG)[1],np.shape(nu_overG)[0],3
image = np.ones([s1,s2,s3])
#ws = np.arange(0.,70.,0.25)
image[nu_overG[:,:,0].T>0,:] = [226/255,132/255, 133/255]#(204/255,102/255,0/255)
image[nu_overG[:,:,0].T==0,:] =[250/255,238/255,220/255]#(1.,204/255,153/255)
image[nu_overG[:,:,2].T>0,:] =[231./255,231./255,231./255]#(204./255,229/255,255/255)

# Load numerical data
#Gs_num = [1.1,1.0,0.8,0.5,0.2]#np.arange(0.1,1.2,0.1)#[::-1]

#Gs =[1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.2] #@np.arange(0.1,1.2,0.1)
Gs_num =[1.1,1.0,0.95,0.9,0.8,0.7,0.6,0.5,0.2] #@np.arange(0.1,1.2,0.1)
#ResBic = np.load('ResBic2.npy',allow_pickle=True)
ResBic = np.load('ResBicFull.npy',allow_pickle=True)
wNum = np.zeros([2,len(Gs_num)])
wErrors = np.zeros([2,len(Gs_num)])
for i,res in enumerate(ResBic):
    if len(res[0]):
        minw_ =[np.min(r) for r in res[0]]
        maxw_ =[np.max(r) for r in res[0]]
        max_w=np.mean(maxw_)
        min_w=np.mean(minw_)
        wNum[0,i] = max_w#np.nanmax(np.hstack([np.nan,[r[-1] for r in res[0]]]))
        wNum[1,i] = min_w
        wErrors[0,i] =np.std(maxw_)
        wErrors[1,i] =np.std(minw_)

plt.figure()
plt.imshow(image, aspect = 'auto',origin = 'lower', extent =
           (np.max(Gs),np.min(Gs),np.min(ws/20),np.max(ws/20)),interpolation =
          'bilinear')
plt.errorbar(Gs_num,wNum[0,:],wErrors[0,i],fmt='-o')
plt.errorbar(Gs_num,wNum[1,:],wErrors[1,i],fmt='-o')
#plt.yscale('symlog')
plt.show(block=False)


unNu_overW,unw = unfoldNu_overW(nu_overW[0],ws)
unNu_overW2,unw2 = unfoldNu_overW(nu_overG[5],ws)

# Figure 5
sns.set_context('talk')
fig, ax1 = plt.subplots(figsize = (10,8))
# These are in unitless percentages of the figure size. (0,0 is bottom left)

left, bottom, width, height = [.6, 0.6, 0.3, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])

ax1.plot(unw/20,unNu_overW,color = 'k',label ='analytics',linewidth = 4)
ax3 = ax1.twinx()
ax1.plot(unw2/20,unNu_overW2,color = 'C0',label ='analytics',linewidth = 4)
ax1.set_xlim([0.1,1.5])
ax1.set_xlim([0.1,1.5])
#ax1.set_ylim([0,150])
#ax3.set_ylim([0,400])
#ax3.spines['right'].set_edgecolor('C0')
#ax3.set_ylabel('',color = 'C0')
#meanW,meanSC = returnMeanTraj(trajs)
#ax1.plot(meanW,meanSC,color = 'k')

i0,i1 =0,20#len(meanW_burst)-20
meanW_burst,meanSC_burst = ResBic[1]
[ax1.plot(mw,sc_,'-k',alpha =0.1) for mw,sc_ in zip(meanW_burst[i0:i1],meanSC_burst[i0:i1])]
wBic,scBic = ResBic[3]
[ax3.plot(mw,sc_,'C0',alpha =0.1) for mw,sc_ in zip(wBic[i0:i1],scBic[i0:i1])]

#sns.kdeplot(np.hstack(meanW_burst),np.hstack(meanSC_burst),cmap="Reds",
#            shade=True, shade_lowest=False,bw=(0.02,2.1), ax = ax1)
#ax1.imshow(ys.T,cmap ='Reds',aspect='auto',origin='lower',extent = (np.min(w_all),np.max(w_all),np.min(sc_all),np.max(sc_all)))
#ax1.contour(ys.T,aspect='auto',origin='lower',extent =
#            (np.min(w_all),np.max(w_all),np.min(sc_all),np.max(sc_all)), cmap
#            ='Reds')
ax1.set_xlabel('adaptation current')
ax1.set_ylabel('Rate (Hz)')
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.legend()
#ax1.set_xlim([0.35,np.max(np.hstack(meanW_burst))])
#ax1.set_xlim([0.35,4])
#ax1.set_ylim([-2,110])

#plt.axvline(np.max(FPs)*params['b']/ 20)
#sns.despine(trim =True)
ax2.imshow(image, aspect = 'auto',origin = 'lower', extent =
           (np.max(Gs),np.min(Gs),np.min(ws/20),np.max(ws/20)))
ax2.plot(Gs_num,wNum[0,:],'o',color = 'gray')
ax2.plot(Gs_num,wNum[1,:],'s',color ='gray')
ax2.set_yscale('log')
ax2.set_ylabel('adaptation current')
ax2.set_xlabel('Inhibitory strength')
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.spines['top'].set_visible(False)
plt.show(block = False)

plt.savefig('Fig6v4.pdf',fmt='pdf')


ResTraj = np.load('ResTraj4.npy', allow_pickle = True)
#ResTraj = np.load('ResTraj_b0.01t8000.npy', allow_pickle = True)
#np.save('ResTraj3.npy',ResTraj)
meanW_burst,meanSC_burst = [],[]
for res in ResTraj:
    for i in np.arange(0,len(res[0])):
        meanW_burst.append(res[0][i])
        meanSC_burst.append(res[1][i])
trajs = [meanW_burst,meanSC_burst]
#meanW_burst,meanSC_burst = ResBic[0]
# Conditional KDE
w_all, sc_all = np.hstack(meanW_burst),np.hstack(meanSC_burst)
cond_ws = np.arange(np.min(w_all),np.max(w_all),0.1)
_,cond_boundaries = np.histogram(w_all, bins = 10)
xs = np.arange(np.min(sc_all)-10,np.max(sc_all),0.1)
ys = np.zeros([len(cond_boundaries),len(xs)])
for w_i,w in enumerate(cond_boundaries[:-1]):
    mask = (w_all<cond_boundaries[w_i+1])*(w_all>=w)
    y = sc_all[mask]
    kde= stats.gaussian_kde(y,bw_method = 0.9)
    ys[w_i,:] = kde.evaluate(xs)*np.log(len(y))

plt.contour(ys.T,aspect='auto',origin='lower',extent = (np.min(w_all),np.max(w_all),np.min(sc_all),np.max(sc_all)))
plt.show()


#mean trajectory
def returnMeanTraj(res):
    trajectories = np.zeros([len(res[0]),2000,2])*np.nan
    for i in range(len(res[0])):
        w = res[0][i]
        sc = res[1][i]
        trajectories[i,0:len(w),0] = w
        trajectories[i,0:len(sc),1] = sc
    meanW = np.nanmean(trajectories[:,:,0],0)
    meanSC = np.nanmean(trajectories[:,:,1],0)
    return meanW,meanSC
#meanSC = meanSC[meanW>0]
#meanW = meanW[meanW>0]

ws = np.arange(0.,30.,0.1)
all_kde = stats.gaussian_kde([np.hstack(meanW_burst),np.hstack(meanSC_burst)])#,
                      #bw_method= 0.2)
nus = np.arange(0,np.max(np.hstack(meanSC_burst)),1)
fake_dist = np.zeros([len(ws),len(nus)])
ev_par = []
for w_i,w in enumerate(ws/20):
    for nu_i,nu in enumerate(nus):
        ev_par.append([w,nu])
        #fake_dist[w_i,nu_i]= all_kde.evaluate([w,nu])

dist = all_kde.evaluate(na(ev_par).T)

meanW_burst,meanSC_burst = meanTraj(meanW,st,gid,params,bin_size =20, primer=300,
                                    interp = False, interp_dt=0.5)
#Extract the axes from the axis list
fake_dist = np.reshape(dist,[len(ws),len(nus)])
plt.plot(unw,unNu_overW,color = 'k',label ='analytics',linewidth = 4)
#plt.imshow(fake_dist.T,aspect = 'auto',origin='lower',extent = [np.min(ws/20),np.max(ws/20),0,np.max(unNu_overW)])
plt.contour(ws/20,nus,fake_dist.T)
plt.xlim([0,0.6])
plt.show(block=False)

from KDEpy import FFTKDE
import KDEpy
kde = FFTKDE(bw=0.05 ,kernel='exponential')

grid_points = 100,200
grid,points=kde.fit(na([np.hstack(meanW_burst),np.hstack(meanSC_burst)]).T).evaluate(grid_points)

x, y = np.unique(grid[:, 0]), np.unique(grid[:, 1])
z = points.reshape(grid_points[0], grid_points[1]).T

        # Plot the kernel density estimate
N = 2
plt.plot(np.hstack(meanW_burst),np.hstack(meanSC_burst),'.',alpha =0.2)
plt.contour(x, y, z, N, linewidths=0.8, colors='k')
plt.show(block = False)

kde = KernelDensity(bandwidth=0.03)
kde.fit(na([np.hstack(meanW_burst),np.hstack(meanSC_burst)]).T)
samples = kde.score_samples(na(ev_par))

#Extract the axes from the axis list
fake_dist = np.reshape(samples,[len(ws),len(nus)])
plt.plot(unw/20,unNu_overW,color = 'k',label ='analytics',linewidth = 4)
#plt.imshow(fake_dist.T,aspect = 'auto',origin='lower',extent = [np.min(ws/20),np.max(ws/20),0,np.max(unNu_overW)])
plt.contour(ws/20,nus,fake_dist.T)
plt.xlim([0,0.6])
plt.show(block=False)

plt.figure()
plt.plot(st,gid,'.',markersize = 1.0);
plt.xlabel('sim time')
plt.ylabel('neuron id')
plt.show(block=False)
#plt.xlim([50000,100000])


plt.figure()
plt.plot(ws,nu_overW[0,:,0],'-',color = 'r', label ='analytics')
plt.plot(ws,FRs[2,:],'.',label = 'numerics');
plt.xlabel('Input offset (no adaptation dynamics)')
plt.ylabel('Rate (Hz)')
plt.legend()
#plt.plot(meanW,sc,'-',alpha = 0.3); plt.show(block =False)
#plt.axvline(np.max(FPs)*params['b']/ 20)
plt.show(block = False)


primer = 0#50000
print('file:', name)
st,gid = read_gdf(params['directory'],name,(primer,params['simtime']),threads = params['threads'])
sc,times =np.histogram(st,np.arange(0,params['simtime'],1.))
plt.figure()
plt.plot(st,gid,'.',markersize = 1.0);plt.show(block=False)

#plt.plot(times[senders==2],V[senders==2],'C0',linewidth = 1);
#plt.plot(times[senders==2],V[senders==2],'r',linewidth = 1);
plt.plot(w[senders==2][times[senders==2]>100000],sc[times[senders==2]>100000]/10,'.')

plt.show(block=False)

name = get_hash(params)
print('file:', name)
st,gid = read_gdf(params['directory'],name,(5000,params['simtime']+40000),threads = params['threads'])

plt.show()
plt.plot(st,gid,'.');
plt.show()


t =np.arange(1,1000,0.01)
conds = np.arange(0,40.,.1)#[0.,1.,5.,10.,100.,]
sol = [odeint(ps_dyn, cond, t, args =(params,)) for cond in conds]
#Solving in parallel
#psDpart=partial(psD_wrapper, t=t ,w=w,params=params)
#pool = Pool(processes=20)
#sol= pool.map(psDpart, conds)
#pool.close()
plt.show()
[plt.plot(t[:1500],s[:1500]) for s in sol]
plt.yscale('symlog')
u_fp = np.unique([np.round(s[-1],2) for s in sol])
plt.show(block = False)
if len(u_fp)>1:
    interInd = np.where(np.diff(np.hstack(np.round(na(sol)[:,-1],2))))[0]
    unstableFP = conds[interInd]

fr_stat = firing_rate(params)

plt.plot(np.ones(len(u_fp)),u_fp,'.')
plt.plot([1],unstableFP,'o',color = 'gray')
plt.plot(1, fr_stat[1],'^')
plt.errorbar(1,np.mean(fr_stat[0]),np.std(fr_stat[0]), fmt = 'o-')
plt.show(block = False)

Res = prepare_data(params)
plt.plot(Res[-2],Res[-1],'.');

# Fixed adaptation 

sns.set_context('talk')
sns.set_style('ticks')
mus  = np.arange(0,60)
plt.plot([F(i,1,params) for i in mus ],label ='$\\sigma=1 $')
plt.plot([F(i,3,params) for i in mus],label ='$\\sigma=3 $')
plt.plot([F(i,10,params) for i in mus],label ='$\\sigma=10 $')
#plt.plot([0,50],[0,50],'-',color = 'gray')
plt.xlim([-5,55])
plt.ylim([-5,65])
plt.xlabel('input $\\mu$')
plt.ylabel('output rate $\\nu$')
plt.legend()
sns.despine(trim=1)
plt.tight_layout()
plt.show(block=False)

t =np.arange(0,1000,0.01)    
sol = odeint(ps_dyn, (0.2), t, args =(params,))
# scipy.integrate.RK45(ps_dyn, 0., [10.], t_bound= 100.)
# nu,a = sol.T
nu_overG = np.zeros((3,len(np.arange(0.,6.,0.1))))
nus = [450.,500.,600.]
for eta_i, eta in enumerate(nus):
    params['p_rate']= eta#eta*nu_thr
    for g_i,g in enumerate(np.arange(0.,6.,0.1)):
        params['g']= g
        sol = odeint(ps_dyn, (10.2), t, args =(params,))[-1]
        nu_overG[eta_i, g_i ]= sol
        
for i in range(len(nus)):    
    plt.plot(np.arange(0.,6.,0.1), nu_overG[i],label="Input per Ke= %s Hz"%(nus[i]) )
plt.legend()
plt.xlabel('g')
plt.ylabel('nu_0 (Hz)')
plt.tight_layout()
plt.show(block = False)
