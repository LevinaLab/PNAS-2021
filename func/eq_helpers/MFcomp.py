
"""Helper functions to compare mean field with fixed adaptation(seperation of timescales) and network simulations"""
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

def prepare_data(params, bin_size =20):
    """ 
    minifunction to read the data, detect bursts and calculate the spike
    counts
    Args:
            params(dict): see adLIF network

    Returns:
            tuple: bursts,spike counts, spike count times
    """
    name = get_hash(params)
    print('file:', name)
    st,gid = read_gdf(params['directory'],name,(5000,params['simtime']),threads = params['threads'])
    bursts = MI_bursts(st)
    sc,times =np.histogram(st,np.arange(0,params['simtime'],bin_size))
    return bursts, sc, times[1:],st,gid

def firing_rate(params,bin_size = 20, primer=50000):
    """ 
    print Firing rate inside and outside of the burst
    Args:
            path (dict): paramters of the simualtions (see adLIF networks)
            primer (float): time to cut at the beginning of simulaitonation


    Returns:
            list: [[mean fr in each burst], mean fr outside of the burst]
    """
    force = False
    try: 
        if force:
            raise Exception('force simulate')
        name = get_hash(params)
        print('file:', name)
        st,gid = read_gdf(params['directory'],name,(primer,params['simtime']),threads = params['threads'])
    except:
        nest.set_verbosity('M_FATAL')
        A = adLIFNet(params)
        A.build()
        A.connect()
        A.run()
    name = get_hash(params)
    print('file:', name)
    N = params['N_rec']
    st,gid = read_gdf(params['directory'],name,(primer,params['simtime']),threads = params['threads'])
    meanfr = (len(st)/((params['simtime']-primer)/1000)/params['N_rec'])
    print('Mean FR %s Hz'%(len(st)/((params['simtime']-primer)/1000)/params['N_rec']))
    bursts = MI_bursts(st)
    if len(bursts)>1:
        fr_in = np.zeros(len(bursts))
        fr_out = 0
        for i,burst in enumerate(bursts):
            mask = (st>burst[0]+20.)*(st<burst[1])
            st_ = st[mask]
            fr_in[i] =len(st_)/((burst[1]-burst[0])/1000)/N
            if i>0:
                mask1 = (st>bursts[i-1][1])*(st<burst[0])
                st_ =st[mask1]
                fr_out +=len(st_)/((burst[0]-bursts[i-1][1])/1000)/N
        mFrout = fr_out/(len(bursts)-1)
        print('mean FR burst',np.mean(fr_in))
        print('mean FR quiescency',fr_out/(len(bursts)-1))
    else:
        fr_in, mFrout = np.nan,np.nan
    return fr_in, mFrout,meanfr


def psDw_wrapper(x0,w,params,t):
    """ 
    wrapper for solving pseudodynamics in parallel
    Args:
            x0 (float):initcond

    Returns:
            list: Burst times 
    """
    sol = odeint(ps_dynW, x0, t, args =(w,params,))
    return sol

def psD_wrapper(x0,params,t):
    """ 
    wrapper for solving pseudodynamics in parallel
    Args:
            x0 (float):initcond

    Returns:
            list: Burst times 
    """
    sol = odeint(ps_dyn, x0, t, args =(params,))
    return sol

def plot_FI(sigma,params):
    mus  = np.arange(-100,100,1)
    plt.figure()
    plt.plot(mus,[F(i,sigma,params) for i in mus ],label ='$\\sigma=%s$'%(sigma))
    plt.plot([np.min(mus),np.max(mus)],[0,200],'-',color = 'gray')
    #plt.plot([0,100],[0,100],color = 'gray')
    plt.show(block = False)

#Fi as a funtion of nu
def plot_FIf(params):
    plt.figure()
    nus = np.arange(0,150,0.1)
    Fs = [F(mu_(nu,params),np.sqrt(sigma_sq(nu,params)),params) for nu in nus]
    plt.plot(nus,Fs)
    plt.plot([np.min(nus),np.max(nus)],[0,np.max(nus)],'-',color = 'gray')
    plt.show(block = False)


def get_fp(params, parallel = False, plot_solutions =False):
    """ 
    get FP numerically
    Args:
            params (dict):

    Returns:
            list: fixed points
    """
    t =np.arange(1,1000,0.01)
    conds = np.arange(0,40.,1.)
    if parallel:
        pool = Pool(processes=20)
        psDpart=partial(psD_wrapper, t=t, params=params)
        sol= pool.map(psDpart, conds)
        pool.close()
    else:
        sol = [odeint(ps_dyn, cond, t, args =(params,)) for cond in conds]
    #Solving in parallel
    #pool = Pool(processes=20)
    #sol= pool.map(psDpart, conds)
    #pool.close()
    #plt.yscale('symlog')
    u_fp = np.unique([np.round(s[-1],2) for s in sol])
    FPs = [fp for fp in u_fp]
    if len(u_fp)>1:
        interInd = np.where(np.diff(np.hstack(np.round(na(sol)[:,-1],2))))[0]
        unstableFP = conds[interInd]
        FPs.append(unstableFP)
    if plot_solutions:
        plt.figure()
        [plt.plot(t[:1500],s[:1500]) for s in sol]
        plt.show(block = False)
    return FPs

def collectVoltage(params,primer = 50000,force = False):
    try:
        if force:
            raise Exception('force simulate')
        name = get_hash(params)
        print('file:', name)
        st,gid = read_gdf(params['directory'],name,(primer,params['simtime']),threads = params['threads'])
        Ws,Vs,time = np.load(params['directory']+name+'/Voltages.npy',allow_pickle = True)
    except:
        nest.set_verbosity('M_FATAL')
        A = adLIFNet(params)
        A.build()
        A.connect()
        A.run()
        name = get_hash(params)
        dmm =nest.GetStatus(A.voltmeter)
        w = dmm[0]['events']['w']
        V = dmm[0]['events']['V_abs']
        times = dmm[0]['events']['times']
        senders = dmm[0]['events']['senders']
        Ws = na([w[senders==i] for i in np.unique(senders)])
        Vs = na([V[senders==i] for i in np.unique(senders)])
        time = na(times[senders==senders[1]])
        np.save(params['directory']+name+'/Voltages.npy',[Ws,Vs,time])
    return time,Vs,Ws

#Gs = np.arange(0.1,1.5,0.1)
def fixedW_stat(params,ws,proc = 20):
    """ 
    Stationary solutions with fixed adaptaiton
    Args:
            params(dict)

    Returns:
            ws (array): values of w
            nu over w (): firing rate
    """
    t =np.arange(1,1000,0.01)
    Gs = [4]
    #ws = np.arange(0.,40.,0.5)
    nu_overW = np.zeros((len(Gs),len(ws),3))*np.nan
    #for g_i, g_mod in enumerate(Gs):
    g_i= 0
    conds = np.arange(0,40.,.1)#[0.,1.,5.,10.,100.,]
    pool = Pool(processes=proc)
    for w_i, w in enumerate(ws):
        #params['g'] = (NE/NI)*g_mod
        psDpart=partial(psDw_wrapper, t=t ,w=w,params=params)
        sol= pool.map(psDpart, conds)
        u_sol = (np.unique([np.round(s[-1],1) for s in sol]))
        if len(u_sol)>1:
            interInd = np.where(np.diff(np.hstack(np.round(na(sol)[:,-1],1))))[0]
            unstableFP = conds[interInd]
            nu_overW[g_i,w_i,2]= unstableFP
        nu_overW[g_i,w_i,0:len(u_sol)]= u_sol
    pool.close()
    return nu_overW

#Gs = np.arange(0.1,1.5,0.1)
def fixedW_blockInh(params,ws,proc = 20):
    """ 
    Stationary solutions with fixed adaptaiton
    Args:
            params(dict)

    Returns:
            ws (array): values of w
            nu over w (): firing rate
    """
    t =np.arange(1,1000,0.01)
    Gs = np.arange(0.1,1.2,0.01)[::-1]
    #ws = np.arange(0.,40.,0.5)
    nu_overW = np.zeros((len(Gs),len(ws),3))*np.nan
    g_i= 0
    conds = np.arange(0,40.,.1)#[0.,1.,5.,10.,100.,]
    pool = Pool(processes=proc)
    for g_i, g_mod in enumerate(Gs):
        for w_i, w in enumerate(ws):
            params['g'] = 4*g_mod
            psDpart=partial(psDw_wrapper, t=t ,w=w,params=params)
            sol= pool.map(psDpart, conds)
            u_sol = (np.unique([np.round(s[-1],1) for s in sol]))
            if len(u_sol)>1:
                interInd = np.where(np.diff(np.hstack(np.round(na(sol)[:,-1],1))))[0]
                unstableFP = conds[interInd]
                nu_overW[g_i,w_i,2]= unstableFP
            nu_overW[g_i,w_i,0:len(u_sol)]= u_sol
    pool.close()
    return nu_overW

def getW_FP(w,params,extra= None, extra_value = None):
    """ 
    Helper to collect fixed points
    Args:
            w (float): fixed adaptation current
            params (str): dict of network params
            extra(string): an extra parameter to set from params
            extra_value (): value of a extra paramter (type same as in params)
    Returns:
            list: Burst times
    """
    fps =np.zeros(3)*np.nan
    t =np.arange(1,1000,0.01)
    conds = [0.,500.]#np.arange(0,40.,10.)#[0.,1.,5.,10.,100.,]
    if extra:
        params[extra] = extra_value
    sol = [odeint(ps_dynW, cond, t, args =(w,params,)) for cond in conds]
    u_sol = (np.unique([np.round(s[-1],1) for s in sol]))
    u_fp = np.nan
    if len(u_sol)>1:
        #find the unstable FP with reverse time
        sol = odeint(ps_dynW, np.mean(u_sol), -t, args =(w,params,))
        u_fp = sol[-1]
    u_sol = np.hstack([u_sol,u_fp])
    fps[0:len(u_sol)] = u_sol
    return fps


def getWG_FP(WG,params,extra= None, extra_value = None):
    """ 
    Helper to collect fixed points
    Args:
            w (float): fixed adaptation current
            params (str): dict of network params
            extra(string): an extra parameter to set from params
            extra_value (): value of a extra paramter (type same as in params)
    Returns:
            list: Burst times
    """
    w,g = WG
    fps =np.zeros(3)*np.nan
    t =np.arange(1,1000,0.01)
    conds = [0.,500.]#np.arange(0,40.,10.)#[0.,1.,5.,10.,100.,]
    params['g'] = g
    sol = [odeint(ps_dynW, cond, t, args =(w,params,)) for cond in conds]
    u_sol = (np.unique([np.round(s[-1],1) for s in sol]))
    u_fp = np.nan
    if len(u_sol)>1:
        #find the unstable FP with reverse time
        sol = odeint(ps_dynW, np.mean(u_sol), -t, args =(w,params,))
        u_fp = sol[-1]
    u_sol = np.hstack([u_sol,u_fp])
    fps[0:len(u_sol)] = u_sol
    return fps


def SimplefixedW_blockInh(params,ws,Gs,proc = 20):
    """ 
    Stationary solutions with fixed adaptaiton using reserse time for the unstable FP
    Args:
            params(dict)
            ws(array)
            Gs(array)

    Returns:
            ws (array): values of w
            nu over w (): firing rate
    """
#    Gs = np.arange(0.1,1.2,0.01)[::-1]
    #ws = np.arange(0.,40.,0.5)
    nu_overW = np.zeros((len(Gs),len(ws),3))*np.nan
    g_i= 0
    pool = Pool(processes=proc)
    for g_i, g_mod in enumerate(Gs):
        params['g'] = 4*g_mod
        fpOverWp = partial(getW_FP,params=params)
        results= pool.map(fpOverWp,ws)
        for w_i, w in enumerate(ws):
            nu_overW[g_i,w_i,:]= results[w_i]
    pool.close()
    return nu_overW


def bif_WsGs(params,ws,Gs,proc = 20):
    """ 
    Stationary solutions with fixed adaptaiton using reserse time for the unstable FP
    Args:
            params(dict)
            ws(array)
            Gs(array)

    Returns:
            ws (array): values of w
            nu over w (): firing rate
    """
#    Gs = np.arange(0.1,1.2,0.01)[::-1]
    #ws = np.arange(0.,40.,0.5)
    nu_overW = np.zeros((len(Gs),len(ws),3))*np.nan
    g_i= 0
    pool = Pool(processes=proc)
    wsgs = []
    for g_i,g_mod in enumerate(Gs):
        for w_i, w in enumerate(ws):
            g= 4*g_mod
            wsgs.append((w,g))

    fpOverWp = partial(getWG_FP,params=params,)
    results= pool.map(fpOverWp,wsgs)
    #for g_i,g_mod in enumerate(Gs):
    #    for w_i, w in enumerate(ws):
    #        nu_overW[g_i,w_i,:]= results[w_i]
    pool.close()
    return results

def numericalnu_overW(params,ws,primer = 50000):
    """ 
    Check how fixed W affects the firing rate
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """
    #ws = np.arange(0,40,0.5)
    FRs = np.zeros([3,len(ws)])
    for i,w in enumerate(ws):
        params['constantI'] = w/params['tauMem']
        fr = firing_rate(params, primer = primer)
        FRs[0,i]=np.mean(fr[0])
        FRs[1,i] = fr[1]
        FRs[2,i] = fr[2]
    return FRs

from scipy import interpolate
from scipy.signal import savgol_filter
def meanTraj(meanW,st,gid,params,dt_v=1.0,bin_size =20, primer =(0,0), interp
            = False, smooth= False,smooth_par = (10,3),interp_dt = 0.1):
    
    """ 
    Get the mean trajectories of adaptation and firing rate over bursts
    Args:
        meanW(array): mean adaptation over time
        st(list): spike times
        gid(list): neurons ids
        dt_v(float): time resolution of adaptation
        bin_size(float): time bin for the mean trajectories (ms)
        primer(float): ms to add before and after the burst Returns:
        interpolate(bool): upsample for better kde 
        smooth(bool): smooth with savgol filter
        smooth_par (tuple)
        interp_dt(float): time steps for interpolations (ms)
    Return:
            meanW_burst: mean adaptation variale over bursts
            meanSC_burst: mean spike count over bursts
    """

    bursts = MI_bursts(st)
    sc,times =np.histogram(st,np.arange(0,params['simtime'],bin_size))
    sc = sc/params['N_rec']/(bin_size/1000)
    times = times[1:]
    rec_bin = np.int(bin_size/dt_v)
    mw_stack = [np.mean(meanW[i:i+rec_bin]) for i in range(0,len(meanW)-rec_bin,rec_bin)]
    mw_stack = [meanW[i:i+rec_bin][-1] for i in range(0,len(meanW)-rec_bin,rec_bin)]
    mw_stack = na(mw_stack)
    meanW_burst = []
    meanSC_burst =[]
    primer1,primer2 = primer
    for burst in bursts[1:]:#ignore the first burst

        mask = (times>burst[0]-primer1)*(times<burst[1]+primer2)
        w_burst = mw_stack[mask]
        sc_burst = sc[mask]
        if interp:
            times_ = times[mask]
            f = interpolate.interp1d(times_, sc_burst)
            f1 = interpolate.interp1d(times_, w_burst)
            times_new = np.arange(np.min(times_),np.max(times_),interp_dt);
            sc_new =f(times_new)
            w_new =f1(times_new)
            meanW_burst.append(w_new)
            meanSC_burst.append(sc_new)
        elif smooth:
            timest,deg = smooth_par
            w_new = savgol_filter(w_burst,timest,deg)
            sc_new = savgol_filter(sc_burst, timest,deg)
            meanW_burst.append(w_new)
            meanSC_burst.append(sc_new)
        else:
            meanW_burst.append(w_burst)
            meanSC_burst.append(sc_burst)
    return meanW_burst,meanSC_burst


def unfoldNu_overW(nu_overW,ws):
    """ 
    unfold 3 stationary solutions of nu over W into 1 vector
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """


    segmentA =np.hstack([nu_overW[nu_overW[:,0]>1,0],nu_overW[nu_overW[:,1]>1,1]])
    wsA = np.hstack([ws[nu_overW[:,0]>1],ws[nu_overW[:,1]>0]])
    #segmentBA = nu_overW[0,nu_overW[0,:,2]>0,2]
    #wsBA = ws[nu_overW[0,:,2]>0]
    segmentB = nu_overW[nu_overW[:,2]>0,2][::-1]
    wsB = ws[nu_overW[:,2]>0][::-1]
    segmentC = nu_overW[nu_overW[:,0]<5,0]
    wsC = ws[nu_overW[:,0]<5]
    unNu_overW = np.hstack([segmentA,segmentB,segmentC])
    unw = np.hstack([wsA,wsB,wsC])
    return unNu_overW,unw
    #for i in range(len(nu_overW)):
        

