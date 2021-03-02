from scipy import stats
import sys
import numpy as np
na = np.array
from scipy import stats
import sys
import numpy as np
from scipy.signal import find_peaks
from numba import jit
import pickle
na = np.array
#check for bimodality 
#check for bimodality 

def my_floor(a, precision=0):
    return np.round(a - 0.5 * 10**(-precision), precision)

def save_obj(obj, name ):
    """ 
    general utility to save python files
    Args:
            obj: Name of the object to save
	    name: saving dir and name

    """
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    """ 
    General utility to load objects in python
    Args:
            file name(str): dir and file name

    Returns:
            python object
    """

    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


def check_burstAmps(sc,times,bursts, thr = 3, bin_size = 20):
    """ 
    check if amplitudes of bursts are n sigma [1] above the noise level
    Even single amp above the threshold is enough to mark the simulation as
    bursty
    The function uses exp kernel with 2s timescale (FLUO-4 indicator timescale)
    The use of these thresholds is standard in denoising algorithms based on
    the wavelet transform (Donoho, 1995)
    see also Quiroga et al 2004
    Args:
            sc (array): spike counts
            time (array): signal time
            bursts (list of tuples): burst start and end
            thr  (str): n of sigmas above the threshold needed to accept
            bin_size (int): bin size of the spike count (in ms)
    Returns:
            bursiness (bool): True or False(if non bursty)
    """
    tau =2000/bin_size
    signal  = np.convolve(sc,np.exp(-np.arange(10000/bin_size)/tau),'full')
    signal = signal[:-499]
    sigma = np.median(signal/0.6745)
    burst_amps = get_amplitudes(signal,times,bursts)
    if np.any(burst_amps>sigma*thr):
        return True
    else:
        return False

def b_check(isi):
    """ 
    burstiness check
    Checks if log isi distribution is bimodal
    Args:
            [1] (array) isi (or spike )
    Returns:
            1 if bursty, 0 otherwise
    """
    if len(isi)<1:
        return 0
    data = np.log(isi[isi>0])
    all_kde = stats.gaussian_kde(data,bw_method=0.25)    
    x= np.arange(-1,np.max(data),0.1)
    y_ = all_kde.evaluate(x)
    indexes,_ = find_peaks(y_, distance=5)
    maxind = np.argmax(y_)
    rmaxsc =np.round(x[maxind],2)
    allpeaks = np.round(x[indexes],2) 
    maxsc = allpeaks[rmaxsc==allpeaks]

    if np.any(allpeaks>maxsc):
        return 1
    else:
        return 0

from collections import Counter
#@jit#(nopython=True)
def MI_bursts_old(st,
              maxISIstart=4.5,
              maxISIb=4.5,
              minBdur=10,
              minIBI=20,
              minSburst=100):
    """Min Interval method [1,2] for burst detections
    stable version from before 03.03.2020

    Args:
            [1] spike times (list) (signle neuron or popultation) (ms)
            [2] (float) max ISI at start of burst (ms)
            [3] (float) max ISI in burst (ms)
            [4] (float) min burst duration (ms)
            [5] (float) min inter-burst interval (ms)
            [6] (float) min number of spike in a burst 
    Returns:
            burst (list of tuples): burst start, burst end
    [1] Nex Technologies.NeuroExplorer Manual.  Nex Technologies,2014
    [2] Cotterill, E., and Eglen, S.J. (2018). Burst detection methods. ArXiv:1802.01287.        """
    isi_pop = np.diff(np.sort(st))
    spikes = np.sort(st)
    b_spikes = []
    burst_ = []
    # find in-burst ISIs
    # plt.plot(inBurst)
    b_start = 0
    for i,s in enumerate(spikes[:-1]):
        if isi_pop[i]<maxISIstart or b_start>1:#50:#maxISIstart:
            b_start += 1
            if isi_pop[i]<=maxISIb:
                b_spikes.append(s)
        #         if isi_pop[i]< maxISIstart:
        #             bStart.append(s)
            elif len(b_spikes)>=minSburst:
                #uni,counts= np.unique(np.round(b_spikes,-1),return_counts=True)
                counts = Counter(np.round(b_spikes,-1))
                counts = list(counts.values())
                if np.any([c>50 for c in counts]):# np.any(counts>50):
                    burst_.append((b_spikes[0],b_spikes[-1]))
                    b_spikes = []
                    b_start= 0
                else:
                    b_spikes = []
                    b_start = 0
            else:
                b_spikes = []
                b_start = 0
        else:
            b_spikes= []
            b_start = 0
    bursts = []
    if burst_:
        bursts.append(burst_[0])
        for i,b in enumerate(burst_[1:]):
            if b[1]-b[0]>=minBdur and b[0]-bursts[-1][1] >= minIBI:
                bursts.append(b)
            elif b[0]-bursts[-1][1]<= minIBI:
                bursts[-1] = (bursts[-1][0],b[1])
    return bursts

@jit(nopython=True)
def find_burstlets(spikes,r_spikes,isi_pop,
              maxISIstart=4.5,
              maxISIb=4.5,
              minSburst=100):
    """ 
    Helper to find burstlets
    Args:
        spikes (arr): spike times
        r_spikes(arr): rounded spike times
        isi_pop(arr):isi

    Returns:
            burst_ (list of tuples): Burst start, burst end
    """
    b_spikes = 0
    burst_ = []
    sync_b = False
    b_start = 0
    b_size = 0
    for i,s in enumerate(spikes[:-1]):
        if isi_pop[i]<maxISIstart and b_start==0:#50:#maxISIstart:
            b_size = 0
            b_start += 1
            b_spikes=s
        elif isi_pop[i]<=maxISIb and b_start>0: #start if two conseq init isi
            b_start+=1
            if r_spikes[i-1]==r_spikes[i]:#dynamicaly check burst size
                #equal spike times come from the different electrodes
                b_size+=1
            #    print(b_size)
                if b_size>minSburst:
                    sync_b = True
            else:
                b_size = 0
            #else:#reset burst size if it in a new burstlet
            #    b_size =0
        elif b_start>=minSburst and sync_b:
            burst_.append((b_spikes,s))
            b_spikes =None
            b_size = 0
            sync_b = False
            b_start = 0

        else:
            b_spike =None
            b_size = 0
            sync_b = False
            b_start = 0
    return burst_

#from collections import Counter
#@jit(nopython=True)
def MI_bursts(st,
              maxISIstart=4.5,
              maxISIb=4.5,
              minBdur=40,
              minIBI=40,
              minSburst=50):
    """Min Interval method [1,2] for burst detections
    Optimized version from 03.03.20
    OV

    Args:
            [1] spike times (list) (signle neuron or popultation) (ms)
            [2] (float) max ISI at start of burst (ms)
            [3] (float) max ISI in burst (ms)
            [4] (float) min burst duration (ms)
            [5] (float) min inter-burst interval (ms)
            [6] (float) min number of spike in a burst 
    Returns:
            burst (list of tuples): burst start, burst end
    [1] Nex Technologies.NeuroExplorer Manual.  Nex Technologies,2014
    [2] Cotterill, E., and Eglen, S.J. (2018). Burst detection methods. ArXiv:1802.01287.        """
    spikes = np.sort(st)
    r_spikes = np.round_(spikes,-1)
    isi_pop = np.diff(spikes)
    burst_ =find_burstlets(spikes,r_spikes,isi_pop,maxISIstart,maxISIb,minSburst)
    bursts = []
    if burst_:
        bursts.append(burst_[0])
        for i,b in enumerate(burst_[1:]):
            if b[1]-b[0]>=minBdur and b[0]-bursts[-1][1] >= minIBI:
                bursts.append(b)
            elif b[0]-bursts[-1][1]<= minIBI:
                bursts[-1] = (bursts[-1][0],b[1])
    return bursts

#@jit(nopython=True)
def get_amplitudes(sc,times,bursts):
    """get amplitudes from total spike count using burst indicies detected with adapted MI method 
    (see MI_bursts)

    Args:
            [1] spike counts (list) (signle neuron or popultation) (ms)
            [2] time (list or array) time in ms for spike count (should match the bursts)
            [3] bursts (list) detected bursts: list of tuples (burst_start, burst_end) [output of MI_bursts]
    Returns:
            Amp (list of tuples): burst amplitudes
    """
    if bursts:
        bb = na(bursts)[:,0]
        b1 = na(bursts)[:,1]
        Amp = []

        for i,b in enumerate(bb):
            amp = np.max(sc[na(times>b-100) * na(times<b1[i]+100)])
            Amp.append(amp)
        return Amp 
    else:
        return np.nan



from scipy.optimize import minimize as minimize
from scipy.optimize import basinhopping
 
def find_MAP_2D(data):
    """ 
    find MAP in 2D
    Args:
            path (str): Directory of the simulation with "/".

    Returns:
            list: Burst times 
    """
    x= na([data[0],data[1]])
    all_kde = stats.gaussian_kde(x)
    x_list= np.linspace(np.min(x[0]),np.max(x[0]),20)
    y_list = np.linspace(np.min(x[1]),np.max(x[1]),20)
    x_list, y_list = np.meshgrid(x_list, y_list)
    def neg_kde(x):
        return -all_kde.evaluate(x)
    zs = neg_kde(na([x_list.ravel(), y_list.ravel()]))
    opt = basinhopping(neg_kde,[np.mean(x[0]),np.mean(x[1])],
                       minimizer_kwargs={'method':"Nelder-Mead"})
    return opt.x

def find_MAP(data):
    x= na([data[0],data[1],data[2]])
    all_kde = stats.gaussian_kde(x)
    x_list= np.linspace(np.min(x[0]),np.max(x[0]),20)
    y_list = np.linspace(np.min(x[1]),np.max(x[1]),20)
    z_list = np.linspace(np.min(x[2]),np.max(x[2]),20)
    x_list, y_list,z_list = np.meshgrid(x_list, y_list, z_list)
    def neg_kde(x):
        return -all_kde.evaluate(x)
    zs = neg_kde(na([x_list.ravel(), y_list.ravel(),z_list.ravel()]))
    opt = basinhopping(neg_kde,[np.mean(x[0]),np.mean(x[1]),np.mean(x[2])],
                       minimizer_kwargs={'method':"Nelder-Mead"})
    return opt.x
