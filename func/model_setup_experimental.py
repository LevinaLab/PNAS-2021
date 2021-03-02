#WORKING VERSTION OF THE ABC MODEL 25.02.19
from func.simple_abc_truncatedGaus import Model, basic_abc, pmc_abc
import numpy as np
from scipy import stats
import sys
from func.helpers import *
from func.network import *
na = np.array
# plt.style.use('ggplot')
from func.ABChelpers import *
import sys
import pandas as pd


class LoacError(Exception): pass
class MyModel(Model):
    """ 
    Basis for ABC model fit
    #This method initializes the model object. In this case it does nothing, but you can have you model 
    #do something when the model object is created, like read in a table of something. 
    ___init___
    Args:
            params (dict): paramters of the simualations (see network for
            desciption)
            orig_data (list): IBI and CV of the experimantal data
            bin_size (float): bin_size of for spike count calculation (relevant
            for burstiness measure)
            visual (boolean): if True avoid simulations and allow only loading
            from memory
            alpha (float): scaling of the IBI factor in the objective
            beta (float): scaling of the CV factor in the objective
    """
    verbose = False

    def __init__(self,params,orig_data, bin_size=20,threads = 20, visual =
                 False, alpha =1, beta =1):
        self.params= params
        self.path = params['directory']
        self.threads = params['threads']
        self.sim_time = params['simtime']
        self.bin_size = bin_size
        self.orig_data = orig_data.copy()
        self.visual = visual
        self.alpha = alpha
        self.beta = beta
    #This is the method that draws from you prior. In this example it draws from frozen scipy.stats 
    #distributions that are set with the Model.set_priors method.
    def draw_theta(self):
        theta = []
        for p in self.prior:
            theta.append(p.rvs())
        return theta
    
    #The method that generates the synthetic data sets.
    def generate_data(self, theta):
        #print(theta)
        nu_ext,j_ext,KI = theta[0], theta[1],theta[2]
        self.params['J_ext'] =  j_ext
        self.params['p_rate'] =nu_ext
        self.params['Ks'][1] = int(KI)
        #print(self.params)
        try:
            name = get_hash(self.params)
            st,gid = read_gdf(sef.params['directory'],name,(5000,self.sim_time),threads = self.threads)
        except:
            if self.visual:
                raise LoacError('wont load')
            else:
                A = adLIFNet(self.params)
                A.build()
                A.connect()
                A.run()
                if self.verbose:
                    print('sim...')
                try:
                    name = get_hash(self.params)
                    st,gid = read_gdf(sef.params['directory'],name,(5000,self.sim_time),threads =  self.threads)
                    if self.verbose:
                        print('loaded')
                except:
                    st,gid = [0],[0]

        burstiness = False
        if len(st)>100:
            bin_size = 20
            sc,times =np.histogram(st,np.arange(0,self.params['simtime'],bin_size))
            times = times[1:]
            if np.median(sc)<10:
                burstiness = True
        if burstiness:
            if self.verbose:
                print('bursty')
            bursts = MI_bursts(st)
            #second check of burstiness
            #makes sure that there is at leas one burst 3 sigma above the
            #threshold
            burstiness=  check_burstAmps(sc,times,bursts)
        else:
            bursts = [0]

        if bursts and burstiness:
            ibis =na(bursts)[1:,0]-na(bursts)[:-1,1]
            #dur = np.mean(na(bursts)[:,1]- na(bursts)[:,0])
            mibis = np.mean(ibis)
            cvs = np.std(ibis)/mibis
        else:
            mibis = np.nan
            cvs = np.nan
#
        return na([mibis,cvs])

    def get_origData(self,theta, epsilon):
        j,nu_ext,j_ext,KI = theta[0], theta[1],theta[2],theta[3]
        self.params['J'] = j
        self.params['J_ext'] =  j_ext#theta[2]#j/10
        NI = N*epsilon
        NE = N-NI
        KE = NE*p
#         KI = int(KE/self.params['g'])
        self.params['epsilon'] =epsilon
        self.params['p_rate'] =nu_ext*N#Scaled by Number of neurons #(ex_rate/16)*((N-(epsilon*N))*0.1)
        self.params['Ks'] = (int(KE),int(KI))
        try: 
            name = get_hash(params)
            st,gid = read_gdf(self.params['directory'],name,(5000,sim_time),threads = self.threads)
        except:

            A = adLIFNet(self.params); 
            A.build()
            A.connect()
            A.run()
        
            try:
                st = nest.GetStatus(A.espikes,keys='events')[0]['times']
                gid = nest.GetStatus(A.espikes,keys='events')[0]['senders']
                #name = get_hash(self.params)
                #st,gid = read_gdf(self.path,name,(5000,self.sim_time),threads =  self.threads)
            except:
                st,gid = [0],[0]
        bursts = MI_bursts(st)
        return bursts
    
    #The method that computes your summary statistics, for a Gaussian we just need mu and sigma
    def summary_stats(self, data):
        # normalization
        return (data[0]/self.orig_data[0],data[1]/self.orig_data[1])
    
    #And finally the distance function. We are just going to use the euclidean distance 
    #from our observed summary stats
    def distance_function(self, data, synth_data):
        loss = (self.alpha*np.nanmean((data[0]- synth_data[0])**2)+
                (self.beta*np.nanmean((data[1]- synth_data[1])**2)))
#        print(data)
#        print(synth_data)
        if self.verbose:
            print(loss/2)

        return loss/2

class nuki_fit(MyModel):
    """ 
    ABC fit of adLIF network with nu_ext and K^I as paratmers 
    (expliting the fact the J_ext depends on the nu_ext)
    updates get data method
    """
    def generate_data(self, theta):
        """ 
        generate data for ABC
        Args:
                theta (list): nu_ext an K^I

        Returns:
                bursts (list): Burst times detected with MI
        """

        if self.verbose:
            print(theta)
        nu_ext,KI = theta[0],theta[1]
        self.params['p_rate'] =nu_ext
        self.params['Ks'][1] = int(KI)
        try:
            name = get_hash(self.params)
            st,gid = read_gdf(self.path,name,(5000,self.sim_time),threads = self.threads)
        except:
            if self.visual:
                raise LoacError('wont load')
            else:
                A = adLIFNet(self.params)
                A.build()
                A.connect()
                A.run()
                try:
                    name = get_hash(self.params)
                    st,gid = read_gdf(self.path,name,(5000,self.sim_time),threads =  self.threads)
                except:
                    st,gid = [0],[0]
        burstiness = False
        if len(st)>100:
            bin_size = 20
            sc,times =np.histogram(st,np.arange(0,self.params['simtime'],bin_size))
            times = times[1:]
            if np.median(sc)<10:
                burstiness = True

        if burstiness:
            bursts = MI_bursts(st)
            #second check of burstiness
            #makes sure that there is at leas one burst 3 sigma above the
            #threshold
            burstiness=  check_burstAmps(sc,times,bursts)
        else:
            bursts = [0]

        if bursts and burstiness:
            ibis =na(bursts)[1:,0]-na(bursts)[:-1,1]
            dur = np.mean(na(bursts)[:,1]- na(bursts)[:,0])
        else:
            ibis = 0
            dur = 0
        mibis = np.mean(ibis)
        cvs = np.std(ibis)/mibis
        if self.verbose:
            print('data:',mibis,cvs)

        return na([mibis,cvs])
