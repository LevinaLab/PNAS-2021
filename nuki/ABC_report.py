# Generate a report for fitted adLIF
import sys
sys.path.insert(0,'/home/ovinogradov/Projects/EI_project/ABC_bursting/')
from func_.report_helpers import *
from matplotlib.backends.backend_pdf import PdfPages

force =False#Force recalculation of the all data

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
Eps = [0.05,0.1,0.2,0.3,0.5,0.7,0.8,0.92]
for epsilon in Eps:
    name =  'ABC_results/nuki_fit_%s.npy'%epsilon
    posterior_= np.load(name, allow_pickle = True)
    posterior_ = [posterior_[na([p[4] for p in posterior_])>0]]
    Post.append(posterior_[0])

# Get posteriors 
X=[]
selected_post = []
for post in Post:
    if False:
        acc_rate = [50/i['n total'] for i in post]
        acc_rate[0] = 100#avoid the first step
        i = np.argmin(acc_rate)
        selected_post.append(i)
    else:
        i =-1
        selected_post.append(i)
    nu = post[i][0][0]
    ki= post[i][0][1]
    x1 = pd.Series(nu, name="$nu^{ext}$")
    x2 = pd.Series(ki, name="$K^I$")
    X.append([x1,x2])

print(selected_post)
data = {}
data['DatamIBI']=DatamIBI
data['DatamCV'] =DatamCV
data['origEps2']=origEps2
data['Post'] = Post
data['X'] = X
data['Eps'] = Eps
data['i']= selected_post
# Get accepted samples from the last ABC step
try:
    AllData = np.load('ABC_results/AllDataReport.npy')
    if force:
        raise Warning('forced recalculation')
    
except:
    AllData = np.zeros([2,50,len(Eps)])
    for i,eps in enumerate(Eps):
        AllData[:,:,i] = get_AllFeatures(eps,data)
    np.save('ABC_results/AllDataReport.npy',AllData)

with PdfPages('/home/ovinogradov/Projects/EI_project/ABC_bursting/figs/report3.pdf') as pdf:

    sns.plotting_context('paper')
    # Plot convergence
    plt.figure(figsize=(12,4))
    sns.set_context('talk')
    labels =Eps#[0.05,0.2,0.5,0.8,0.95]
    ind__ = -1

    for ind,pmc_posterior in enumerate(Post):
        plt.subplot(1,2,1)
        eps = [i['epsilon'] for i in pmc_posterior]
        plt.plot(eps,'-',label=' %s%% network, final $\\epsilon$ = %s' % (np.int(labels[ind]*100),np.round(eps[ind__],5)))
        plt.xlabel('Step', fontsize=24)
        plt.axhline(0.05,linestyle = '--')
        plt.ylabel(r'$\epsilon$', fontsize=30)
        plt.xlim(-1,30)

        plt.subplot(1,2,2)
        acc_rate = [50/i['n total'] for i in pmc_posterior]
        plt.plot(acc_rate,'-',label=' %s%% network ' % (np.int(labels[ind]*100)))
        plt.ylabel('acceptance rate')
        plt.xlabel('Step', fontsize=24)
    
    plt.legend(loc='best',fontsize =10)#loc='center right', bbox_to_anchor=(1.2, -0.6))
    sns.despine(trim=1)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

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
    plt.tight_layout()
    pdf.savefig()
    plt.close()
    for i,x in enumerate(X):
        plot_2Dposterior(x,Eps[i])
        plt.tight_layout()
        pdf.savefig()
        plt.close()

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
    fig.tight_layout()
    pdf.savefig()

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
        res= get_activity(eps,data,sim_type=sim_type, index=
                          0,in_degs=in_deg,plot=False)
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
    allAmps = np.zeros([len(Eps),50])
    for i,eps in enumerate(Eps):
        for index in range(50):
            res= get_activity(eps,data,sim_type='ML', index=index,plot=False, verbosity =False)
            allAmps[i,index] = np.mean(res[-2])

    # show parameters as a table
    sns.set_context('paper')
    fig, ax = plt.subplots()
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')
    df = pd.DataFrame(par, columns=['%','$\\nu_{ext}$','$K^I$','deviation'])
    ax.table(cellText=df.values, colLabels=df.columns, loc='center',edges='open', fontsize =7)
    #fig.tight_layout()
    pdf.savefig()

    # Plot MAP
    n_samples= [np.sum(~np.isnan(IBIs[:,i])) for i in range(len(origEps2))]
    fig, ax = plt.subplots(2,2)
    fig.suptitle(sim_type)
    #Plot the model
    im='/home/ovinogradov/Projects/EI_project/ABC_bursting/figs/model.png'
    img=plt.imread(im)
    ax[0][0].imshow(img)
    ax[0][0].set_xticks([])
    ax[0][0].set_yticks([])
    ax[0][0].spines['left'].set_visible(False)
    ax[0][0].spines['right'].set_visible(False)
    ax[0][0].spines['top'].set_visible(False)
    ax[0][0].spines['bottom'].set_visible(False)

    for i,eps in enumerate(Eps):
        #IBI
        data_index = np.where(my_floor(na(origEps2),1)==my_floor(eps,1))[0][0]
        xs = np.random.normal(eps,0.01,size=(50,))
        ax[0][1].plot(xs,AllData[0,:,i]/1000,'.',markersize = 2)
        xs_real = np.ones(len(IBIs[:,data_index]))*eps
        ax[0][1].plot(xs_real,IBIs[:,data_index],'o',markersize= 3)
        ax[0][1].set_ylabel('IBI(s)')
        #plt.xlim([0.5,1.5])
        ax[1][0].plot(xs,AllData[1,:,i],'.',markersize = 2)
        xs_real = np.ones(len(CVs[:,data_index]))*eps
        ax[1][0].plot(xs_real,CVs[:,data_index],'o',markersize = 3)
        ax[0][0].set_ylabel('CV')

    ax[1][1].plot(Eps,mapAmp/mapAmp[1],'or')
    ax[0][1].errorbar(origEps2,np.nanmean(IBIs,0),np.nanstd(IBIs,0))
    ax[0][1].plot(Eps,mapData[0,:]/1000,'o')
    #ax[0].plot(Eps,AllData[0,:,:].T/1000,'.r',markersize = 2)
    ax[1][0].errorbar(origEps2,np.nanmean(CVs,0),np.nanstd(CVs,0))
    ax[1][0].plot(Eps,mapData[1,:],'o')
    #ax[1].plot(Eps,AllData[1,:,:].T,'.r',markersize = 2)
    ax[1][1].errorbar(origEps2,np.nanmean(Amps/Amps[0.1].mean(),0),np.nanstd(Amps/Amps[0.1].mean(),0))
    ax[1][1].set_ylabel('Mean burst amplitude \n (normalized to 10% network)')
    fig.tight_layout()
    pdf.savefig()
    plt.close()

    #Plot acitivity Traces
    import math
    def sigmoid(x):
       return 7 / (1 + np.exp(-x))
    plt.figure(figsize=(8,5))
    # bax = brokenaxes(ylims=((-1,shift*len(Signal)+2), (shift*len(Signal)+2, shift*len(Signal)+4)), hspace=.05)
    labels = ['%s%%'%(i*100) for i in Eps]
    position =[]
    for i,signal in enumerate(Signal):
        signal=signal.copy()
        signal/=np.max(Signal[1])#np.mean(signal)
        signal= sigmoid(signal)
        shift = 1.5#0.5#1.1
        plt.plot(times/1000,signal-(shift*i),label=labels[i])
        position.append(np.mean(signal-(shift*i)))
    plt.yticks(na(position),labels)
    sns.despine(left= True)
    plt.xlabel('time (s)')
    plt.xlim([50,500])
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    sns.set_context('talk')
    plt.figure(figsize=(10, 5))
    plt.subplot(121)
    plt.errorbar(na(origEps2)*100,np.nanmean(IBIs,0),np.nanstd(IBIs,0),fmt = 'bo-',label='Data',zorder=0)
    for i,eps in enumerate(Eps):
        plt.plot([eps*100]*50, AllData[0,:,i]/1000,'.',color = 'salmon',alpha =0.2)
        plt.errorbar([eps*100], np.mean(AllData[0,:,i]/1000),np.std(AllData[0,:,i]/1000),fmt='-or',alpha =0.9,zorder=10)
    plt.errorbar([eps*100], np.mean(AllData[0,:,i]/1000),np.std(AllData[0,:,i]/1000),fmt='-or',alpha =0.9,zorder=10,label ='Model with adaptation')
    plt.plot(na(Eps)*100,mapData[0,:]/1000,'*',label = 'MAP')
    #plt.plot(na(StaticData['perc'])*100 ,na(StaticData['mIBI']),'o-',color = 'dimgray', label ='Model without adaptation')
    plt.ylabel('IBI(s)')
    plt.xlabel('Inhibitory percentage (%%)')
    plt.legend(fontsize =14)
    sns.despine()
    plt.subplot(122)
    plt.errorbar(na(origEps2)*100,np.nanmean(CVs,0),np.nanstd(CVs,0),fmt = 'bo-',label='Data',zorder=0)
    plt.plot()
    for i,eps in enumerate(Eps):
        plt.plot([eps*100]*50, AllData[1,:,i],'.',color = 'salmon',alpha =0.2)
        plt.errorbar([eps*100], np.mean(AllData[1,:,i]),np.std(AllData[1,:,i]),fmt='-or',alpha =0.9,zorder=10)
    plt.errorbar([eps*100], np.mean(AllData[1,:,i]),np.std(AllData[1,:,i]),fmt='-or',alpha =0.9,zorder=10,label ='Model with adaptation')
    plt.plot(na(Eps)*100,mapData[1,:],'*', label ='MAP')
    #plt.plot(na(StaticData['perc'])*100 ,na(StaticData['cvIBI']),'o-',color = 'dimgray', label ='Model with adaptation')
    plt.ylabel('CV')
    plt.xlabel('Inhibitory percentage (%)')
    sns.despine()
    plt.tight_layout()
    pdf.savefig()

    sns.set_context('talk')
    # plt.figure(figsize=(5, 5))
    # EPS = Eps
    # plt.errorbar(na(origEps2)*100,na(means),stds,fmt = 'bo-',label='Data',zorder=0)
    plt.figure(figsize=(5,5))

    for i,eps in enumerate(Eps):
        plt.plot([eps*100]*50, allAmps[i,:]/mapAmp[1],'.',color = 'salmon',alpha =0.2)
    plt.plot(na(Eps)*100,mapAmp/mapAmp[1],'-or',alpha =0.9,zorder=10)
    plt.errorbar(na(origEps2)*100,np.nanmean(Amps/Amps[0.1].mean(),0),np.nanstd(Amps/Amps[0.1].mean(),0))
    plt.xlabel('Inhibitory percentage (%)')
    # amps = na(pd.read_excel('../ED-simulations/Figures/data/amp_all_export.xlsx'))
    # meanAmps = [np.mean(amps) for amps in Amps]
    # meanAmps/=meanAmps[1]
    # plt.plot(Eps,meanAmps,'-o')
    plt.ylabel('Mean amplitude  \n (normalized to 10% network)')
    sns.despine()
    plt.tight_layout()
    pdf.savefig()

    #plot spiking activity
    plt.figure()
    position=[]
    for i,st_gid in enumerate(ST):
        st,gid = st_gid
        shift = 200#1.1
        plt.plot(st,gid-(shift*i),'|',markersize = 0.25)
        position.append(np.mean(-(-100+shift*i)))
    plt.yticks(na(position),labels)
    sns.despine(left= True)
    plt.xlabel('time (s)')
    plt.xlim([50,np.max(st)])
    plt.tight_layout()
    pdf.savefig()
    plt.close()


