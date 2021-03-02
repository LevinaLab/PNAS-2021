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
# MAP 
sim_type = 'sample'#'MAP'
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
                      0,sample = True,in_degs=False,plot=False)
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
#pdf.savefig()
plt.show(block=False)
#plt.close()

