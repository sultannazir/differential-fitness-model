import numpy as np
import itertools
from joblib import Parallel, delayed
from pathlib import Path
import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import general_functions as gf

data_folder = Path("Data_ParStudy/")
dataName = "data_HV_inv.csv"

Parameters = {'H' : 1,             # Number of hosts: each host is an independent simulation
              'K1' : 5000,       # Mean value of 'within-host fitness' of microbes of type 1 in the environment
              'K2' : 1000,       # Mean value of 'within-host fitness' of microbes of type 2 in the environment
              'stdK1' : 500,      # Standard deviation of 'within-host fitness' of microbes in the environment
              'stdK2' : 500,
              'stdI' : 1,          # Standard deviation of 'within-host interaction coefficients'...
                                   # of microbes in the environment
              'env_rat1' : 0.9,    # Relative abundance of type 1 microbes in the environment = K1'/(K1' + K2')
              'init_size' : 100,   # Initial population size of each microbe type in the host(s)
              'K_min' : 100,       # Minimum value of within-host fitness any microbe can attain
              'K1_max' : 10000,  # Maximum value of within-host fitness type 1 microbes can attain
              'K2_max' : 10000,   # Maximum value of within-host fitness type 2 microbes can attain
              'd' : 0.0,          # Probability of death of a microbe in host at each time step
              'muK' : 50,           # mutation parameters fro K and Alpha
              'muI' : 0.01,
              'w' : 0.5,           # Relative effect of intraspecific interactions to interspecific interactions in
                                   # birth and death of a microbe
              'm' : 100,           # Size of colonizing microbe population at each time step
              'sign1' : 1,        # Nature of effect of Microbe type 2 on Microbe type 1 (choose from -1,0,1)
              'sign2' : 1,        # Nature of effect of Microbe type 1 on Microbe type 2 (choose from -1,0,1)
              'b' : 0.05,         # Bottleneck ratio - fraction of number of parent's microbes inherited by offspring
              'T' : 500,          # Host generation time - time before next bottleneck event
              'sim_time' : 1999    # Simulation time
              }


# Set parameter ranges
par1 = 'b'
par2 = 'm'
par3 = 'T'
range1 = np.array([0.01])
range2 = np.array([10, 100, 1000])
range3 = np.array([5, 50, 500])

def set_parameters(p1, p2, p3):
    Parameters_local = Parameters.copy()
    Parameters_local[par1] = p1
    Parameters_local[par2] = p2
    Parameters_local[par3] = p3

    return Parameters_local

def run_simulation_finalstate(Parameters):
    gf.parameters_to_global_variables(Parameters)

    K1val, K2val = gf.initialize_Kval()  # initial K distributions
    I12val, I21val = gf.initialize_Ival()  # initial alpha distributions

    t = 1
    while t <= gf.sim_time:
        K1val, K2val, I12val, I21val = gf.update_microbes_new(K1val, K2val, I12val, I21val)

        if t % gf.T == 0:  # bottleneck event at every Tth time step
            K1val, K2val, I12val, I21val = gf.bottleneck(K1val, K2val, I12val, I21val)

        t += 1

    K1all = list(itertools.chain.from_iterable(K1val))
    K1mean, K1med, K1std = np.mean(K1all), np.median(K1all), np.std(K1all)
    K2all = list(itertools.chain.from_iterable(K2val))
    K2mean, K2med, K2std = np.mean(K2all), np.median(K2all), np.std(K2all)
    I1all = list(itertools.chain.from_iterable(I12val))
    I1mean, I1med, I1std = np.mean(I1all), np.median(I1all), np.std(I1all)
    I2all = list(itertools.chain.from_iterable(I21val))
    I2mean, I2med, I2std = np.mean(I2all), np.median(I2all), np.std(I2all)

    M1all = [*map(len, K1val)]
    M2all = [*map(len, K2val)]
    M1mean, M1med, M1std = np.mean(M1all), np.median(M1all), np.std(M1all)
    M2mean, M2med, M2std = np.mean(M2all), np.median(M2all), np.std(M2all)

    OPmean = ['mean', Parameters[par1], Parameters[par2], Parameters[par3], K1mean, K2mean, I1mean, I2mean, M1mean, M2mean]
    OPmed = ['med', Parameters[par1], Parameters[par2], Parameters[par3], K1med, K2med, I1med, I2med, M1med, M2med]
    OPstd = ['std', Parameters[par1], Parameters[par2], Parameters[par3], K1std, K2std, I1std, I2std, M1std, M2std]
    print(par1,': {} , '.format(Parameters[par1]),par2,': {} , '.format(Parameters[par2]),par3,': {} complete'.format(Parameters[par3]))
    return(OPmean, OPmed, OPstd)



def run_model():
    # set modelpar list to run
    modelParList = [set_parameters(*x)
                    for x in itertools.product(*(range1, range2, range3))]

    # run model
    nJobs = min(len(modelParList), -1)
    print('starting with %i jobs' % len(modelParList))
    results = Parallel(n_jobs=nJobs, verbose=9, timeout=1.E9)(
        delayed(run_simulation_finalstate)(par) for par in modelParList)

    # store output
    Output = zip(*results)
    statData = np.vstack(Output)
    print(statData)

    saveName = data_folder / dataName
    with open(saveName, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(statData)


#run_model()


#### reading data

cmap = sns.color_palette('coolwarm', as_cmap=True)

cols = ['stat', par1, par2, par3, 'K1', 'K2', 'I1', 'I2', 'M1', 'M2']
data = pd.read_csv(data_folder / dataName, names=cols)

data = data[data['T'] != 1]

data_mean = data[cols[1:]].loc[data['stat'] == 'mean']
data_med = data[cols[1:]].loc[data['stat'] == 'med']
data_std = data[cols[1:]].loc[data['stat'] == 'std']

data_mean['rat'] = data_mean['M1']/(data_mean['M1']+data_mean['M2'])
d0 = data_mean[[cols[i] for i in [2,3]]+['rat']]
d0 = d0.pivot('m','T','rat')
sns.heatmap(d0, vmin=0, vmax=1, cmap=cmap)
plt.title('M1/M1+M2')
plt.show()
#### plotting data

cmap = sns.color_palette('Spectral', as_cmap=True)
fig, ax = plt.subplots(ncols=6,nrows=2, figsize=(15.4,5.5))

#plt.subplot(2,3,1)
d1 = data_mean[[cols[i] for i in [2,3,4]]]
d1 = d1.pivot('m','T','K1')
g = sns.heatmap(d1, vmin=0, vmax=10000, ax=ax[0][0], cmap=cmap)
ax[0][0].set_title('K1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,4)
d2 = data_mean[[cols[i] for i in [2,3,5]]]
d2 = d2.pivot('m','T','K2')
g = sns.heatmap(d2, vmin=0, vmax=10000, ax=ax[0][1], cmap=cmap)
ax[0][1].set_title('K2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,2)
d5 = data_mean[[cols[i] for i in [2,3,6]]]
d5 = d5.pivot('m','T','I1')
g = sns.heatmap(d5, vmin=0, vmax=1, ax=ax[0][2], cmap=cmap)
ax[0][2].set_title('Alpha1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,5)
d6 = data_mean[[cols[i] for i in [2,3,7]]]
d6 = d6.pivot('m','T','I2')
g = sns.heatmap(d6, vmin=0, vmax=1, ax=ax[0][3], cmap=cmap)
ax[0][3].set_title('Alpha2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,3)
d9 = data_mean[[cols[i] for i in [2,3,8]]]
d9 = d9.pivot('m','T','M1')
g = sns.heatmap(d9, vmin=0, vmax=10000, ax=ax[0][4], cmap=cmap)
ax[0][4].set_title('M1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,6)
d10 = data_mean[[cols[i] for i in [2,3,9]]]
d10 = d10.pivot('m','T','M2')
g = sns.heatmap(d10, vmin=0, vmax=10000, ax=ax[0][5], cmap=cmap)
ax[0][5].set_title('M2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,1)
d1 = data_std[[cols[i] for i in [2,3,4]]]
d1 = d1.pivot('m','T','K1')
g = sns.heatmap(d1, vmin=0, vmax=1000, ax=ax[1][0], cmap=cmap)
ax[1][0].set_title('K1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,4)
d2 = data_std[[cols[i] for i in [2,3,5]]]
d2 = d2.pivot('m','T','K2')
g = sns.heatmap(d2, vmin=0, vmax=1000, ax=ax[1][1], cmap=cmap)
ax[1][1].set_title('K2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,2)
d5 = data_std[[cols[i] for i in [2,3,6]]]
d5 = d5.pivot('m','T','I1')
g = sns.heatmap(d5, vmin=0, vmax=0.5, ax=ax[1][2], cmap=cmap)
ax[1][2].set_title('Alpha1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,5)
d6 = data_std[[cols[i] for i in [2,3,7]]]
d6 = d6.pivot('m','T','I2')
g = sns.heatmap(d6, vmin=0, vmax=0.5, ax=ax[1][3], cmap=cmap)
ax[1][3].set_title('Alpha2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,3)
d9 = data_std[[cols[i] for i in [2,3,8]]]
d9 = d9.pivot('m','T','M1')
g = sns.heatmap(d9, vmin=0, vmax=1000, ax=ax[1][4], cmap=cmap)
ax[1][4].set_title('M1')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

#plt.subplot(2,3,6)
d10 = data_std[[cols[i] for i in [2,3,9]]]
d10 = d10.pivot('m','T','M2')
g = sns.heatmap(d10, vmin=0, vmax=1000, ax=ax[1][5], cmap=cmap)
ax[1][5].set_title('M2')
g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)
g.set_xticklabels(g.get_xticklabels(), fontsize = 8)

plt.suptitle('Mean value (above) and standard deviation (below) of evolvable parameters at t=1999\nK1 = 5000, K2 = 1000, K1max = 10000, '
             'K2max = 10000\nstdK = 500, b = 0.01, muK = 50, muI = 0.01, K10/K20 = 9\n Interaction: -/-')
fig.tight_layout()
plt.show()