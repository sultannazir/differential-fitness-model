import numpy as np
import itertools
from joblib import Parallel, delayed
import simulation_manyhost as sm
from pathlib import Path
import csv
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data_folder = Path("Data_ParStudy/")
dataName = "data.csv"

Parameters = {'H' : 15,             # Number of hosts: each host is an independent simulation
              'K1' : 5000,       # Mean value of 'within-host fitness' of microbes of type 1 in the environment
              'K2' : 1000,       # Mean value of 'within-host fitness' of microbes of type 2 in the environment
              'stdK1' : 100,      # Standard deviation of 'within-host fitness' of microbes in the environment
              'stdK2' : 1000,
              'stdI' : 1,          # Standard deviation of 'within-host interaction coefficients'...
                                   # of microbes in the environment
              'env_rat1' : 0.1,    # Relative abundance of type 1 microbes in the environment = K1'/(K1' + K2')
              'init_size' : 100,   # Initial population size of each microbe type in the host(s)
              'K_min' : 100,       # Minimum value of within-host fitness any microbe can attain
              'K1_max' : np.inf,  # Maximum value of within-host fitness type 1 microbes can attain
              'K2_max' : np.inf,   # Maximum value of within-host fitness type 2 microbes can attain
              'd' : 0.01,          # Probability of death of a microbe in host at each time step
              'w' : 0.5,           # Relative effect of intraspecific interactions to interspecific interactions in
                                   # birth and death of a microbe
              'm' : 100,           # Size of colonizing microbe population at each time step
              'sign1' : -1,        # Nature of effect of Microbe type 2 on Microbe type 1 (choose from -1,0,1)
              'sign2' : -1,        # Nature of effect of Microbe type 1 on Microbe type 2 (choose from -1,0,1)
              'b' : 0.05,         # Bottleneck ratio - fraction of number of parent's microbes inherited by offspring
              'T' : 500,          # Host generation time - time before next bottleneck event
              'sim_time' : 9999    # Simulation time
              }

num_bins = 2

# Ste parameter ranges
par1 = 'stdK2'
par2 = 'm'
par3 = 'env_rat1'
range1 = np.array([100, 1000])
range2 = np.array([100, 1000])
range3 = np.array([0.1, 0.01])

def set_parameters(p1, p2, p3):
    Parameters_local = Parameters.copy()
    Parameters_local[par1] = p1
    Parameters_local[par2] = p2
    Parameters_local[par3] = p3

    return Parameters_local


def run_model():
    # set modelpar list to run
    modelParList = [set_parameters(*x)
                    for x in itertools.product(*(range1, range2, range3))]

    # run model
    nJobs = min(len(modelParList), -1)
    print('starting with %i jobs' % len(modelParList))
    results = Parallel(n_jobs=nJobs, verbose=9, timeout=1.E9)(
        delayed(sm.run_simulation_finalstate)(par) for par in modelParList)

    # process and store output
    Output = zip(*results)
    statData = np.vstack(Output)
    print(statData)

    saveName = data_folder / dataName
    with open(saveName, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(statData)


run_model()


#### reading data

cols = ['stat', 'stdK2', 'm', 'env_rat1', 'K1', 'K2', 'I1', 'I2', 'M1', 'M2']
data = pd.read_csv(data_folder / dataName, names=cols)

data_mean = data[cols[1:]].loc[data['stat'] == 'mean']
data_med = data[cols[1:]].loc[data['stat'] == 'med']
data_std = data[cols[1:]].loc[data['stat'] == 'std']

#### plotting data

fig, ax = plt.subplots(3,2, figsize=(5.5,7.7))

#plt.subplot(2,3,1)
d1 = data_mean[[cols[i] for i in [1,2,4]]].loc[data_mean['env_rat1'] == 0.1]
d1 = d1.pivot('m','stdK2','K1')
sns.heatmap(d1, vmin=0, vmax=8000, ax=ax[0][0])
ax[0][0].set_title('K1')

#plt.subplot(2,3,4)
d2 = data_mean[[cols[i] for i in [1,2,5]]].loc[data_mean['env_rat1'] == 0.1]
d2 = d2.pivot('m','stdK2','K2')
sns.heatmap(d2, vmin=0, vmax=8000, ax=ax[0][1])
ax[0][1].set_title('K2')

#plt.subplot(2,3,2)
d5 = data_mean[[cols[i] for i in [1,2,6]]].loc[data_mean['env_rat1'] == 0.1]
d5 = d5.pivot('m','stdK2','I1')
sns.heatmap(d5, vmin=0, vmax=0.3, ax=ax[1][0])
ax[1][0].set_title('Alpha1')

#plt.subplot(2,3,5)
d6 = data_mean[[cols[i] for i in [1,2,7]]].loc[data_mean['env_rat1'] == 0.1]
d6 = d6.pivot('m','stdK2','I2')
sns.heatmap(d6, vmin=0, vmax=0.3, ax=ax[1][1])
ax[1][1].set_title('Alpha2')

#plt.subplot(2,3,3)
d9 = data_mean[[cols[i] for i in [1,2,8]]].loc[data_mean['env_rat1'] == 0.1]
d9 = d9.pivot('m','stdK2','M1')
sns.heatmap(d9, vmin=0, vmax=8000, ax=ax[2][0])
ax[2][0].set_title('M1')

#plt.subplot(2,3,6)
d10 = data_mean[[cols[i] for i in [1,2,9]]].loc[data_mean['env_rat1'] == 0.1]
d10 = d10.pivot('m','stdK2','M2')
sns.heatmap(d10, vmin=0, vmax=8000, ax=ax[2][1])
ax[2][1].set_title('M2')

plt.suptitle('Mean value of microbes from 15 replicates at t=9999\nK1 = 5000, K2 = 1000, K1max = inf, K2max = inf\nstdK1 = 100, K10/K20 = 1/9, d = 0.01, b = 0.05, T = 500')
fig.tight_layout()
plt.show()