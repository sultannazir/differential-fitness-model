import general_functions as gf
import csv
import numpy as np
from statistics import mean
import itertools

Parameters = {'H' : 15,             # Number of hosts: each host is an independent simulation
              'K1' : 5000,       # Mean value of 'within-host fitness' of microbes of type 1 in the environment
              'K2' : 1000,       # Mean value of 'within-host fitness' of microbes of type 2 in the environment
              'stdK1' : 500,      # Standard deviation of 'within-host fitness' of microbes in the environment
              'stdK2' : 500,
              'stdI' : 1,          # Standard deviation of 'within-host interaction coefficients'...
                                   # of microbes in the environment
              'env_rat1' : 0.1,    # Relative abundance of type 1 microbes in the environment = K1'/(K1' + K2')
              'init_size' : 100,   # Initial population size of each microbe type in the host(s)
              'K_min' : 100,       # Minimum value of within-host fitness any microbe can attain
              'K1_max' : 10000,  # Maximum value of within-host fitness type 1 microbes can attain
              'K2_max' : 10000,   # Maximum value of within-host fitness type 2 microbes can attain
              'mu' : 10,
              'd' : 0.0,          # Probability of death of a microbe in host at each time step
              'w' : 0.5,           # Relative effect of intraspecific interactions to interspecific interactions in
                                   # birth and death of a microbe
              'm' : 1000,           # Size of colonizing microbe population at each time step
              'sign1' : -1,        # Nature of effect of Microbe type 2 on Microbe type 1 (choose from -1,0,1)
              'sign2' : -1,        # Nature of effect of Microbe type 1 on Microbe type 2 (choose from -1,0,1)
              'b' : 0.05,         # Bottleneck ratio - fraction of number of parent's microbes inherited by offspring
              'T' : 500,          # Host generation time - time before next bottleneck event
              'sim_time' : 10000    # Simulation time
              }

def run_simulation_getdat(Parameters):
    gf.parameters_to_global_variables(Parameters)

    data1 = []  # tracks mean K1 values in time in each host
    data2 = []  # tracks mean K2 values in time in each host
    datai1 = []  # tracks mean alpha values of type 1 microbes in each host
    datai2 = []  # tracks mean alpha values of type 2 microbes in each host
    dataM1 = []  # tracks pop size if microbe 1 in each host
    dataM2 = []  # tracks pop size of microbe 2 in each host

    K1val, K2val = gf.initialize_Kval()  # initial K distributions
    I12val, I21val = gf.initialize_Ival()  # initial alpha distributions

    t = 1
    while t <= gf.sim_time:
        K1val, K2val, I12val, I21val = gf.update_microbes_new(K1val, K2val, I12val, I21val)

        if t%gf.T == 0: # bottleneck event at every Tth time step
            data1.append(['K1', t] + [*map(mean,K1val)])
            data2.append(['K2', t] + [*map(mean,K2val)])
            datai1.append(['I1', t] + [*map(mean,I12val)])
            datai2.append(['I2', t] + [*map(mean,I21val)])
            dataM1.append(['M1', t] + [*map(len, K1val)])
            dataM2.append(['M2', t] + [*map(len,K2val)])

            K1val, K2val, I12val, I21val = gf.bottleneck(K1val, K2val, I12val, I21val)
        print(t)
        t += 1

    return(data1, data2, datai1, datai2, dataM1, dataM2)

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

    OPmean = ['mean', gf.b, gf.m, gf.T, K1mean, K2mean, I1mean, I2mean, M1mean, M2mean]
    OPmed = ['med', gf.b, gf.m, gf.T, K1med, K2med, I1med, I2med, M1med, M2med]
    OPstd = ['std', gf.b, gf.m, gf.T, K1std, K2std, I1std, I2std, M1std, M2std]
    print('b: {} , m: {} , T: {} complete'.format(gf.b, gf.m, gf.T))
    return(OPmean, OPmed, OPstd)


#data1, data2, datai1, datai2, dataM1, dataM2 = run_simulation_getdat(Parameters)

#dat = np.vstack([data1, data2, datai1, datai2, dataM1, dataM2])

#with open("Data_manyhost_new.csv", "w", newline="") as f:
#    writer = csv.writer(f)
#    writer.writerows(dat)
