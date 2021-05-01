import general_functions as gf
import csv
import numpy as np

Parameters = {'H' : 1,             # Number of hosts
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
              'muK' : 50,
              'muI' : 0.01,
              'd' : 0.0,          # Probability of death of a microbe in host at each time step
              'w' : 0.5,           # Relative effect of intraspecific interactions to interspecific interactions in
                                   # birth and death of a microbe
              'm' : 5000,           # Size of colonizing microbe population at each time step
              'sign1' : -1,        # Nature of effect of Microbe type 2 on Microbe type 1 (choose from -1,0,1)
              'sign2' : -1,        # Nature of effect of Microbe type 1 on Microbe type 2 (choose from -1,0,1)
              'b' : 0.01,         # Bottleneck ratio - fraction of number of parent's microbes inherited by offspring
              'T' : 12,          # Host generation time - time before next bottleneck event
              'sim_time' : 1999    # Simulation time
              }

def run_simulation_getdist(Parameters):
    gf.parameters_to_global_variables(Parameters)

    data1 = []  # tracks all K1 values in time in the first host (for now there is only one host)
    data2 = []  # tracks all K2 values in time in the first host
    datai1 = []  # tracks all alpha values of type 1 microbes in the first host
    datai2 = []  # tracks all alpha values of type 2 microbes in the first host
    K1val, K2val = gf.initialize_Kval()  # initial K distributions
    I12val, I21val = gf.initialize_Ival()  # initial alpha distributions

    t = 1
    while t <= gf.sim_time:
        K1val, K2val, I12val, I21val = gf.update_microbes_new(K1val, K2val, I12val, I21val)

        if t%gf.T == 0: # bottleneck event at every Tth time step
            K1val, K2val, I12val, I21val = gf.bottleneck(K1val, K2val, I12val, I21val)

        data1.append(K1val[0])
        data2.append(K2val[0])
        datai1.append(I12val[0])
        datai2.append(I21val[0])
        print(t)
        t += 1

    return(data1, data2, datai1, datai2)

data1, data2, datai1, datai2 = run_simulation_getdist(Parameters)

# write data
# can change data filenames here

with open("Data/K1.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(data1)

with open("Data/K2.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(data2)

with open("Data/I1.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(datai1)

with open("Data/I2.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerows(datai2)