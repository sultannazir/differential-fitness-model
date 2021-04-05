import general_functions as gf
import pandas as pd
import statistics as st

Parameters = {'H' : 1,             # Number of hosts
              'K1' : 500000,       # Mean value of 'within-host fitness' of microbes of type 1 in the environment
              'K2' : 100000,       # Mean value of 'within-host fitness' of microbes of type 2 in the environment
              'stdK' : 50000,      # Standard deviation of 'within-host fitness' of microbes in the environment
              'stdI' : 1,          # Standard deviation of 'within-host interaction coefficients'...
                                   # of microbes in the environment
              'env_rat1' : 0.2,    # Relative abundance of type 1 microbes in the environment = K1'/(K1' + K2')
              'init_size' : 100,   # Initial population size of each microbe type in the host(s)
              'K_min' : 100,       # Minimum value of within-host fitness any microbe can attain
              'K1_max' : 1000000,  # Maximum value of within-host fitness type 1 microbes can attain
              'K2_max' : 200000,   # Maximum value of within-host fitness type 2 microbes can attain
              'd' : 0.01,          # Probability of death of a microbe in host at each time step
              'w' : 0.5,           # Relative effect of intraspecific interactions to interspecific interactions in
                                   # birth and death of a microbe
              'm' : 100,           # Size of colonizing microbe population at each time step
              'sign1' : -1,        # Nature of effect of Microbe type 2 on Microbe type 1 (choose from -1,0,1)
              'sign2' : -1,        # Nature of effect of Microbe type 1 on Microbe type 2 (choose from -1,0,1)
              'b' : 0.001,         # Bottleneck ratio - fraction of number of parent's microbes inherited by offspring
              'T' : 1000,          # Host generation time - time before next bottleneck event
              'sim_time' : 999    # Simulation time
              }


def run_simulation_getmean(Parameters):
    gf.parameters_to_global_variables(Parameters)
    data = []

    K1val, K2val = gf.initialize_Kval()
    I12val, I21val = gf.initialize_Ival()

    t = 1
    while t <= gf.sim_time:
        K1val, K2val, I12val, I21val = gf.update_microbes(K1val, K2val, I12val, I21val)

        if t%gf.T == 0:
            K1val, K2val, I12val, I21val = gf.bottleneck(K1val, K2val, I12val, I21val)

        if t%10 == 1:
            meanK1 = st.mean(K1val[0])
            meanK2 = st.mean(K2val[0])
            stdK1 = st.stdev(K1val[0])
            stdK2 = st.stdev(K2val[0])
            meanI1 = st.mean(I12val[0])
            meanI2 = st.mean(I21val[0])
            stdI1 = st.stdev(I12val[0])
            stdI2 = st.stdev(I21val[0])
            M1 = len(K1val[0])
            M2 = len(K2val[0])

            data.append((M1,meanK1, stdK1, meanI1, stdI1, M2,meanK2, stdK2, meanI2, stdI2))
            print(t, M1, M2)
        t += 1

    cols = ['M1','meanK1', 'stdK1', 'meanI1', 'stdI1', 'M2', 'meanK2', 'stdK2', 'meanI2', 'stdI2']
    df = pd.DataFrame(data, columns=cols)

    return(df)

Data = run_simulation_getmean(Parameters)

Data.to_csv('sample3.csv')