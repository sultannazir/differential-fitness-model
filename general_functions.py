import numpy as np
from scipy import stats
import random

def parameters_to_global_variables(Parameters):
# reads input parameters from dictionary and assigns them as global variables with variable names same as the dict keys
    keys = list(Parameters.keys())
    for i in keys:
        globals()[i] = Parameters[i]

def trunc_samples(mu, sigma, lower, upper, num_samples):
# samples some num_sumples = (n,m) sized list of lists of numbers from normal(mu, sigma), bounded below and above,
# and renormalized
    n = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    samples = [list((n.rvs(num_samples[0])).flatten()) for i in range(num_samples[1])]
    return samples

def initialize_Kval():
# inital K1 and K2 value distributions in H hosts
    K1val = trunc_samples(K1, stdK, K_min, K1_max, (init_size, H))
    K1val = [[int(num) for num in j] for j in K1val]
    K2val = trunc_samples(K2, stdK, K_min, K2_max, (init_size, H))
    K2val = [[int(num) for num in j] for j in K2val]

    return(K1val, K2val)

def initialize_Ival():
# inital alpha value distribution
    I12val = trunc_samples(0, stdI, 0, 1, (init_size, H))
    I12val = [[round(num, 3) for num in j] for j in I12val]
    I21val = trunc_samples(0, stdI, 0, 1, (init_size, H))
    I21val = [[round(num, 3) for num in j] for j in I21val]

    return(I12val, I21val)


def death(Kvals, Iijvals):
# inputs are K and alpha distributions of a single microbe type of a single host
    num_mic = len(Kvals)
    probs = np.random.random(num_mic)
    selected = [i for i in range(len(Kvals)) if probs[i] < d] # select microbes to die

    Kvals = [v for i,v in enumerate(Kvals) if i not in frozenset(set(selected))]
    # redifine microbe values excluding dead microbes
    Iijvals = [v for i, v in enumerate(Iijvals) if i not in frozenset(set(selected))]

    return(Kvals, Iijvals)

def birth_propensity(Kvals, Iijval, Mi, Mj, sign):
    birth = w*(1 - len(Kvals)/np.array(Kvals))
    interaction = (1 - w) * sign * np.array(Iijval) * Mj / (Mi + Mj)
    return(birth + interaction)

def colonize(K1val, K2val, I12val, I21val):

    n1 = int(np.random.binomial(m, env_rat1, 1))  # number of type 1 microbes in colonizing population
    n2 = m - n1  # number of type 2 microbes in colonizing population

    #samples K and alpha values from corresponding distributions
    newK1 = trunc_samples(K1, stdK, K_min, K1_max, (n1, 1))[0]
    newK1 = [int(num) for num in newK1]
    newI12 = trunc_samples(0, stdI, 0, 1, (n1, 1))[0]
    newI12 = [round(num, 3) for num in newI12]
    newK2 = trunc_samples(K2, stdK, K_min, K2_max, (n2, 1))[0]
    newK2 = [int(num) for num in newK2]
    newI21 = trunc_samples(0, stdI, 0, 1, (n2, 1))[0]
    newI21 = [round(num, 3) for num in newI21]

    # add the new microbes to the host
    K1val = K1val + newK1
    I12val = I12val + newI12
    K2val = K2val + newK2
    I21val = I21val + newI21

    return(K1val, K2val, I12val, I21val)

def update_microbes(K1val, K2val, I12val, I21val):

    for h in range(len(K1val)):

        # colonization
        K1val[h], K2val[h], I12val[h], I21val[h] = colonize(K1val[h], K2val[h], I12val[h], I21val[h])

        # random death
        K1val[h], I12val[h] = death(K1val[h], I12val[h])
        K2val[h], I21val[h] = death(K2val[h], I21val[h])

        # birth + interaction values
        prop1 = birth_propensity(K1val[h],I12val[h], len(I12val), len(I21val), sign1)
        prop2 = birth_propensity(K2val[h],I21val[h], len(I21val), len(I12val), sign2)

        # death by antagonistic interactions outweighing birth propensity
        selected_die1 = [i for i in range(len(prop1)) if prop1[i] < 0]
        selected_die2 = [i for i in range(len(prop2)) if prop2[i] < 0]

        K1val[h] = [v for i, v in enumerate(K1val[h]) if i not in frozenset(set(selected_die1))]
        I12val[h] = [v for i, v in enumerate(I12val[h]) if i not in frozenset(set(selected_die1))]
        prop1 = [v for i, v in enumerate(prop1) if i not in frozenset(set(selected_die1))]
        K2val[h] = [v for i, v in enumerate(K2val[h]) if i not in frozenset(set(selected_die2))]
        I21val[h] = [v for i, v in enumerate(I21val[h]) if i not in frozenset(set(selected_die2))]
        prop2 = [v for i, v in enumerate(prop2) if i not in frozenset(set(selected_die2))]

        # birth
        probs1 = np.random.random(len(prop1))
        probs2 = np.random.random(len(prop2))

        selected_birth1 = [i for i in range(len(prop1)) if prop1[i] > probs1[i]]
        selected_birth2 = [i for i in range(len(prop2)) if prop2[i] > probs2[i]]

        # create replicates of reproducing microbes
        for i in selected_birth1:
            K1val[h].append(K1val[h][i])
            I12val[h].append(I12val[h][i])
        for i in selected_birth2:
            K2val[h].append(K2val[h][i])
            I21val[h].append(I21val[h][i])

    return(K1val, K2val, I12val, I21val)

def bottleneck(K1val, K2val, I12val, I21val):

    N = int(np.ceil(b*(len(K1val) + len(K2val))))  # total number of microbes transmitted vertically
    n1 = int(np.random.binomial(N, len(K1val)/(len(K1val) + len(K2val))))  # number of type 1 microbes
    n2 = N - n1  # number of type 2 microbes

    # select microbes from parent microbiome
    selected1 = random.sample(range(0, len(K1val)), n1)
    selected2 = random.sample(range(0, len(K2val)), n2)

    K1val = [v for i, v in enumerate(K1val) if i in frozenset(set(selected1))]
    I12val = [v for i, v in enumerate(I12val) if i in frozenset(set(selected1))]
    K2val = [v for i, v in enumerate(K2val) if i in frozenset(set(selected2))]
    I21val = [v for i, v in enumerate(I21val) if i in frozenset(set(selected2))]

    return(K1val, K2val, I12val, I21val)