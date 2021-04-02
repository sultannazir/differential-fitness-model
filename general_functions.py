import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def trunc_samples(mu, sigma, lower, upper, num_samples):
    n = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    samples = [list((n.rvs(num_samples[0])).flatten()) for i in range(num_samples[1])]
    return samples

def initialize_Kval(K1,K2,stdK,init_size, K_min_clip, K1_max_clip, K2_max_clip, M):
# inital K1 and K2 value distributions in M hosts
    K1val = trunc_samples(K1, stdK, K_min_clip, K1_max_clip, (M, init_size))
    K2val = trunc_samples(K2, stdK, K_min_clip, K2_max_clip, (M, init_size))

    return(K1val, K2val)

def initialize_Ival(stdI, init_size):
# inital alpha value distribution I12mat = effects of microbe 1 on microbe 2
    I12val = trunc_samples(0, stdI, -1, 1, (M, init_size))
    I21val = trunc_samples(0, stdI, -1, 1, (M, init_size))

    return(I12val, I21val)


def death(Kvals, Iijvals, d):

    num_mic = len(Kvals)
    probs = np.random.random(num_mic)
    selected = [i for i in range(len(Kvals)) if probs[i] < d] # select microbes to die

    Kvals = [v for i,v in enumerate(Kvals) if i not in frozenset(set(selected))]
    Iijvals = [v for i, v in enumerate(Iijvals) if i not in frozenset(set(selected))]

    return(Kvals, Iijvals)

def birth_propensity(Kvals, w):
    prop = w*(1 - len(Kvals)/np.array(Kvals)) # birth propensities
    return(prop)

def colonize(K1val, K2val, I12val, I21val, env_rat1, m, K1, K2, stdK, stdI, K_min_clip, K1_max_clip, K2_max_clip):
    colprob = np.random.random()
    if colprob < m:
        prob = np.random.random()
        if prob < env_rat1:
            K1val.append(trunc_samples(K1,stdK,K_min_clip, K1_max_clip,(1,1))[0][0])
            I12val.append(trunc_samples(0, stdI, -1, 1, (1, 1))[0][0])
        else:
            K2val.append(trunc_samples(K2, stdK, K_min_clip, K2_max_clip, (1, 1))[0][0])
            I21val.append(trunc_samples(0, stdI, -1, 1, (1, 1))[0][0])
    return(K1val, K2val, I12val, I21val)

def update_microbes(K1val, K2val, I12val, I21val, d, w, m, env_rat1, K1, K2, stdK, stdI, K_min_clip, K1_max_clip,
                    K2_max_clip):

    for h in range(len(K1val)):

        # random death
        K1val[h], I12val[h] = death(K1val[h], I12val[h], d)
        K2val[h], I21val[h] = death(K2val[h], I21val[h], d)

        prop1 = birth_propensity(K1val[h], w)
        prop2 = birth_propensity(K2val[h], w)

        interaction1 = (1-w)*sum(I21val[h])/len(I21val[h])
        interaction2 = (1-w)*sum(I12val[h])/len(I12val[h])

        prop1 += interaction1
        prop2 += interaction2

        # death by antagonistic interactions
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

        for i in selected_birth1:
            K1val[h].append(K1val[h][i])
            I12val[h].append(I12val[h][i])
        for i in selected_birth2:
            K2val[h].append(K2val[h][i])
            I21val[h].append(I21val[h][i])

        # colonization
        K1val[h], K2val[h], I12val[h], I21val[h] = colonize(K1val[h], K2val[h], I12val[h], I21val[h], env_rat1,
                                                            m, K1, K2, stdK, stdI, K_min_clip, K1_max_clip, K2_max_clip)

    return(K1val, K2val, I12val, I21val)
