import numpy as np
import random
import collections
import pandas as pd
import matplotlib.pyplot as plt

Parameters = {'K1' : 10000,
              'K2': 1000,
              'a12' : 0,
              'a21' : 0,
              'm' : 1,
              'w' : 0.5,
              'b' : 0.01,
              'T' : 1000,
              'd' : 0.1,
              'K01' : 10000,
              'K02' : 100000,
              'M' : 100}

sim_time = 10000

w = Parameters['w']
K1 = Parameters['K1']
K2 = Parameters['K2']
a12 = Parameters['a12']
a21 = Parameters['a21']
env_ratio = Parameters['K01']/(Parameters['K02']+Parameters['K01'])

A = np.ones((Parameters['M'],2))

data = []

for t in range(sim_time):

    dfrow = []

    for i in range(Parameters['M']):
        probs = np.random.uniform(0,1,6)

        Hdprob = probs[0]
        if Hdprob < 1/Parameters['T']:
            mics = int(A[i][0])*[1] + int(A[i][1])*[2]
            offspring = random.sample(mics,int(round(Parameters['b']*len(mics),0)))
            num = collections.Counter(offspring)
            A[i][0] = num[1]
            A[i][1] = num[2]
            continue

        d1 = probs[1]
        if d1 < Parameters['d']:
            A[i][0] = max(0,A[i][0]-1)
        d2 = probs[2]
        if d2 < Parameters['d']:
            A[i][1] = max(0, A[i][1] - 1)

        b1 = w*(1 - A[i][0]/K1)
        I1 = (1-w)*(a21*A[i][1]/(A[i][0]+A[i][1])) if (A[i][0]+A[i][1]) else 0

        b2 = w * (1 - A[i][1] / K2)
        I2 = (1 - w) * (a12 * A[i][0] / (A[i][0] + A[i][1])) if (A[i][0]+A[i][1]) else 0

        if 0 < b1 + I1 <= probs[3]:
            A[i][0] += 1

        if 0 < b2 + I2 <= probs[4]:
            A[i][1] += 1

        if probs[5]<Parameters['m']:
            p = random.random()
            if p<env_ratio:
                A[i][0] += 1
            else:
                A[i][1] += 1

        dfrow.append(A[i][0] / (A[i][0] + A[i][1]) if (A[i][0]+A[i][1]) else 0)

    data.append(dfrow)
    print(t)

df = pd.DataFrame(data)
print(df)

plt.plot(df, alpha=0.05)
plt.show()
