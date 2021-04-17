import pandas as pd
import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

file = 'Data_manyhost_sel2_nom.csv'

with open(file) as csvfile:
    rows = csv.reader(csvfile)
    dat = list(rows)

H = 15  # number of hosts
cols = ['Par', 'Time'] + ['H'+str(i) for i in range(1,H+1,1)]
df = pd.DataFrame(dat, columns=cols)

#print(df[[cols[1:]]] )
for i in cols[1:]:
    df[i] = pd.to_numeric(df[i])
print(df)

dfK = df.loc[df['Par'].isin(['K1','K2'])]

fig, a = plt.subplots(1,3)

a[0].grid()
lines = [Line2D([0], [0], color='blue'), Line2D([0], [0], color='red')]
labels = ['Microbe 1', 'Microbe 2']
for i in range(H):
    sns.lineplot(data=dfK, x='Time', y='H'+str(i+1), hue='Par', alpha=0.3, ax=a[0], palette=['b','r'])
a[0].legend(lines, labels)
a[0].set_xlabel('Time')
a[0].set_ylabel('Fitness value')
a[0].set_title('Evolution of mean K value in 15 replicates\n K1 = 5000, K2 = 1000, stdK1 = 100, stdK2 = 1000 \n K10/K20 = 0.1, m = 100, b = 0.05, T = 500')

dfI = df.loc[df['Par'].isin(['I1','I2'])]

a[1].grid()
lines = [Line2D([0], [0], color='blue'), Line2D([0], [0], color='red')]
labels = ['Microbe 1', 'Microbe 2']
for i in range(H):
    sns.lineplot(data=dfI, x='Time', y='H'+str(i+1), hue='Par', alpha=0.3, ax=a[1], palette=['b','r'])
a[1].legend(lines, labels)
a[1].set_xlabel('Time')
a[1].set_ylabel('Alpha value')
a[1].set_title('Evolution of mean Alpha values in 15 replicates')

dfM = df.loc[df['Par'].isin(['M1','M2'])]

a[2].grid()
lines = [Line2D([0], [0], color='blue'), Line2D([0], [0], color='red')]
labels = ['Microbe 1', 'Microbe 2']
for i in range(H):
    sns.lineplot(data=dfM, x='Time', y='H'+str(i+1), hue='Par', alpha=0.3, ax=a[2], palette=['b','r'])
a[2].legend(lines, labels)
a[2].set_xlabel('Time')
a[2].set_ylabel('Population size')
a[2].set_title('Evolution of Population sizes in 15 replicates')

plt.show()