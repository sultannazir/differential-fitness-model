import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns

file = 'Data_manyhost.csv'

with open(file) as csvfile:
    rows = csv.reader(csvfile)
    dat = list(rows)

H = 100  # number of hosts
cols = ['Par', 'Time'] + ['H'+str(i) for i in range(1,H+1,1)]
df = pd.DataFrame(dat, columns=cols)

#print(df[[cols[1:]]] )
for i in cols[1:]:
    df[i] = pd.to_numeric(df[i])
print(df)

df1 = df.loc[df['Par'].isin(['K1','K2','M1','M2'])]
for i in range(H):
    sns.lineplot(data=df1, x='Time', y='H'+str(i+1), hue='Par', alpha=0.2)
plt.show()