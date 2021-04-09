import seaborn as sns
import matplotlib.pyplot as plt
import csv
import random
import pandas as pd
import numpy as np

with open('Data/BD.csv') as csvfile:
    rows = csv.reader(csvfile)
    dat1 = list(rows)

dat = [[np.log10(int(i)) for i in j] for j in dat1]

col = ['death1', 'death2', 'antdeath1', 'antdeath2', 'birth1', 'birth2']
df = pd.DataFrame(dat, columns=col)

df.iloc[170:400, 2:6].plot()
plt.title('Number of births and deaths due to antagonistic interspecific interactions')
plt.ylabel('Number in log10')
plt.xlabel('Time step')
plt.grid()
plt.show()
print(df)
