import seaborn as sns
import matplotlib.pyplot as plt
import csv
import random
import pandas as pd
import numpy as np

with open('Data/K1dist.csv') as csvfile:
    rows = csv.reader(csvfile)
    dat1 = list(rows)

newdat1 = []
for i in range(len(dat1)):
    if i%10 == 0:
        newdat1 += random.sample(dat1[i], 100)

#newdat1 = [float(num) for num in newdat1]
#sns.kdeplot(newdat1)
#plt.show()

#g = np.array(sum([[str(i)]*len(dat1[10*i]) for i in range(0,10)], []))
g = np.array(sum([[str(i)]*100 for i in range(0,10)], []))
print(len(g), len(newdat1))

df = pd.DataFrame(dict(x=newdat1, g=g))

print(df)
# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
g = sns.FacetGrid(df, row="g", hue="g", aspect=15, height=.5, palette=pal)
df["x"] = pd.to_numeric(df["x"], downcast="float")

#g.map(sns.histplot, "x")
g.map(sns.kdeplot, "x",
      bw_adjust=.5, clip_on=False,
      fill=True, alpha=1, linewidth=1.5)
g.map(plt.axhline, y=0, lw=2, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


g.map(label, "x")

# Set the subplots to overlap
g.fig.subplots_adjust(hspace=-0)

# Remove axes details that don't play well with overlap
g.set_titles("")
g.set(yticks=[])
g.despine(bottom=True, left=True)
plt.tight_layout()

plt.show()

