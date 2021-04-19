import numpy as np
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.lines import Line2D

### plots csv files of data already written after simulations

fileK1 = 'Data/K1.csv'  # filename of csv file of K1 dist timeseries data
fileK2 = 'Data/K2.csv'  # filename of csv file of K2 dist timeseries data
fileI1 = 'Data/I1.csv'  # filename of csv file of alpha1 dist timeseries data
fileI2 = 'Data/I2.csv'  # filename of csv file of alpha2 dist timeseries data

binsK = np.linspace(0,10000, 21)
binsI = np.linspace(0,1,21)

def plot_dist(filename, a, bins):
    with open(filename) as csvfile:
        rows = csv.reader(csvfile)
        dat = list(rows)

    hist = []
    for i in dat:
        i = [float(j) for j in i]
        hist.append(np.log10(1 + np.histogram(i, bins)[0]).tolist())
    hist= np.flip(np.array(hist),1)
    print(hist)
    im = a.imshow(hist.T, cmap='Greys', aspect='auto')
    ticks = [0,4,8,12,16]
    a.set_yticks(ticks)
    a.set_yticklabels([np.flip(np.around(bins, decimals=2))[i] for i in ticks])
    divider = make_axes_locatable(a)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')

def plot_track_stats(filename, col, a):

    with open(filename) as csvfile:
        rows = csv.reader(csvfile)
        dat = list(rows)

    means = []
    medians = []
    stds = []

    for i in dat:
        i = [float(j) for j in i]
        means.append(np.mean(i))
        medians.append((np.median(i)))
        stds.append(np.std(i))

    meanline = np.array(means)
    overline = np.array(means) + np.array(stds)
    underline = np.array(means) - np.array(stds)
    medline = np.array(medians)

    x = np.array([i for i in range(len(meanline))])

    a.plot(x, meanline, linewidth=2, color = col, alpha = 0.5)
    a.plot(x, medline, linewidth=1, color = col, linestyle='dashed')
    a.fill_between(x, underline, overline, color = col, alpha = 0.1)

def plot_track_popsize(filename, col, a):
    with open(filename) as csvfile:
        rows = csv.reader(csvfile)
        dat = list(rows)

    M = []
    for i in dat:
        M.append(len(i))

    M = np.array(M)
    a.plot(M, linewidth=2, color=col)


fig, a = plt.subplots(2,2, figsize=(9,5))

#a[0].grid()
#lines = [Line2D([0], [0], color='black'), Line2D([0], [0], color='black', linestyle='--'), Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
#labels = ['Mean', 'Median', 'Microbe 1', 'Microbe 2']
#a[0].legend(lines, labels)
plot_dist(fileK1, a[0][0], binsK)

a[0][0].set_xlabel('Time')
a[0][0].set_ylabel('K1 value')
#a[0].set_title('Within-host evolution of K value distribution \n K1 = 5000, K2 = 1000, stdK = 500, muK = 10, muI = 0.01\n K10/K20 = 1/9, m = 500, b = 0.05, T = 500')

plot_dist(fileK2, a[0][1], binsK)

a[0][1].set_xlabel('Time')
a[0][1].set_ylabel('K2 value')

plot_dist(fileI1, a[1][0], binsI)

a[1][0].set_xlabel('Time')
a[1][0].set_ylabel('Alpha1 value')

plot_dist(fileI2, a[1][1], binsI)

a[1][1].set_xlabel('Time')
a[1][1].set_ylabel('Alpha2 value')
#a[1].grid()
#plot_track_mean_std(fileI1, 'b', a[1])
#plot_track_mean_std(fileI2, 'r', a[1])
#a[1].set_xlabel('Time')
#a[1].set_ylabel('Alpha value')
#a[1].set_title('Within-host evolution of Alpha value distribution')

#a[2].grid()
#plot_track_popsize(fileI1, 'b', a[2])
#plot_track_popsize(fileI2, 'r', a[2])
#a[2].set_xlabel('Time')
#a[2].set_ylabel('Population size')
#a[2].set_title('Within-host evolution of population sizes')

plt.suptitle('Frequency distribution of evolvable parameters in a single host (in log10) \n with vertical (T=200, b=0.05) and horizontal drift (m=500)'
             '\n K1 = 5000, K2 = 1000, stdK = 500, K1max = K2max = 10000, muK = 10, muI = 0.01, K10/K20 = 1/9')

plt.tight_layout()
plt.show()
