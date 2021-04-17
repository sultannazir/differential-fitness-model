import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

### plots csv files of data already written after simulations

fileK1 = 'Data/K1.csv'  # filename of csv file of K1 dist timeseries data
fileK2 = 'Data/K2.csv'  # filename of csv file of K2 dist timeseries data
fileI1 = 'Data/I1.csv'  # filename of csv file of alpha1 dist timeseries data
fileI2 = 'Data/I2.csv'  # filename of csv file of alpha2 dist timeseries data

def plot_track_mean_std(filename, col, a):

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


fig, a = plt.subplots(1,3)

a[0].grid()
lines = [Line2D([0], [0], color='black'), Line2D([0], [0], color='black', linestyle='--'), Line2D([0], [0], color='b'), Line2D([0], [0], color='r')]
labels = ['Mean', 'Median', 'Microbe 1', 'Microbe 2']
a[0].legend(lines, labels)
plot_track_mean_std(fileK1, 'b', a[0])
plot_track_mean_std(fileK2, 'r', a[0])

a[0].set_xlabel('Time')
a[0].set_ylabel('Fitness value')
a[0].set_title('Within-host evolution of K value distribution \n K1 = 5000, K2 = 1000, stdK1 = 100, stdK2 = 1000 \n K10/K20 = 0.001, m = 1000, b = 0.05, T = 500')

a[1].grid()
plot_track_mean_std(fileI1, 'b', a[1])
plot_track_mean_std(fileI2, 'r', a[1])
a[1].set_xlabel('Time')
a[1].set_ylabel('Alpha value')
a[1].set_title('Within-host evolution of Alpha value distribution')

a[2].grid()
plot_track_popsize(fileI1, 'b', a[2])
plot_track_popsize(fileI2, 'r', a[2])
a[2].set_xlabel('Time')
a[2].set_ylabel('Population size')
a[2].set_title('Within-host evolution of population sizes')

plt.show()
