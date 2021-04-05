import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

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

    a.plot(x, meanline, linewidth=2, color = col)
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
plot_track_mean_std('Data/K1.csv', 'b', a[0])
plot_track_mean_std('Data/K2.csv', 'r', a[0])

a[0].set_xlabel('Time')
a[0].set_ylabel('Fitness value')
a[0].set_title('Within-host evolution of K value distribution \n K1 = 500000, K2 = 100000, K1max = 1000000, K2max = 200000 \n K10/K20 = 0.2, m = 100')

a[1].grid()
plot_track_mean_std('Data/I1.csv', 'b', a[1])
plot_track_mean_std('Data/I2.csv', 'r', a[1])
a[1].set_xlabel('Time')
a[1].set_ylabel('Alpha value')
a[1].set_title('Within-host evolution of Alpha value distribution')

a[2].grid()
plot_track_popsize('Data/I1.csv', 'b', a[2])
plot_track_popsize('Data/I2.csv', 'r', a[2])
a[2].set_xlabel('Time')
a[2].set_ylabel('Population size')
a[2].set_title('Within-host evolution of population sizes')

plt.show()
