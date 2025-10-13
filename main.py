import csv
import numpy as np
import matplotlib.pyplot as plt

result1 = np.empty(shape=(102, 6))
result2 = np.empty(shape=(102, 6))
result3 = np.empty(shape=(102, 6))
result4 = np.empty(shape=(102, 6))

with open('Build\Debug\RungeKutta.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result1[j][i] = cell
            i += 1
        j += 1

with open('Build\Debug\Adams_Bashforth_Moulton.csv', newline='') as csvfile2:
    reader = csv.reader(csvfile2, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result2[j][i] = cell
            i += 1
        j += 1

with open('Build\Debug\Milne_Simpson.csv', newline='') as csvfile3:
    reader = csv.reader(csvfile3, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result3[j][i] = cell
            i += 1
        j += 1

with open('Build\Debug\Analytical_solution.csv', newline='') as csvfile4:
    reader = csv.reader(csvfile4, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result4[j][i] = cell
            i += 1
        j += 1

t = np.arange(0., 10.2, 0.1)
plt.plot(t, result1[:, 3], 'r', label='X1_RungeKutta')
plt.plot(t, result1[:, 4], 'b', label='X2_RungeKutta')
plt.plot(t, result1[:, 5], 'g', label='X3_RungeKutta')

plt.plot(t, result2[:, 3], 'r', label='X1_Adams_Bashforth_Moulton')
plt.plot(t, result2[:, 4], 'b', label='X2_Adams_Bashforth_Moulton')
plt.plot(t, result2[:, 5], 'g', label='X3_Adams_Bashforth_Moulton')

plt.plot(t, result3[:, 3], 'r', label='X1_Milne_Simpson')
plt.plot(t, result3[:, 4], 'b', label='X2_Milne_Simpson')
plt.plot(t, result3[:, 5], 'g', label='X3_Milne_Simpson')

plt.plot(t, result4[:, 3], 'r', label='X1_Analytical')
plt.plot(t, result4[:, 4], 'b', label='X2_Analytical')
plt.plot(t, result4[:, 5], 'g', label='X3_Analytical')

plt.legend()
plt.show()
