import csv
import numpy as np
import matplotlib.pyplot as plt

result1 = np.empty(shape=(102, 6))
result2 = np.empty(shape=(102, 6))
result3 = np.empty(shape=(102, 6))
result1_ = np.empty(shape=(102, 6))
result2_ = np.empty(shape=(102, 6))
result3_ = np.empty(shape=(102, 6))
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

# with open('Build\Debug\Analytical_solution.csv', newline='') as csvfile4:
#     reader = csv.reader(csvfile4, delimiter=',')
#     j = 0
#     for row in reader:
#         i = 0
#         for cell in row:
#             result4[j][i] = cell
#             i += 1
#         j += 1

with open('Build\Debug\RungeKutta+.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result1_[j][i] = cell
            i += 1
        j += 1

with open('Build\Debug\Adams_Bashforth_Moulton+.csv', newline='') as csvfile2:
    reader = csv.reader(csvfile2, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result2_[j][i] = cell
            i += 1
        j += 1

with open('Build\Debug\Milne_Simpson+.csv', newline='') as csvfile3:
    reader = csv.reader(csvfile3, delimiter=',')
    j = 0
    for row in reader:
        i = 0
        for cell in row:
            result3_[j][i] = cell
            i += 1
        j += 1

t = np.arange(0., 10.2, 0.1)

#________________Рунге-Кутта___________________
# plt.plot(t, result1[:, 3], 'g', label='X1_M-')
# # plt.plot(t, result1_[:, 3], 'r', label='X1_Runge_Kutta')
# plt.plot(t, result1[:, 4], 'g', label='X2_M-')
# # plt.plot(t, result1_[:, 4], 'r', label='X2_Runge_Kutta')
# plt.plot(t, result1[:, 5], 'g', label='X3_M-')
# # plt.plot(t, result1_[:, 5], 'r', label='X3_Runge_Kutta')

# plt.plot(t, result1[:, 0], 'r', label='Y1_M-')
# plt.plot(t, result1[:, 1], 'g', label='Y2_M-')
# plt.plot(t, result1[:, 2], 'b', label='Y3_M-')
plt.plot(t, result1[:, 3], '#FF5733', label='X1_M-')
plt.plot(t, result1[:, 4], '#FF57FF', label='X2_M-')
plt.plot(t, result1[:, 5], '#57FFFF', label='X3_M-')

# plt.plot(t, result1_[:, 0], 'r', label='Y1_M+')
# plt.plot(t, result1_[:, 1], 'g', label='Y2_M+')
# plt.plot(t, result1_[:, 2], 'b', label='Y3_M+')
# plt.plot(t, result1_[:, 3], '#FF5733', label='X1_M+')
# plt.plot(t, result1_[:, 4], '#FF57FF', label='X2_M+')
# plt.plot(t, result1_[:, 5], '#57FFFF', label='X3_M+')

#__________ Адамс Башфорт Моултон___________
# plt.plot(t, result2_[:, 3], 'r', label='X1_Adams_Bashforth_Moulton')
plt.plot(t, result2[:, 3], 'b', label='X2_Adams_Bashforth_Moulton')
# plt.plot(t, result2_[:, 4], 'g', label='X3_Adams_Bashforth_Moulton')
plt.plot(t, result2[:, 4], 'r', label='X1_Adams_Bashforth_Moulton')
# plt.plot(t, result2_[:, 5], 'b', label='X2_Adams_Bashforth_Moulton')
plt.plot(t, result2[:, 5], 'g', label='X3_Adams_Bashforth_Moulton')


#_______________Милн-Симпсон__________________
# plt.plot(t, result3_[:, 3], 'r', label='X1_Milne_Simpson')
# plt.plot(t, result3_[:, 4], 'b', label='X2_Milne_Simpson')
# plt.plot(t, result3_[:, 5], 'g', label='X3_Milne_Simpson')

plt.plot(t, result3[:, 3], 'r', label='X1_Milne_Simpson')
plt.plot(t, result3[:, 4], 'b', label='X2_Milne_Simpson')
plt.plot(t, result3[:, 5], 'g', label='X3_Milne_Simpson')

# plt.plot(t, result4[:, 3], 'r', label='X1_Analytical')
# plt.plot(t, result4[:, 4], 'b', label='X2_Analytical')
# plt.plot(t, result4[:, 5], 'g', label='X3_Analytical')

plt.legend()
plt.show()
