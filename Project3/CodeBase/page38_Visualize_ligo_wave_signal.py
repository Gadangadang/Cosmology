import numpy as np
import matplotlib.pyplot as plt

Observed_data = open("observed.txt", "r")

time_array_observed = []
measurment_observed = []

Observed_data.readline()

for line in Observed_data:
    time, data = line.split()
    time_array_observed.append(float(time))
    measurment_observed.append(float(data))

Numerical_data = open("numerical_sol.txt", "r")

time_array_numerical = []
numerical_data = []

Numerical_data.readline()

for line in Numerical_data:
    time, data = line.split()
    time_array_numerical.append(float(time))
    numerical_data.append(float(data))

plt.plot(time_array_observed, measurment_observed, label="Observed")
plt.plot(time_array_numerical, numerical_data, label="Numerical")
plt.xlabel("t [s]")
plt.ylabel("strain")
plt.title("Data multiplied with 1e21")
plt.legend()
plt.savefig("page38_straindata.jpeg")

plt.show()
