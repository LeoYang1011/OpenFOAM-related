import csv
import numpy as np

timeUV = np.loadtxt("timeUV.csv",delimiter=",")
UV = timeUV[:,1:]

m = UV.shape[0]
n = UV.shape[1]
meanUV = np.zeros(shape=(n,1))
for i in range(n):
    meanUV[i] = np.mean(UV[:,i])

UfVf = np.zeros(shape=UV.shape)
meanUUVV = np.zeros(shape=(n,1))
for i in range(n):
    UUVV = []
    for j in range(m):
        UfVf[j,i] = UV[j,i] - meanUV[i]
        UUVV.append(UfVf[j,i]**2)
    meanUUVV[i] = np.mean(np.array(UUVV))


np.savetxt("meanUV.csv", meanUV, delimiter=',')
np.savetxt("meanUUVV.csv", meanUUVV, delimiter=',')
np.savetxt("UfVf.csv", UfVf, delimiter=',')



