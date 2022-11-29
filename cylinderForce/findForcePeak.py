#Find the peak moment of lift coefficient
#Copyright (c) 2022.11.29, Leo Yang.

import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

def findClPeak(clFilePath,heightRange,pointDis):
    data = np.loadtxt(clFilePath,delimiter=',')
    time = data[:,0]
    cl = data[:,2]

    [peakId,_] = signal.find_peaks(cl, height=heightRange, distance=pointDis)
    peakTime = time[peakId]
    peakCl = cl[peakId]

    return (time,cl,peakTime,peakCl)

def writeClPeak(peakFilePath,peakTime,peakCl):

    peak = np.vstack((peakTime,peakCl))
    np.savetxt(peakFilePath,peak.T,delimiter=',',fmt='%.7f')

if __name__ == '__main__':
    pwdRoot = os.getcwd()
    dirPath = os.path.join(pwdRoot, '072/forces')
    clFilePath = os.path.join(dirPath, 'cdxyz.csv')
    peakFilePath = os.path.join(dirPath, 'clPeak.csv')

    [time,cl,peakTime,peakCl] = findClPeak(clFilePath,-0.05,2000)

    writeClPeak(peakFilePath,peakTime,peakCl)

    plt.figure()
    plt.plot(time,cl)
    plt.scatter(peakTime,peakCl,marker='o',c='r',edgecolors='r')
    plt.show()