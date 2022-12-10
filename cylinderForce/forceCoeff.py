#Calculate the force coefficient based on the force data.
#The result of the force coefficient in each direction will be written into the corresponding csv file.
#The spectrum of the lift coefficient is also calculated and written to a csv file.
#The 'forceCombin_all.py' needs to be run before running this program.
#Copyright (c) 2022.11.29, Leo Yang.

import csv
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

def writeFile(time,cdx,cdy,cdz,st,cdxMean,cdyMean,cdzMean,cdyRMS,meanName,cdxyzName,powerName):
    eachMean = ["st","cdxMean","cdyMean","cdzMean","cdyRMS"]
    meanData = [str(st), str(cdxMean), str(cdyMean), str(cdzMean), str(cdyRMS)]
    rows = zip(eachMean, meanData)
    with open(meanName, "w", newline='') as f:
        writer = csv.writer(f)
        for row in rows:
            writer.writerow(row)

    cdxyz = np.vstack((time,cdx,cdy,cdz))
    power = np.vstack((frq,fftCdy))

    np.savetxt(cdxyzName,cdxyz.T,delimiter=',',fmt='%.9f')
    np.savetxt(powerName, power.T, delimiter=',')

def powerCl(time,cl,dt,fftStart,fftEnd):
    clOld = cl.copy()
    timeOld = time.copy()
    time = timeOld[np.where(timeOld >= fftStart)]
    cl = clOld[np.where(time >= fftStart)]

    clOld = cl.copy()
    timeOld = time.copy()
    cl = clOld[np.where(timeOld <= fftEnd)]

    n = len(cl)
    k = np.arange(n)
    T = n * dt
    frqFull = k / T
    frq = frqFull[range(int(n / 2))]

    fftClFull = np.fft.fft(cl - np.mean(cl))
    fftCl = np.abs(fftClFull)/(n / 2.0)
    fftCl = fftCl[range(int(n / 2))]

    frq = frq[2:]
    fftCl = fftCl[2:]

    clMean = np.mean(cl)
    clRMS = 0

    for cli in cl:
        clRMS += (cli - clMean) ** 2

    clRMS = clRMS / len(cl)
    clRMS = np.sqrt(clRMS)

    return (frq,fftCl,clRMS)

def medfiltCDxyz(cd,n):
    filtCd  = signal.medfilt(cd, n)
    # filtCd = cd
    # for i in range(len(cd)):
    #     if i > 0 and i < len(cd) - 1:
    #         mean = 0.5*(cd[i-1] + cd[i+1])
    #         rel = abs(mean - cd[i])/(1.0*abs(cd[i])+0.000001)
    #         print(rel)
    #         if rel > 0.1:
    #             filtCd[i] = mean

    return (filtCd)

def calcCDxyz(inputName,startTime,endTime,UInf,SpLength,D, rho):
    cdx = list()
    cdy = list()
    cdz = list()

    print(inputName)
    data = np.loadtxt(inputName,comments='#')

    time = data[:,0]
    totalx = data[:,1]
    totaly = data[:,2]
    totalz = data[:,3]

    # time = np.array(time)
    # totalx = np.array(totalx)
    # totaly = np.array(totaly)
    # totalz = np.array(totalz)

    timeOld = time.copy()
    totalxOld = totalx.copy()
    totalyOld = totaly.copy()
    totalzOld = totalz.copy()

    time = timeOld[np.where(timeOld >= startTime)]
    totalx = totalxOld[np.where(timeOld >= startTime)]
    totaly = totalyOld[np.where(timeOld >= startTime)]
    totalz = totalzOld[np.where(timeOld >= startTime)]

    timeOld = time.copy()
    totalxOld = totalx.copy()
    totalyOld = totaly.copy()
    totalzOld = totalz.copy()

    time = timeOld[np.where(timeOld <= endTime)]
    totalx = totalxOld[np.where(timeOld <= endTime)]
    totaly = totalyOld[np.where(timeOld <= endTime)]
    totalz = totalzOld[np.where(timeOld <= endTime)]

    for i in range(len(time)):
        cdxi =  2 * totalx[i] / (UInf ** 2 * SpLength * D * rho)
        cdyi =  2 * totaly[i] / (UInf ** 2 * SpLength * D * rho)
        cdzi =  2 * totalz[i] / (UInf ** 2 * SpLength * D * rho)
        cdx.append(cdxi)
        cdy.append(cdyi)
        cdz.append(cdzi)

    return (time,cdx,cdy,cdz)

if __name__ == '__main__':
    pwd_root = os.getcwd()
    D = 1.0              #Cylinder diameter
    UInf = 1.0           #Flow velocity
    SpLength = 4         #Cylinder length
    rho = 1.0            #Flow density
    startTime = 150
    endTime = 300
    dt = 0.002           #Delta time for fft
    fftStart = 156.618
    fftEnd = 294.54

    dir_path = os.path.join(pwd_root, '012/forces')
    inputName = os.path.join(dir_path, 'force.dat')
    cdxyzName = os.path.join(dir_path, 'cdxyz.csv')
    powerName = os.path.join(dir_path, 'power.csv')
    meanName = os.path.join(dir_path, 'mean.csv')

    [time,cdx,cdy,cdz] = calcCDxyz(inputName, startTime, endTime, UInf, SpLength, D, rho)
    cdxNew = medfiltCDxyz(cdx,5) #Numerical filtering to remove outlier data points (optional)
    cdyNew = medfiltCDxyz(cdy,5)
    [frq,fftCdy,cdyRMS] = powerCl(time,cdyNew,dt,fftStart,fftEnd)

    st = frq[np.argmax(fftCdy)]
    cdxMean = np.mean(cdx)
    cdyMean = np.mean(cdy)
    cdzMean = np.mean(cdz)

    writeFile(time, cdxNew, cdyNew, cdz, st, cdxMean, cdyMean, cdzMean, cdyRMS, meanName, cdxyzName, powerName)

    print("st: "+str(st))
    print("cdxMean: "+str(cdxMean))
    print("cdyMean: "+str(cdyMean))
    print("cdzMean: "+str(cdzMean))
    print("cdyRMS: "+str(cdyRMS))

    plt.figure()
    plt.plot(frq[:50], fftCdy[:50])

    plt.figure()
    plt.plot(time, cdxNew)
    plt.plot(time, cdyNew)
    plt.show()
