#Calculate the force coefficient of each Bin based on the forceBin data.
#The result of the force coefficient in each direction will be written into the corresponding csv file.
#Writing data to tecplot format is also optinal.
#The 'forceCombin_all.py' needs to be run before running this program.
#The file 'forceBinLocate' stores the length of each Bin and the location of the face center contained in the Bin.
#This file can be obtained by running the post-processing program 'postForceBinLength'
#Copyright (c) 2022.11.29, Leo Yang.

import csv
import os
import re
import numpy as np
import linecache
from scipy import signal
import matplotlib.pyplot as plt

def writeBinData2Tecplot(binLocate, time, cdxBin, cdyBin, cdzBin, tecName):

    with open(tecName,'w') as f:
        f.write('TITLE     = "Force bin"\n')
        f.write('Variables = "time","Z/D","cdx","cdy","cdz"\n\n')
        f.write('ZONE I=' + str(len(time)) + ', J=' + str(len(binLocate)) + ' F=POINT\n')
        f.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )\n')

        for j in range(len(binLocate)):
            locatej = '{:.9E}'.format(binLocate[j])
            for i in range(len(time)):
                timei = '{:.9E}'.format(time[i])
                cdxij = '{:.9E}'.format(cdxBin[i,j])
                cdyij = '{:.9E}'.format(cdyBin[i,j])
                cdzij = '{:.9E}'.format(cdzBin[i,j])
                f.write(str(timei) + ' ' + str(locatej) + ' ')
                f.write(str(cdxij) + ' ' + str(cdyij) + ' ' + str(cdzij))
                f.write('\n')

def writeBinData2CSV(binLocate, time, cdxBin, cdyBin, cdzBin, cdxMeanBin, cdyMeanBin, cdzMeanBin,cdxName,cdyName,cdzName,meanName):

    meanTitle = ["locate(m)","cdxMean","cdyMean","cdzMean"]
    with open(meanName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(meanTitle)
        for i in range(len(binLocate)):
            meanData = [binLocate[i],str(cdxMeanBin[i]),str(cdyMeanBin[i]),str(cdzMeanBin[i])]
            writer.writerow(meanData)

    first = [" "]
    first.extend(list(binLocate))

    with open(cdxName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(time)):
            row = [time[i]]
            row.extend(list(cdxBin[i,:]))
            writer.writerow(row)

    with open(cdyName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(time)):
            row = [time[i]]
            row.extend(list(cdyBin[i,:]))
            writer.writerow(row)

    with open(cdzName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(time)):
            row = [time[i]]
            row.extend(list(cdzBin[i,:]))
            writer.writerow(row)

def calcCDxyz(binData,UInf,length,D,rho):
    cdx = list()
    cdy = list()
    cdz = list()

    forceX = binData[:,0]
    forceY = binData[:,1]
    forceZ = binData[:,2]

    for i in range(len(forceX)):
        cdxi =  2 * forceX[i] / (UInf ** 2 * length * D * rho)
        cdyi =  2 * forceY[i] / (UInf ** 2 * length * D * rho)
        cdzi =  2 * forceZ[i] / (UInf ** 2 * length * D * rho)
        cdx.append(cdxi)
        cdy.append(cdyi)
        cdz.append(cdzi)

    return (cdx,cdy,cdz)

def calcBinCdxyz(binNum, data, UInf, binLength, D, rho):
    cdxBin = np.zeros((len(time),binNum))
    cdyBin = np.zeros((len(time),binNum))
    cdzBin = np.zeros((len(time),binNum))
    cdxMeanBin = np.zeros(binNum)
    cdyMeanBin = np.zeros(binNum)
    cdzMeanBin = np.zeros(binNum)

    for i in range(binNum):
        j = 9*i
        binData = data[:,j:j+3]
        [cdx, cdy, cdz] = calcCDxyz(binData, UInf, binLength[i], D, rho)
        cdxBin[:,i] = cdx
        cdyBin[:,i] = cdy
        cdzBin[:,i] = cdz
        cdxMeanBin[i] = np.mean(cdx)
        cdyMeanBin[i] = np.mean(cdy)
        cdzMeanBin[i] = np.mean(cdz)

    return (cdxBin,cdyBin,cdzBin,cdxMeanBin,cdyMeanBin,cdzMeanBin)

def readBinLocate(inputTitle,binLocateName):

    binNum = linecache.getline(inputTitle,2)
    binNum = re.findall(r'\d+',binNum)
    binNum = eval(binNum[0])
    print("bins   : "+ str(binNum))

    binStart = linecache.getline(inputTitle,3)
    binStart = re.findall(r'[\\+|-]?\d+\.\d+[Ee][\\+|-]\d+',binStart)
    binStart = eval(binStart[0])
    print("start  : " + str(binStart))

    binDelta = linecache.getline(inputTitle,5)
    binDelta = re.findall(r'[\\+|-]?\d+\.\d+[Ee][\\+|-]\d+',binDelta)
    binDelta = eval(binDelta[0])
    print("delta  : " + str(binDelta))

    binDir = linecache.getline(inputTitle,6)
    binDir = re.findall(r'[\\+|-]?\d+\.\d+[Ee][\\+|-]\d+',binDir)
    for i in range(len(binDir)):
        binDir[i] = eval(binDir[i])
    binDir = np.array(binDir)
    print("binDir : " + "( " + str(binDir[0]) + " " + str(binDir[1]) + " " + str(binDir[2])+ " )")

    binCenter = np.zeros(shape=(binNum,3))
    for i in range(binNum):
        binCenter[i] = (binStart + (i + 0.5)*binDelta)*binDir

    binLength = list()
    with open(binLocateName, "r") as f:
        for line in f.readlines():
            li = line.lstrip()
            if not li.startswith("#"):
                ls = line.rstrip('\n')
                first = ls.split('\t')[0]
                if first:
                    binLength.append(float(first))

    binLength = np.array(binLength)

    return (binNum,binCenter,binDir,binLength)

def readData(inputName,startTime,endTime):

    print(inputName)
    data = np.loadtxt(inputName,comments='#')

    time = data[:,0]
    data = data[:,1:]

    time = np.array(time)

    timeOld = time.copy()
    dataOld = data.copy()

    time = timeOld[np.where(timeOld >= startTime)]
    data = dataOld[np.where(timeOld >= startTime)]

    timeOld = time.copy()
    dataOld = data.copy()

    time = timeOld[np.where(timeOld <= endTime)]
    data = dataOld[np.where(timeOld <= endTime)]

    return (time,data)


def medfiltCDxyz(data,n):
    dataShape = data.shape
    filtData  = np.zeros(shape=dataShape)

    for i in range(dataShape[1]):
        filtData[:,i]  = signal.medfilt(data[:,i], n)

    return (filtData)

if __name__ == '__main__':
    pwd_root = os.getcwd()
    D = 1.0         #Cylinder diameter
    UInf = 1.0      #Flow velocity
    rho = 1.0       #Flow density
    startTime = 200.1
    endTime = 300
    writeTec = False

    dirPath = os.path.join(pwd_root, '012/forces')

    inputData = os.path.join(dirPath, 'forceBin.dat')
    inputTitle = os.path.join(dirPath, 'forceBinTitle.dat')
    binLocateName = os.path.join(dirPath, 'forceBinLocate')

    cdxName = os.path.join(dirPath, 'cdxBin.csv')
    cdyName = os.path.join(dirPath, 'cdyBin.csv')
    cdzName = os.path.join(dirPath, 'cdzBin.csv')
    meanName = os.path.join(dirPath, 'meanBin.csv')
    tecName = os.path.join(dirPath, 'cdxyzBinRaw.dat')

    [binNum,binCenter,binDir,binLength] = readBinLocate(inputTitle,binLocateName)
    [time,data] = readData(inputData,startTime,endTime)

    [cdxBin,cdyBin,cdzBin,cdxMeanBin,cdyMeanBin,cdzMeanBin] = calcBinCdxyz(binNum,data,UInf,binLength,D,rho)
    cdxBin = medfiltCDxyz(cdxBin,9)
    cdyBin = medfiltCDxyz(cdyBin,5)
    cdzBin = medfiltCDxyz(cdzBin,5)

    binLocate = np.dot(binCenter,binDir)
    writeBinData2CSV(binLocate, time, cdxBin, cdyBin, cdzBin, cdxMeanBin, cdyMeanBin, cdzMeanBin, cdxName, cdyName, cdzName,meanName)

    if writeTec:
        writeBinData2Tecplot(binLocate, time, cdxBin, cdyBin, cdzBin, tecName)

    # x,y = np.meshgrid(time,binLocate)
    # plt.figure()
    # plt.contourf(x,y,cdyBin.T)
    # plt.show()
