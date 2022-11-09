import csv
import os
import numpy as np
import linecache

def writeFile(binLocate, time, cdxBin, cdyBin, cdzBin, cdxMeanBin, cdyMeanBin, cdzMeanBin,cdxName,cdyName,cdzName,meanName):
    meanTitle = ["height(m)","cdxMean","cdyMean","cdzMean"]
    with open(meanName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(meanTitle)
        for i in range(len(binLocate)):
            meanData = [binLocate[i],str(cdxMeanBin[i]),str(cdyMeanBin[i]),str(cdzMeanBin[i])]
            writer.writerow(meanData)

    first = ["height(m)"," "]
    first.extend(list(time))

    with open(cdxName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(binLocate)):
            row = [binLocate[i]," "]
            row.extend(list(cdxBin[i,:]))
            writer.writerow(row)

    with open(cdyName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(binLocate)):
            row = [binLocate[i]," "]
            row.extend(list(cdyBin[i,:]))
            writer.writerow(row)

    with open(cdzName, "w", newline='') as f:
        writer = csv.writer(f)
        writer.writerow(first)
        for i in range(len(binLocate)):
            row = [binLocate[i]," "]
            row.extend(list(cdzBin[i,:]))
            writer.writerow(row)

def calcCDxyz(binData,UInf,SpLength,D,binNum):
    cdx = list()
    cdy = list()
    cdz = list()

    forceX = binData[:,0]
    forceY = binData[:,1]
    forceZ = binData[:,2]

    for i in range(len(forceX)):
        cdxi = -2 * forceX[i] / (UInf ** 2 * SpLength/binNum * D * 1000)
        cdyi =  2 * forceY[i] / (UInf ** 2 * SpLength/binNum * D * 1000)
        cdzi =  2 * forceZ[i] / (UInf ** 2 * SpLength/binNum * D * 1000)
        cdx.append(cdxi)
        cdy.append(cdyi)
        cdz.append(cdzi)

    return (cdx,cdy,cdz)

def readBinLocate(inputName):

    binLocate = linecache.getline(inputName,9)
    binLocate = binLocate.rstrip('\n')
    binLocate = binLocate.split("\t")
    del(binLocate[0])

    return (binLocate)

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

if __name__ == '__main__':
    pwd_root = os.getcwd()
    D = 2.3         #圆柱直径
    UInf = 0.4     #入流速度
    SpLength = 10   #圆柱长度
    startTime = 40
    endTime = 140
    binNum = 20

    dir_path = os.path.join(pwd_root, 'forces36\\23')
    inputName = os.path.join(dir_path, 'forceBin.dat')
    cdxName = os.path.join(dir_path, 'cdx.csv')
    cdyName = os.path.join(dir_path, 'cdy.csv')
    cdzName = os.path.join(dir_path, 'cdz.csv')
    meanName = os.path.join(dir_path, 'mean.csv')

    binLocate = readBinLocate(inputName)
    [time,data] = readData(inputName,startTime,endTime)

    cdxBin = np.zeros((binNum,len(time)))
    cdyBin = np.zeros((binNum,len(time)))
    cdzBin = np.zeros((binNum,len(time)))
    cdxMeanBin = np.zeros(binNum)
    cdyMeanBin = np.zeros(binNum)
    cdzMeanBin = np.zeros(binNum)

    for i in range(len(binLocate)):
        j = 3*i
        binData = data[:,j:j+3]
        [cdx, cdy, cdz] = calcCDxyz(binData, UInf, SpLength, D, binNum)
        cdxBin[i,:] = cdx
        cdyBin[i,:] = cdy
        cdzBin[i,:] = cdz
        cdxMeanBin[i] = np.mean(cdx)
        cdyMeanBin[i] = np.mean(cdy)
        cdzMeanBin[i] = np.mean(cdz)

    writeFile(binLocate, time, cdxBin, cdyBin, cdzBin, cdxMeanBin, cdyMeanBin, cdzMeanBin, cdxName, cdyName, cdzName,meanName)
