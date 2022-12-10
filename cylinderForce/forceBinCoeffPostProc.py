#Further process the force coefficients of each Bin, including de-averaging,
#interpolating to the face center of each Bin and outputting a file in Tecpolt format
#Copyright (c) 2022.11.29, Leo Yang.

import os
import numpy as np
import math

def writeBinData2Tecplot(binCenter, geoMin, geoMax, time, cdxBin, cdyBin, cdzBin, tecName):

    if dataTFace:
        [row,col] = binCenter.shape

        notNan = 0
        for i in range(row):
            for j in range(col):
                if not np.isnan(binCenter[i,j]):
                    notNan += 1

        with open(tecName,'w') as f:
            f.write('TITLE     = "Force bin"\n')
            f.write('Variables = "time","Z/D","cdx","cdy","cdz"\n\n')
            f.write('ZONE I=' + str(len(time)) + ', J=' + str(notNan) + ' F=POINT\n')
            f.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )\n')

            for k in range(row):
                for j in range(col):
                    if not math.isnan(binCenter[k,j]):
                        locatej = '{:.9E}'.format(binCenter[k,j])
                        for i in range(len(time)):
                            timei = '{:.9E}'.format(time[i])
                            cdxik = '{:.9E}'.format(cdxBin[i,k])
                            cdyik = '{:.9E}'.format(cdyBin[i,k])
                            cdzik = '{:.9E}'.format(cdzBin[i,k])
                            f.write(str(timei) + ' ' + str(locatej) + ' ')
                            f.write(str(cdxik) + ' ' + str(cdyik) + ' ' + str(cdzik))
                            f.write('\n')
    else:
         binCenter = np.append(geoMin,binCenter)
         binCenter = np.append(binCenter,geoMax)
         cdxBinNew = np.zeros(shape=(cdzBin.shape[0],len(binCenter)))
         cdyBinNew = np.zeros(shape=(cdzBin.shape[0],len(binCenter)))
         cdzBinNew = np.zeros(shape=(cdzBin.shape[0],len(binCenter)))
         cdxBinNew[:,0] = cdxBin[:,0]
         cdxBinNew[:,1:-1] = cdxBin
         cdxBinNew[:,-1] = cdxBin[:,-1]
         cdyBinNew[:,0] = cdyBin[:,0]
         cdyBinNew[:,1:-1] = cdyBin
         cdyBinNew[:,-1] = cdyBin[:,-1]
         cdzBinNew[:,0] = cdzBin[:,0]
         cdzBinNew[:,1:-1] = cdzBin
         cdzBinNew[:,-1] = cdzBin[:,-1]

         with open(tecName,'w') as f:
            f.write('TITLE     = "Force bin"\n')
            f.write('Variables = "time","Z/D","cdx","cdy","cdz"\n\n')
            f.write('ZONE I=' + str(len(time)) + ', J=' + str(len(binCenter)) + ' F=POINT\n')
            f.write('DT=(SINGLE SINGLE SINGLE SINGLE SINGLE )\n')

            for j in range(len(binCenter)):
                locatej = '{:.9E}'.format(binCenter[j])
                for i in range(len(time)):
                    timei = '{:.9E}'.format(time[i])
                    cdxik = '{:.9E}'.format(cdxBinNew[i,j])
                    cdyik = '{:.9E}'.format(cdyBinNew[i,j])
                    cdzik = '{:.9E}'.format(cdzBinNew[i,j])
                    f.write(str(timei) + ' ' + str(locatej) + ' ')
                    f.write(str(cdxik) + ' ' + str(cdyik) + ' ' + str(cdzik))
                    f.write('\n')


def deAverage(data,mean):
    dataShape = data.shape
    dataNormal = np.zeros(shape=dataShape)

    for i in range(dataShape[1]):
        dataNormal[:,i] = data[:,i] - mean[i]

    return dataNormal

def readData(meanName,cdxName,cdyName,cdzName,binLocateName):

    meanData = np.loadtxt(meanName,skiprows=1,delimiter=',')
    cdxMeanData = meanData[:,1]
    cdyMeanData = meanData[:,2]
    cdzMeanData = meanData[:,3]

    cdxData = np.loadtxt(cdxName,skiprows=1,delimiter=',')
    time = cdxData[:,0]
    cdxData = cdxData[:,1:]

    cdyData = np.loadtxt(cdyName,skiprows=1,delimiter=',')
    cdyData = cdyData[:,1:]

    cdzData = np.loadtxt(cdzName,skiprows=1,delimiter=',')
    cdzData = cdzData[:,1:]

    binCenter = meanData[:,0]

    if dataTFace:
        centerList = list()
        with open(binLocateName, "r") as f:
            for line in f.readlines():
                li = line.lstrip()
                if not li.startswith("#"):
                    ls = line.rstrip('\n')
                    lineSp = ls.split('\t')[1:]
                    centerList.append(lineSp)

        row = len(centerList)
        colume = -1
        for each in centerList:
            colume = max(len(each),colume)

        binCenter = np.full((row,colume),np.nan)
        for i in range(row):
            for j in range(len(centerList[i])):
                binCenter[i,j] = float(centerList[i][j])

    return (cdxMeanData,cdyMeanData,cdzMeanData,cdxData,cdyData,cdzData,binCenter,time)

if __name__ == '__main__':
    pwd_root = os.getcwd()

    dataTFace = False
    geoMin = 0
    geoMax = 4
    dirPath = os.path.join(pwd_root, '012/forces')

    binLocateName = os.path.join(dirPath, 'forceBinLocate')
    cdxName = os.path.join(dirPath, 'cdxBin.csv')
    cdyName = os.path.join(dirPath, 'cdyBin.csv')
    cdzName = os.path.join(dirPath, 'cdzBin.csv')
    meanName = os.path.join(dirPath, 'meanBin.csv')
    tecName = os.path.join(dirPath, 'cdxyzBinPost.dat')

    [cdxMean,cdyMean,cdzMean,cdx,cdy,cdz,binCenter,time] = readData(meanName,cdxName,cdyName,cdzName,binLocateName)

    cdy = deAverage(cdy,cdyMean)

    writeBinData2Tecplot(binCenter, geoMin, geoMax,time, cdx, cdy, cdz, tecName)
