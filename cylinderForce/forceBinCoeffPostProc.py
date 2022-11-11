import os
import numpy as np
from scipy import interpolate
import linecache
import matplotlib.pyplot as plt
from scipy import signal

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

def dataNormalization(data,mean):
    dataShape = data.shape
    dataNormal = np.zeros(shape=dataShape)

    for i in range(dataShape[1]):
        dataNormal[:,i] = data[:,i] - mean[i]

    dataNormal[dataNormal < 0] = -1
    dataNormal[dataNormal > 0] = 1

    return dataNormal

def dataInterpSpan(data,time,binLocate):
    interpFunc = interpolate.interp2d(time, binLocate, data.T, kind='cubic')

    binLocateNew = np.linspace(min(binLocate),max(binLocate),100)
    dataNew = interpFunc(time,binLocateNew)

    return dataNew.T

def readData(meanName,cdxName,cdyName,cdzName):

    meanData = np.loadtxt(meanName,skiprows=1,delimiter=',')
    binLocate = meanData[:,0]
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

    return (cdxMeanData,cdyMeanData,cdzMeanData,cdxData,cdyData,cdzData,binLocate,time)

if __name__ == '__main__':
    pwd_root = os.getcwd()

    dirPath = os.path.join(pwd_root, '072/forces')

    cdxName = os.path.join(dirPath, 'cdxBin.csv')
    cdyName = os.path.join(dirPath, 'cdyBin.csv')
    cdzName = os.path.join(dirPath, 'cdzBin.csv')
    meanName = os.path.join(dirPath, 'meanBin.csv')
    tecName = os.path.join(dirPath, 'cdxyzBinPost.dat')

    [cdxMean,cdyMean,cdzMean,cdx,cdy,cdz,binLocate,time] = readData(meanName,cdxName,cdyName,cdzName)

    #spanwiseSmooth()
    # cdx = dataNormalization(cdx,cdxMean)
    # cdy = dataNormalization(cdy,cdyMean)
    # cdz = dataNormalization(cdz,cdzMean)

    cdy = dataInterpSpan(cdy,time,binLocate)

    writeBinData2Tecplot(binLocate, time, cdx, cdy, cdz, tecName)
