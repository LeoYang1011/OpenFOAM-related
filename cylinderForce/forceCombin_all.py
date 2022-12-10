#Extract the data output by the forces function in OpenFOAM and store it in the .dat file.
#If there are multiple time folders, the data will be automatically merged and sorted by time.
#Can be used to process force data or forceBin data.
#Copyright (c) 2022.11.29, Leo Yang.

import os
import numpy as np
import linecache

def tryReadData(forceFile):
    readed = False
    data = list()
    while not readed:
       try:
            data = np.loadtxt(forceFile, comments='#')
            readed = True
       except ValueError:
            os.system(r"sed -i '$d' %s" % (forceFile))

    return (data)

def writeData(dataSavePath,forceFile,timeDir):
    forceData = list()

    if timeDir:
        forceData = tryReadData(forceFile[0])

        if len(timeDir) > 1:
            time = forceData[:,0]
            forceData = forceData[np.where(time<=timeDir[1])]

            for i in range(len(timeDir)):
                file = forceFile[i]

                if i > 0 and i < len(timeDir)-1:
                    data = tryReadData(file)
                    time = data[:,0]
                    forceData = np.vstack((forceData,data[np.where(time<=timeDir[i+1])]))
                elif i == len(timeDir)-1:
                    data = tryReadData(file)
                    forceData = np.vstack((forceData,data))

    np.savetxt(dataSavePath,forceData,delimiter="\t",fmt='%.6f')

def readAndWriteTitle(titleSavePath,forceFile,titleNum):
    lines = linecache.getlines(forceFile)[0:titleNum]

    with open(titleSavePath,'w') as f:
        for line in lines:
            f.write(line)

def findForceFile(dirPath,bin):
    forceFile = list()
    timeDir = list()

    for dir in os.listdir(dirPath):
        currentPath = os.path.join(dirPath,dir)
        if os.path.isdir(currentPath):
            forceEachDir = list()
            for file in os.listdir(currentPath):
                if "force" in file:
                    if bin:
                        if "Bin" in file:
                            forceEachDir.append(file)
                    else:
                        if "Bin" not in file:
                            forceEachDir.append(file)

            if forceEachDir:
                timeDir.append(float(dir))
                if len(forceEachDir) > 1:
                    fileSize = list()
                    for i in range(len(forceEachDir)):
                        fileSize.append(os.path.getsize(os.path.join(currentPath,forceEachDir[i])))
                    forceFile.append(os.path.join(currentPath,forceEachDir[fileSize.index(max(fileSize))]))
                else:
                    forceFile.append(os.path.join(currentPath,forceEachDir[0]))

    forceFile = sorted(forceFile,key=lambda x: timeDir[forceFile.index(x)])
    timeDir.sort()

    for file in forceFile:
        print(file)

    return (forceFile,timeDir)

if __name__ == '__main__':
    pwdRoot = os.getcwd()

    bin = True
    dirPath = os.path.join(pwdRoot, '012/forces')

    if bin:
        titleNum = 11
        forceFileName = "forceBin.dat"
        titleFileName = "forceBinTitle.dat"
    else:
        titleNum = 4
        forceFileName = "force.dat"
        titleFileName = "forceTitle.dat"

    dataSavePath = os.path.join(dirPath,forceFileName)
    titleSavePath = os.path.join(dirPath,titleFileName)

    [forceFile, timeDir] = findForceFile(dirPath,bin)

    for file in forceFile:
        os.system(r'sed -i "s/(/ /g;s/)/ /g" %s' % (file))

    readAndWriteTitle(titleSavePath,forceFile[0],titleNum)
    writeData(dataSavePath,forceFile,timeDir)