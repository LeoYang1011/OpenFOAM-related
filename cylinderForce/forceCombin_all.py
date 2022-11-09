import os
import numpy as np
import linecache

def readAndWriteForceData(dataSavePath,forceFile,timeDir):
    forceData = np.loadtxt(forceFile[0], comments='#')

    if len(timeDir) > 1:
        time = forceData[:,0]
        forceData = forceData[np.where(time<=timeDir[1])]

        for i in range(len(timeDir)):
            file = forceFile[i]
            data = np.loadtxt(file, comments='#')
            if i > 0 and i < len(timeDir)-1:
                time = data[:,0]
                forceData = np.vstack((forceData,data[np.where(time<=timeDir[i+1])]))
            elif i == len(timeDir)-1:
                forceData = np.vstack((forceData,data))

    np.savetxt(dataSavePath,forceData,delimiter="\t",fmt='%.9f')

def readAndWriteFileTitle(titleSavePath,forceFile,titleNum):
    lines = linecache.getlines(forceFile)[0:titleNum]

    with open(titleSavePath,'w') as f:
        for line in lines:
            f.write(line)

def findForceFile(dirPath,titleNum):
    forceFile = list()
    timeDir = list()

    for dir in os.listdir(dirPath):
        currentPath = os.path.join(dirPath,dir)
        if os.path.isdir(currentPath):
            timeDir.append(float(dir))
            forceEachDir = list()
            for file in os.listdir(currentPath):
                if "force" in file:
                    if titleNum == 4:
                        if "Bin" not in file:
                            forceEachDir.append(file)
                    elif titleNum > 4:
                        if "Bin" in file:
                            forceEachDir.append(file)

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
    titleNum = 4

    dirPath = os.path.join(pwdRoot, '062/forces')
    dataSavePath = os.path.join(dirPath,"force.dat")
    titleSavePath = os.path.join(dirPath,"title.dat")

    [forceFile, timeDir] = findForceFile(dirPath,titleNum)

    for file in forceFile:
        os.system(r'sed -i "s/(/ /g;s/)/ /g" %s' % (file))

    readAndWriteFileTitle(titleSavePath,forceFile[0],titleNum)
    readAndWriteForceData(dataSavePath,forceFile,timeDir)