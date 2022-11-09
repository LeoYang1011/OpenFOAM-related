import os
import numpy as np

def readUV(file_path,startTime,endTime):
    totalUV = np.loadtxt(file_path)
    time = totalUV[:,0]
    if(np.max(time) < startTime or np.min(time) > endTime):
        empty = []
        return np.array(empty)
    else:
        reUV = np.zeros(shape=totalUV.shape)
        j = 0
        for i in range(len(time)):
            if (time[i] >= startTime and time[i] <= endTime):
                reUV[j,:] = totalUV[i,:]
                j += 1
        return reUV

if __name__ == '__main__':
    pwd_root = os.getcwd()
    pwd_root = os.path.join(pwd_root, 'probes')
    startTime = 100
    endTime = 200
    deltaT = 0.001
    probeN = 880
    TimeUV = np.zeros(shape=(int((endTime-startTime)/deltaT),probeN*3+1))
    j = 0
    times = ['0']

    for timedir in times:
        dir_path = os.path.join(pwd_root, timedir)
        if os.path.isdir(dir_path):
            file_path = os.path.join(dir_path, "U")
            print(file_path)
            fileUV = readUV(file_path,startTime,endTime)
            if (fileUV.size != 0):
                time = fileUV[:,0]
                for i in range(len(time)):
                    if time[i] not in TimeUV[:,0]:
                        TimeUV[j,:] = fileUV[i,:]
                        j += 1
                    if j >= int((endTime-startTime)/deltaT):
                        break

    np.savetxt("timeUV.csv", TimeUV, delimiter=',')



