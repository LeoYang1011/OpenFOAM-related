#Copyright (c) 2023.04.16, Leo Yang.
import numpy as np
import matplotlib.pyplot as plt
import os

def writeData(writePath,real,imag,frq,k,nUsed,solver):
    fileId = open(writePath,'w')
    n = nUsed

    if(solver == "waveFoam"):
        omg = 2.0*np.pi*frq
        fileId.write("realPart           nonuniform List<scalar>\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            fileId.write("\t" + str(real[i]) + "\n")
        fileId.write(");\n")

        fileId.write("imagPart           nonuniform List<scalar>\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            fileId.write("\t" + str(imag[i]) + "\n")
        fileId.write(");\n")

        fileId.write("frequency          nonuniform List<scalar>\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            fileId.write("\t" + str(omg[i]) + "\n")
        fileId.write(");\n")

        fileId.write("waveNumber           nonuniform List<vector>\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            fileId.write("\t(" + str(k[i][0]) + " 0" + " 0)" + "\n")
        fileId.write(");\n")

    elif(solver == "olaFlow"):

        fileId.write("waveReals\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            j = i + 1
            fileId.write(str(real[j]) + "\n")
        fileId.write(");\n")

        fileId.write("waveImags\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            j = i + 1
            fileId.write(str(imag[j]) + "\n")
        fileId.write(");\n")

        fileId.write("wavePeriods\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            j = i + 1
            fileId.write(str(1.0/frq[j]) + "\n")
        fileId.write(");\n")

        fileId.write("waveDirs\n")
        fileId.write(str(n) + "\n")
        fileId.write("(\n")
        for i in range(n):
            fileId.write(str(0) + "\n")
        fileId.write(");\n")

    fileId.close()

def waveNum(frq,depth):
    n = len(frq)//2
    k = np.ones(shape=(n,1))
    omega = 2.0*np.pi*frq
    for i in range(n):
        if (i > 0):
            tmpK = 0
            while (abs(k[i] - tmpK) > 0.00001):
                tmpK = k[i]
                period = 2.0 * np.pi / omega[i]
                waveLen = 9.81 * period ** 2 / (2.0 * np.pi) * np.tanh(tmpK * depth)
                k[i] = 2.0 * np.pi / waveLen
    return k

def readData(dataPath):

    allData = np.loadtxt(dataPath)

    return (allData[:,0],allData[:,1])

if __name__ == '__main__':
    pwdRoot = os.getcwd()
    depth = 0.75
    nUsed = 1500
    startTime = 20
    endTime = 900
    solver = 'olaFlow'
    checkSpectrum = True #Adujst the value of nUsed by examing the spectrum and reconstructed time series data.
    checkEta = True
    #writeToFile = False
    writeToFile = True

    dataPath = os.path.join(pwdRoot,'C4_Gauge01.txt')
    writePath = os.path.join(pwdRoot,'inletWaveCoeff')

    [time,data] = readData(dataPath)
    #data = data - np.mean(data)
    data /= 100.0       # Convert to the International System of Units

    dataOld = data.copy()
    timeOld = time.copy()
    time = timeOld[np.where(timeOld >= startTime)]
    data = dataOld[np.where(timeOld >= startTime)]

    timeOld = time.copy()
    dataOld = data.copy()
    time = timeOld[np.where(timeOld < endTime)]
    data = dataOld[np.where(timeOld < endTime)]

    N = len(time)
    if (N % 2 == 1):
        N -= 1

    time = time[:N]
    data = data[:N]

    dt = time[1] - time[0]
    frq = np.arange(0,N)/(N*dt)

    K = waveNum(frq,depth)

    fftValue = np.fft.fft(data)/(N/2)
    realPart = fftValue.real
    imagPart = fftValue.imag

    if (checkSpectrum):
        amp = abs(fftValue)
        amp[0] /= 2
        plt.figure()
        plt.plot(frq[:nUsed],amp[:nUsed])

    if (checkEta):
        eta = np.zeros(shape=data.shape)
        omega = 2.0*np.pi*frq
        for i in range(N):
            t = i * dt
            for j in range(nUsed):
                if j > 0:
                    eta[i] += realPart[j] * np.cos(omega[j]*t) - imagPart[j] * np.sin(omega[j]*t)
        plt.figure()
        plt.plot(time,eta,'b',label="reconstructed")
        plt.plot(time,data,'r',label="measured")
        plt.legend()

    if (checkEta or checkSpectrum):
        plt.show()

    if (writeToFile):
        writeData(writePath,realPart,imagPart,frq,K,nUsed,solver)