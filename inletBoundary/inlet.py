import pandas as pd
import numpy as np
import math

u0 = 0.37
kappa = 0.42
d = 0.36/1000
y0 = 2.5*d/30.0
delta = 0.4
ymin = -0.05
cmu = 0.09

uStar = kappa*u0/math.log(delta/y0)
print(uStar)
inletPoint = pd.read_csv("inlet0.csv")
xPoint = inletPoint["Points:0"]
yPoint = inletPoint["Points:1"]
zPoint = inletPoint["Points:2"]
yrel = yPoint - ymin
yfirst = min(filter(lambda yrel: yrel > 0, yrel))

uy = np.zeros(len(yrel))
for i in range(len(yrel)):
    if yrel[i] == 0.0:
        uy[i] = 0
        yrel[i] = yfirst
    else:
        uy[i] = min(uStar/kappa*math.log(yrel[i]/y0),u0)


ky = np.zeros(len(yrel))
omegay = np.zeros(len(yrel))
epsilony = np.zeros(len(yrel))
ly = np.zeros(len(yrel))
for i in range(len(yrel)):
    ly[i] = min(kappa*yrel[i]/(1.0 + 1.5*yrel[i]/delta),cmu*delta)
    ky[i] = max(uStar ** 2/np.sqrt(cmu)*(1.0 - yrel[i]/delta) ** 2,0.0005*u0 ** 2)
    omegay[i] = np.sqrt(ky[i])/np.power(cmu,0.25)/ly[i]
    epsilony[i] = np.power(ky[i], 1.5) * np.power(cmu, 0.75) / ly[i]

#print(ky)
pointHeader = '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\    /   O peration     | Version:     4.1                                |
|   \\\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorField;
    object      points;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'''

pointId = open('./boundaryData/inlet/points', 'w')
pointId.write(pointHeader)
pointId.write("\n")
pointId.write("(\n")
pointId.write("\n")
for i in range(len(yPoint)):
    point = " " + "(" + str(xPoint[i]) + " " + str(yPoint[i]) + " " + str(zPoint[i]) + ")" + "\n"
    pointId.write(point)

pointId.write("\n")
pointId.write(")\n\n\n")
pointId.write("// ************************************************************************* //\n")
pointId.close()

UHeader = '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\    /   O peration     | Version:     4.1                                |
|   \\\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       vectorAverageField;
    object      values;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'''

UId = open('./boundaryData/inlet/0/U', 'w')
UId.write(UHeader)
UId.write("\n")
UId.write("// Average\n")
UId.write("(0 0 0)\n")
UId.write("\n")
UId.write("// Data on points\n")
UId.write(str(len(uy)))
UId.write("\n")
UId.write("(\n\n")

for i in range(len(uy)):
    U = "(" + str(uy[i]) + " " + str(0) + " " + str(0) + ")" + "\n"
    UId.write(U)

UId.write("\n")
UId.write(")\n\n\n")
UId.write("// ************************************************************************* //\n")
UId.close()

kOmegaHeader = '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\    /   O peration     | Version:     4.1                                |
|   \\\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       scalarAverageField;
    object      values;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n'''

kId = open('./boundaryData/inlet/0/k', 'w')
kId.write(kOmegaHeader)
kId.write("\n")
kId.write("// Average\n")
kId.write("0.0\n")
kId.write("\n")
kId.write("// Data on points\n")
kId.write(str(len(ky)))
kId.write("\n")
kId.write("(\n\n")

for i in range(len(uy)):
    kId.write(str(ky[i]) + "\n")

kId.write("\n")
kId.write(")\n\n\n")
kId.write("// ************************************************************************* //\n")
kId.close()

omegaId = open('./boundaryData/inlet/0/omega', 'w')
omegaId.write(kOmegaHeader)
omegaId.write("\n")
omegaId.write("// Average\n")
omegaId.write("0.0\n")
omegaId.write("\n")
omegaId.write("// Data on points\n")
omegaId.write(str(len(ky)))
omegaId.write("\n")
omegaId.write("(\n\n")

for i in range(len(uy)):
    omegaId.write(str(omegay[i]) + "\n")

omegaId.write("\n")
omegaId.write(")\n\n\n")
omegaId.write("// ************************************************************************* //\n")
omegaId.close()

epsilonId = open('./boundaryData/inlet/0/epsilon', 'w')
epsilonId.write(kOmegaHeader)
epsilonId.write("\n")
epsilonId.write("// Average\n")
epsilonId.write("0.0\n")
epsilonId.write("\n")
epsilonId.write("// Data on points\n")
epsilonId.write(str(len(ky)))
epsilonId.write("\n")
epsilonId.write("(\n\n")

for i in range(len(uy)):
    epsilonId.write(str(epsilony[i]) + "\n")

epsilonId.write("\n")
epsilonId.write(")\n\n\n")
epsilonId.write("// ************************************************************************* //\n")
epsilonId.close()