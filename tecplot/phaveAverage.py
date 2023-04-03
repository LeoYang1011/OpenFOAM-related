#Convert Openfoam data to ensight.
#foamToEnsight -fields "(vorticity)" -time 162,167,174,180,186,192,198
#Requires pre-installed python (>=3.7) and pytecplot tools
#The program is developed by Dr. Bo Yang at OUC 2022-11-13

import os
import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

def saveAveragedData():

    print('\nSaving the averaged data to ' + outputName)

    varSave = [dataSet.variable(V)
               for V in ('x','y','z',dataName,resultName)]

    zoneSave = dataSet.zone(aveZones[0])

    tp.data.save_tecplot_plt(outputName,dataset=dataSet,variables=varSave,zones=[zoneSave])

def fieldAverage():
    aveTimes = list()
    for i in aveZones:
        aveTimes.append(dataSet.zone(i).solution_time)

    timesPrint = ''
    for time in aveTimes:
        timesPrint = timesPrint + str(time) + '; '

    print('\nPerform field averaging of ' +  dataName + ' at ' + str(len(aveZones)) +' solution times.')
    print('Averaged zone name: ' + zoneName + '.')
    print('Averaged solution time: ' + timesPrint)
    print('The result will be shown at the solution time: ' + str(aveTimes[0]))

    tp.data.operate.execute_equation(equation='{sum_data} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[aveZones[0]])
    tp.data.operate.execute_equation(equation='{' + resultName + '} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[aveZones[0]])

    for i in aveZones:
        tp.data.operate.execute_equation(equation='{sum_data} = {sum_data} + {' + dataName + '}[' + str(i+1) + ']',
            ignore_divide_by_zero=True,
            value_location=ValueLocation.CellCentered,
            zones=[aveZones[0]])

    tp.data.operate.execute_equation(equation= '{' + resultName + '} = {sum_data}/' + str(len(aveZones)),
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[aveZones[0]])

def findAverageZones():
    allZones = dataSet.zone_names
    aveZonesIndex = list()

    for i in range(len(allZones)):
        if allZones[i] == zoneName:
            aveZonesIndex.append(i)

    return aveZonesIndex

if __name__ == '__main__':
    tp.session.connect()

    pwdRoot = os.getcwd()

    dataName = 'vorticityZ'
    resultName = 'avgVorticityZ'
    zoneName = 'Surf: pipeline'

    inputName = os.path.join(pwdRoot, '2_pimpleFoam.case')
    outputName = os.path.join(pwdRoot, resultName + '.plt')

    dataSet = tp.data.load_ensight(inputName, read_data_option=ReadDataOption.ReplaceInActiveFrame)

    aveZones = findAverageZones()
    fieldAverage()
    saveAveragedData()
