import os
import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

def saveAveragedData():

    varSave = [dataSet.variable(V)
               for V in ('x','y','z',dataName,resultName)]

    zoneSave = dataSet.zone(sliceZones[0])

    tp.data.save_tecplot_plt(outputName,dataset=dataSet,variables=varSave,zones=[zoneSave])

def sliceAverage():

    tp.data.operate.execute_equation(equation='{sum_data} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[sliceZones[0]])
    tp.data.operate.execute_equation(equation='{' + resultName + '} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[sliceZones[0]])

    for i in sliceZones:
        tp.data.operate.execute_equation(equation='{sum_data} = {sum_data} + {' + dataName + '}[' + str(i+1) + ']',
            ignore_divide_by_zero=True,
            value_location=ValueLocation.CellCentered,
            zones=[sliceZones[0]])

    tp.data.operate.execute_equation(equation= '{' + resultName + '} = {sum_data}/' + str(len(sliceZones)),
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered,
        zones=[sliceZones[0]])


def extractSlice():

    if spanwiseDir == 'X':
        tp.active_frame().plot().slice(0).orientation = SliceSurface.XPlanes
    elif spanwiseDir == 'Y':
        tp.active_frame().plot().slice(0).orientation = SliceSurface.XPlanes
    elif spanwiseDir == 'Z':
        tp.active_frame().plot().slice(0).orientation = SliceSurface.ZPlanes

    tp.active_frame().plot().slice(0).contour.show = False
    tp.active_frame().plot().slice(0).show_primary_slice = False
    tp.active_frame().plot().slice(0).show_start_and_end_slices = True
    tp.active_frame().plot().slice(0).start_position.z = sliceStart
    tp.active_frame().plot().slice(0).end_position.z = sliceEnd
    tp.active_frame().plot().slice(0).show_intermediate_slices = True
    tp.active_frame().plot().slice(0).num_intermediate_slices = sliceNum
    tp.active_frame().plot().slice(0).show = False
    tp.active_frame().plot().slices(0).extract()

def findSliceZones():
    allZones = dataSet.zone_names
    aveZonesIndex = list()

    for i in range(len(allZones)):
        if zoneName in allZones[i]:
            aveZonesIndex.append(i)

    return aveZonesIndex

if __name__ == '__main__':
    tp.session.connect()

    pwdRoot = os.getcwd()

    sliceStart = 0.02
    sliceEnd = 3.98
    sliceNum = 98
    solutionTime = 266.183
    spanwiseDir = 'Z'
    dataName = 'vorticityZ'
    resultName = 'avgVorticityZ'
    zoneName = 'Slice'

    inputName = os.path.join(pwdRoot, '2_pimpleFoam.case')
    outputName = os.path.join(pwdRoot, resultName + '.plt')

    dataSet = tp.data.load_ensight(inputName, read_data_option=ReadDataOption.ReplaceInActiveFrame)

    extractSlice()
    sliceZones = findSliceZones()
    sliceAverage()
    saveAveragedData()
