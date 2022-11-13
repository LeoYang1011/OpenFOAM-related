#Convert Openfoam data to insight. One zone contains data at one time
#foamToEnsight -no-boundary -no-cellZones -fields "(vorticity)" -time 162,167,174,180,186,192,198
#Requires pre-installed python (>=3.7) and pytecplot tools
#The program is developed by Dr. Bo Yang at OUC 2022-11-13

import os
import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *

def phaseAverage():

    tp.data.operate.execute_equation(equation='{sum_data} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered)
    tp.data.operate.execute_equation(equation='{' + resultName + '} = 0',
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered)

    for i in range(dataSet.num_zones):
        tp.data.operate.execute_equation(equation='{sum_data} = {sum_data} + {' + dataName + '}[' + str(i+1) + ']',
            ignore_divide_by_zero=True,
            value_location=ValueLocation.CellCentered)

    tp.data.operate.execute_equation(equation= '{' + resultName + '} = {sum_data}/' + str(dataSet.num_zones),
        ignore_divide_by_zero=True,
        value_location=ValueLocation.CellCentered)

if __name__ == '__main__':
    tp.session.connect()

    pwdRoot = os.getcwd()

    inputName = os.path.join(pwdRoot, '2_pimpleFoam.case')
    dataName = 'vorticityZ'
    resultName = 'avgVorticityZ'

    dataSet = tp.data.load_ensight(inputName,read_data_option=ReadDataOption.ReplaceInActiveFrame)

    phaseAverage()

    varSave = [dataSet.variable(V)
               for V in ('x','y','z',dataName,resultName)]

    print(dataSet.variable('y').values(-1).shape)
    zoneSave = dataSet.zone(dataSet.num_zones-1)

    tp.data.save_tecplot_plt(resultName + '.plt',dataset=dataSet,variables=varSave,zones=[zoneSave])
