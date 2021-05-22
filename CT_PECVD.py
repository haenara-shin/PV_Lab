# coding=utf-8

from veusz.plugins import *
from datetime import *
import os
from os import path

class PECVD_Log(ImportPlugin):
    """Plugin for importing log files from Centrotherm PECVD"""

    name = 'PECVD Log Plugin'
    description = 'Reads in Centrotherm PECVD log file'
    promote = 'CT PECVD Log'

    file_extensions = set(['.txt'])

    def __init__(self):
        ImportPlugin.__init__(self)

        self.fields = [
            ImportFieldCombo('useheader',
                             items=('0_No_prefix',
                                    '1_User-specified'),
                             editable=False,
                             default='0_No_prefix'),
            ImportFieldText('header', descr='User-specified prefix for dataset names', default=''),
            ]

    def doImport(self, params):

        Status = {'Step': '', 'DateTime0': '', 'T_Paddle_Set': {'LZ': 0., 'LCZ': 0., 'CZ': 0., 'CGZ': 0., 'GZ': 0.}}

        Step = []
        DateTime = []
        TimeDelta = []

        T_Spike = {'LZ': [], 'LCZ': [], 'CZ': [], 'CGZ': [], 'GZ': []}
        T_Paddle = {'LZ': [], 'LCZ': [], 'CZ': [], 'CGZ': [], 'GZ': []}
        T_Paddle_Set = {'LZ': [], 'LCZ': [], 'CZ': [], 'CGZ': [], 'GZ': []}
        T_Paddle_Err = {'LZ': [], 'LCZ': [], 'CZ': [], 'CGZ': [], 'GZ': []}
        ActPwr = {'LZ': [], 'LCZ': [], 'CZ': [], 'CGZ': [], 'GZ': []}
        Gas = {'N2': [], 'SiH4': [], 'NH3': [], 'ScrubbOK': []}
        Pressure = {'TubePres': [], 'AtmPres': [], 'ButFlyPos': [], 'ToxGasFree': []}
        HF = {'Power': [], 'Voltage': [], 'Current': []}

        DataFile = params.openFileWithEncoding()
        readlines = []
        for line in DataFile:
            readlines.append(line.strip())
        DataFile.close()

        i = -1
        while i < len(readlines) - 1:
            i += 1
            line = readlines[i]
            if line[:1].isdigit():
                timeStr = line[:19]
                dataStr = line[20:]

                # find process start line
                Header = 'L Memory            Text1      process started'
                if dataStr.startswith(Header):
                    print timeStr
                    Status['DateTime0'] = datetime.strptime(timeStr, "%m/%d/%Y %H:%M:%S")

                # read current step
                Header = 'L Memory            Text1      '
                if dataStr.startswith(Header):
                    Status['Step'] = dataStr[len(Header):]
                
                # read temperature setpoints
                Header = 'L T_Paddle LoadZone Setpoint   '
                if dataStr.startswith(Header):
                    Status['T_Paddle_Set']['LZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                    
                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    Header = 'L T_Paddle Cent-LZ  Setpoint   '
                    if dataStr.startswith(Header):
                        Status['T_Paddle_Set']['LCZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                        
                        i += 1
                        line = readlines[i]
                        if line[:1].isdigit():
                            dataStr = line[20:]
                    
                    Header = 'L T_Paddle Center   Setpoint   '
                    if dataStr.startswith(Header):
                        Status['T_Paddle_Set']['CZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                        
                        i += 1
                        line = readlines[i]
                        if line[:1].isdigit():
                            dataStr = line[20:]
                    
                    Header = 'L T_Paddle Cent-GZ  Setpoint   '
                    if dataStr.startswith(Header):
                        Status['T_Paddle_Set']['CGZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                        
                        i += 1
                        line = readlines[i]
                        if line[:1].isdigit():
                            dataStr = line[20:]

                    Header = 'L T_Paddle GasZone  Setpoint   '
                    if dataStr.startswith(Header):
                        Status['T_Paddle_Set']['GZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)

                # read temperature setpoints
                Header = 'L T_Paddle All      Setpoint   '
                if dataStr.startswith(Header):
                    Status['T_Paddle_Set']['LZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                    Status['T_Paddle_Set']['LCZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                    Status['T_Paddle_Set']['CZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                    Status['T_Paddle_Set']['CGZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)
                    Status['T_Paddle_Set']['GZ'] = round(float(dataStr[len(Header):].split(' ')[0]), 2)

                # data reading starts with T_Spike
                Header = 'O T_Spike Load..Gas  '
                if dataStr.startswith(Header):
                    DateTime.append(str(datetime.strptime(timeStr, "%m/%d/%Y %H:%M:%S")))
                    TimeDelta.append((datetime.strptime(timeStr, "%m/%d/%Y %H:%M:%S") - Status['DateTime0']).seconds)
                    Step.append(Status['Step'])
                    C1, C2, C3, C4, C5, C6 = dataStr[len(Header):].split('  ')
                    T_Spike['LZ'].append(round(float(C1), 2))
                    T_Spike['LCZ'].append(round(float(C2), 2))
                    T_Spike['CZ'].append(round(float(C3), 2))
                    T_Spike['CGZ'].append(round(float(C4), 2))
                    T_Spike['GZ'].append(round(float(C5), 2))

                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    # T_Paddle
                    Header = 'O T_Paddle Load..Gas  '
                    if dataStr.startswith(Header):
                        C1, C2, C3, C4, C5, C6 = dataStr[len(Header):].split('  ')
                        T_Paddle['LZ'].append(round(float(C1), 2))
                        T_Paddle['LCZ'].append(round(float(C2), 2))
                        T_Paddle['CZ'].append(round(float(C3), 2))
                        T_Paddle['CGZ'].append(round(float(C4), 2))
                        T_Paddle['GZ'].append(round(float(C5), 2))

                        T_Paddle_Set['LZ'].append(Status['T_Paddle_Set']['LZ'])
                        T_Paddle_Set['LCZ'].append(Status['T_Paddle_Set']['LCZ'])
                        T_Paddle_Set['CZ'].append(Status['T_Paddle_Set']['CZ'])
                        T_Paddle_Set['CGZ'].append(Status['T_Paddle_Set']['CGZ'])
                        T_Paddle_Set['GZ'].append(Status['T_Paddle_Set']['GZ'])

                        T_Paddle_Err['LZ'].append(round(float(C1), 2) - Status['T_Paddle_Set']['LZ'])
                        T_Paddle_Err['LCZ'].append(round(float(C2), 2) - Status['T_Paddle_Set']['LCZ'])
                        T_Paddle_Err['CZ'].append(round(float(C3), 2) - Status['T_Paddle_Set']['CZ'])
                        T_Paddle_Err['CGZ'].append(round(float(C4), 2) - Status['T_Paddle_Set']['CGZ'])
                        T_Paddle_Err['GZ'].append(round(float(C5), 2) - Status['T_Paddle_Set']['GZ'])
                    
                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    # ActPwr
                    Header = 'O ActPwr Load..Gas  '
                    if dataStr.startswith(Header):
                        C1, C2, C3, C4, C5 = dataStr[len(Header):].split('  ')
                        ActPwr['LZ'].append(round(float(C1), 2))
                        ActPwr['LCZ'].append(round(float(C2), 2))
                        ActPwr['CZ'].append(round(float(C3), 2))
                        ActPwr['CGZ'].append(round(float(C4), 2))
                        ActPwr['GZ'].append(round(float(C5), 2))

                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    # Gas
                    Header = 'O N2,SiH4,NH3,ScrubbOK  '
                    if dataStr.startswith(Header):
                        C1, C2, C3, C4 = dataStr[len(Header):].split('  ')
                        Gas['N2'].append(round(float(C1), 2))
                        Gas['SiH4'].append(round(float(C2), 2))
                        Gas['NH3'].append(round(float(C3), 2))
                        Gas['ScrubbOK'].append(round(float(C4), 2))

                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    # Pressure
                    Header = 'O TubePres,AtmPres,ButFlyPos,ToxGasFree:  '
                    if dataStr.startswith(Header):
                        C1, C2, C3, C4 = dataStr[len(Header):].split('  ')
                        Pressure['TubePres'].append(round(float(C1), 2))
                        Pressure['AtmPres'].append(round(float(C2), 2))
                        Pressure['ButFlyPos'].append(round(float(C3), 2))
                        Pressure['ToxGasFree'].append(round(float(C4), 2))

                    i += 1
                    line = readlines[i]
                    if line[:1].isdigit():
                        dataStr = line[20:]

                    # HF
                    Header = 'O HF: Power [W], Voltage [V], Current [A]  '
                    if dataStr.startswith(Header):
                        C1, C2, C3 = dataStr[len(Header):].split('  ')
                        HF['Power'].append(round(float(C1), 2))
                        HF['Voltage'].append(round(float(C2), 2))
                        HF['Current'].append(round(float(C3), 2))

        Prefix = ''
        if (params.field_results['useheader'][:1] == '1'):
            if (params.field_results['header'] != ''):
                Prefix = params.field_results['header'] + '-'

        result = [ImportDatasetText(Prefix + '00_DateTime', DateTime),
                  ImportDataset1D(Prefix + '01_TimeDelta', TimeDelta),
                  ImportDataset1D(Prefix + '02_T_Spike_LZ', T_Spike['LZ']),
                  ImportDataset1D(Prefix + '03_T_Spike_LCZ', T_Spike['LCZ']),
                  ImportDataset1D(Prefix + '04_T_Spike_CZ', T_Spike['CZ']),
                  ImportDataset1D(Prefix + '05_T_Spike_CGZ', T_Spike['CGZ']),
                  ImportDataset1D(Prefix + '06_T_Spike_GZ', T_Spike['GZ']),
                  ImportDataset1D(Prefix + '07_T_Paddle_LZ', T_Paddle['LZ']),
                  ImportDataset1D(Prefix + '08_T_Paddle_LCZ', T_Paddle['LCZ']),
                  ImportDataset1D(Prefix + '09_T_Paddle_CZ', T_Paddle['CZ']),
                  ImportDataset1D(Prefix + '10_T_Paddle_CGZ', T_Paddle['CGZ']),
                  ImportDataset1D(Prefix + '11_T_Paddle_GZ', T_Paddle['GZ']),
                  ImportDataset1D(Prefix + '12_T_Paddle_Set_LZ', T_Paddle_Set['LZ']),
                  ImportDataset1D(Prefix + '13_T_Paddle_Set_LCZ', T_Paddle_Set['LCZ']),
                  ImportDataset1D(Prefix + '14_T_Paddle_Set_CZ', T_Paddle_Set['CZ']),
                  ImportDataset1D(Prefix + '15_T_Paddle_Set_CGZ', T_Paddle_Set['CGZ']),
                  ImportDataset1D(Prefix + '16_T_Paddle_Set_GZ', T_Paddle_Set['GZ']),
                  ImportDataset1D(Prefix + '17_T_Paddle_Err_LZ', T_Paddle_Err['LZ']),
                  ImportDataset1D(Prefix + '18_T_Paddle_Err_LCZ', T_Paddle_Err['LCZ']),
                  ImportDataset1D(Prefix + '19_T_Paddle_Err_CZ', T_Paddle_Err['CZ']),
                  ImportDataset1D(Prefix + '20_T_Paddle_Err_CGZ', T_Paddle_Err['CGZ']),
                  ImportDataset1D(Prefix + '21_T_Paddle_Err_GZ', T_Paddle_Err['GZ']),
                  ImportDataset1D(Prefix + '22_ActPwr_LZ', ActPwr['LZ']),
                  ImportDataset1D(Prefix + '23_ActPwr_LCZ', ActPwr['LCZ']),
                  ImportDataset1D(Prefix + '24_ActPwr_CZ', ActPwr['CZ']),
                  ImportDataset1D(Prefix + '25_ActPwr_CGZ', ActPwr['CGZ']),
                  ImportDataset1D(Prefix + '26_ActPwr_GZ', ActPwr['GZ']),
                  ImportDataset1D(Prefix + '27_Gas_N2', Gas['N2']),
                  ImportDataset1D(Prefix + '28_Gas_SiH4', Gas['SiH4']),
                  ImportDataset1D(Prefix + '29_Gas_NH3', Gas['NH3']),
                  ImportDataset1D(Prefix + '30_Gas_ScrubbOK', Gas['ScrubbOK']),
                  ImportDataset1D(Prefix + '31_Pressure_TubePres', Pressure['TubePres']),
                  ImportDataset1D(Prefix + '32_Pressure_AtmPres', Pressure['AtmPres']),
                  ImportDataset1D(Prefix + '33_Pressure_ButFlyPos', Pressure['ButFlyPos']),
                  ImportDataset1D(Prefix + '34_Pressure_ToxGasFree', Pressure['ToxGasFree']),
                  ImportDataset1D(Prefix + '35_HF_Power', HF['Power']),
                  ImportDataset1D(Prefix + '36_HF_Voltage', HF['Voltage']),
                  ImportDataset1D(Prefix + '37_HF_Current', HF['Current']),
                  ImportDatasetText(Prefix + '38_Step', Step)]
        
        return result
            
importpluginregistry.append(PECVD_Log())
        
