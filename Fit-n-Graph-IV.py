#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pyodbc, glob, os, codecs, sys, numpy, matplotlib, math, traceback
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

Set =   {   'GraphOutputTo':    'Graphs',
            'DataOutputTo':     'SweepData'}
# Fields
ReportFields = ['Title', 'ID', 'Comment',
            'Eta', 'Uoc', 'Isc', 'Jsc', 'FF', 'Rser', 'Rshunt', 'InsolMpp',
            'pFF', 'Area', 'J01', 'J02', 'Jph', 'RserType', 'RshuntType',
            'TestTime', 'TestDate', 'Classification',
            'UocUncorr','IscUncorr', 'Pmpp', 'Umpp', 'Impp', 'Jmpp',
            'URev1', 'IRev1', 'URev2', 'IRev2', 'IRevmax', 'URevmax',
            'RshuntDf', 'RserDf', 'RserLfDf',
            'Insol',  'InsolVoc', 'InsolFlashControl',
            'Tmonicell', 'MonicellParamTkI', 'MonicellMvToInsol',
            'Tcell', 'CellParamTkI', 'CellParamTkU', 'CellParamArea']

DbFields = ['Title', 'Comment', 'TestTime', 'TestDate', 'ID', 'Classification',
            'Uoc', 'UocUncorr', 'Isc', 'Jsc', 'IscUncorr', 'Pmpp', 'Umpp', 'Impp', 'Jmpp', 'FF', 'Eta',
            'URev1', 'IRev1', 'URev2', 'IRev2', 'IRevmax', 'URevmax',
            'Rser', 'Rshunt', 'RshuntDf', 'RserDf', 'RserLfDf',
            'Insol', 'InsolMpp', 'InsolVoc', 'InsolFlashControl', 'Tmonicell', 'MonicellParamTkI', 'MonicellMvToInsol',
            'Tcell', 'CellParamTkI', 'CellParamTkU', 'CellParamArea']

CellFields = DbFields + ['RserType', 'RshuntType']

FloatFields = ['Uoc', 'UocUncorr', 'Isc', 'Jsc', 'IscUncorr', 'Pmpp', 'Umpp', 'Impp', 'Jmpp', 'FF', 'Eta',
               'URev1', 'IRev1', 'URev2', 'IRev2', 'IRevmax', 'URevmax',
               'Rser', 'Rshunt', 'RshuntDf', 'RserDf', 'RserLfDf',
               'Insol', 'InsolMpp', 'InsolVoc', 'InsolFlashControl', 'Tmonicell', 'MonicellParamTkI', 'MonicellMvToInsol',
               'Tcell', 'CellParamTkI', 'CellParamTkU', 'CellParamArea',
               'pFF', 'Area', 'J01', 'J02', 'Jph']

# Matplotlib settings
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['lines.linewidth'] = 0.3 * 72 / 25.4
matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['legend.fontsize'] = 10
matplotlib.rcParams['font.serif'] = ('바른돋움Pro 2', 'Arial', 'serif')
matplotlib.rcParams['text.latex.unicode'] = True
matplotlib.rcParams['axes.linewidth'] = 0.45 * 72 / 25.4
matplotlib.rcParams['xtick.major.size'] = 1.0 * 72 / 25.4
matplotlib.rcParams['xtick.major.width'] = 0.3 * 72 / 25.4
matplotlib.rcParams['ytick.major.size'] = 1.0 * 72 / 25.4
matplotlib.rcParams['ytick.major.width'] = 0.3 * 72 / 25.4
matplotlib.rcParams['mathtext.default'] = 'regular'

# Utility functions
def PrintAndLog(text, logfile):
    print text
    logfile.write(text + "\n")

def ISODateTimeStr(TestDate, TestTime):
    items = TestDate.split('.')
    return "%s-%s-%s %s+09" % (items[2], items[1], items[0], TestTime)
    
# Graphing functions
def GraphSaveToFile(fignumber, filename):
    plt.figure(fignumber)
    plt.savefig(filename, dpi=300)

def PltFigText(x, y, label, valtext, unit):
    plt.figtext(x, y, label)
    plt.figtext(x+0.1, y, '=')
    plt.figtext(x+0.15, y, valtext, family='monospace')
    plt.figtext(x+0.27, y, unit)

def RserGraphDraw(LightI, LightV, ShiftedDarkI, DarkV, CorrDarkV, CellInfo, fignumber):
    RserNote = {'RserLf': 'Light forward slope$^{-1}$ at $V_{OC}$', 'RserLfDfIEC': 'Light/Dark IEC 891', 'RserLfDfIECcorr': 'Light/Dark IEC 891 + Dicker'}
    RshuntNote = {'RshuntLf': 'Light I slope$^{-1}$ at zero voltage', 'RshuntDf': 'Dark I slope$^{-1}$ at zero voltage'}
    
    plt.figure(fignumber)
    fig = plt.gcf()
    plt.clf()

    fig.set_size_inches( 14./2.54, 5.5 / 2.54) # 14cm x 5.5cm
    plt.axes([0.09, 0.2, 0.39, 0.75])
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A]')
    xmin = 0.05 * math.floor(CellInfo['Umpp'] * 0.8 / 0.05)
    xmax = 0.05 * math.ceil(CellInfo['Uoc']/0.05)
    plt.xlim(xmin, xmax)
    plt.ylim(0., 10.5)
    plt.xticks(numpy.arange(xmin, xmax + 0.01, 0.05))
    plt.yticks(numpy.arange(0., 11., 2.))

    plt.scatter(LightV, LightI, c='#004098', marker='o', lw=0.1, s = 2.5, label='Light') # light I-V
    plt.scatter(DarkV, ShiftedDarkI, c='#505050', marker='o', lw=0.1, s = 2.5, label='Dark, Shifted') # Isc-shifted dark I-V

    plt.scatter(CorrDarkV, ShiftedDarkI, c='#008736', marker='o', lw=0.1, s = 2.5, label="Dark, Shifted\n+ Dicker Correction") # Isc-shifted dark I versus Rs-corrected dark V

    plt.legend(loc='lower left', fontsize='x-small', frameon=False)

    plt.hlines([CellInfo['Impp']], 0., 1., colors="#404040", linestyles='dotted', linewidth=0.4)
    plt.vlines([CellInfo['Umpp']], 0., 12., colors="#404040", linestyles='dotted', linewidth=0.4)
    plt.vlines([CellInfo['Umpp'] + CellInfo['Impp'] * CellInfo['RserLfDfIEC']], 0., 12., colors="#404040", linestyles='dotted', linewidth=0.4)
    plt.vlines([CellInfo['Umpp'] + CellInfo['Impp'] * CellInfo['RserLfDfIECcorr']], 0., 12., colors="#404040", linestyles='dotted', linewidth=0.4)

    plt.figtext(0.285, 0.95, ISODateTimeStr(CellInfo['TestDate'], CellInfo['TestTime']), horizontalalignment='center', verticalalignment='top', size='x-small', bbox=dict(edgecolor='black', facecolor='white'), family='monospace')
    plt.figtext(0.52, 0.98, "%s_%s (%.2f $cm^2$)" % (CellInfo['Title'], CellInfo['ID'].zfill(6), CellInfo['Area']), verticalalignment='top')
    plt.figtext(0.52, 0.87, "%s" % (CellInfo['Comment']), verticalalignment='top', size='small')

    propY = 0.70
    propDY = -0.09

    PltFigText(0.55, propY, '$R_{s}$', ("%.3f" % (CellInfo['Rser'] * 1.e3)).rjust(7), '$m\Omega$')
    propY += propDY
    PltFigText(0.55, propY, '$r_{s}$', ("%.3f" % (CellInfo['Rser'] * CellInfo['Area'])).rjust(7), '$\Omega\cdot{}cm^{2}$')
    propY += propDY
    plt.figtext(0.55, propY, "$R_{s}$ type: %s" % RserNote[CellInfo['RserType']], size='small')
    propY += 2. * propDY
    PltFigText(0.55, propY, '$R_{p}$', ("%.2f" % (CellInfo['Rshunt'])).rjust(6), '$\Omega$')
    propY += propDY
    PltFigText(0.55, propY, '$r_{p}$', ("%.3f" % (CellInfo['Rshunt'] * CellInfo['Area'] * 1.e-3)).rjust(7), '$k\Omega\cdot{}cm^{2}$')
    propY += propDY
    plt.figtext(0.55, propY, "$R_{p}$ type: %s" % RshuntNote[CellInfo['RshuntType']], size='small')

def LightIVGraphDraw(I, V, P, CellInfo, fignumber):
    plt.figure(fignumber)
    fig = plt.gcf()
    plt.clf()

    fig.set_size_inches( 14./2.54, 5.5 / 2.54 ) # 14cm x 5.5cm
    plt.axes([0.09, 0.2, 0.39, 0.75])
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current [A], Power [W]')
    plt.xlim(0., 0.72)
    plt.ylim(0., 10.5)
    plt.xticks(numpy.arange(0.0, 0.75, 0.2))
    plt.yticks(numpy.arange(0., 11., 2.))

    plt.scatter(V, P, c='#008736', marker='o', lw=0.1, s = 2.5) #I-P curve
    plt.scatter(V, I, c='#004098', marker='o', lw=0.1, s = 2.5) #I-V curve

    plt.vlines([CellInfo['Umpp'], CellInfo['Uoc']], 0., 12., colors="#404040", linestyles='dotted', linewidth=0.4)
    plt.hlines([CellInfo['Pmpp'], CellInfo['Impp'], CellInfo['Isc']], 0., 1., colors="#404040", linestyles='dotted', linewidth=0.4)

    plt.figtext(0.285, 0.95, ISODateTimeStr(CellInfo['TestDate'], CellInfo['TestTime']), horizontalalignment='center', verticalalignment='top', size='x-small', bbox=dict(edgecolor='black', facecolor='white'), family='monospace')

    plt.figtext(0.52, 0.98, "%s_%s (%.2f $cm^2$)" % (CellInfo['Title'], CellInfo['ID'].zfill(6), CellInfo['Area']), verticalalignment='top')
    plt.figtext(0.52, 0.87, "%s" % (CellInfo['Comment']), verticalalignment='top', size='small')

    propY = 0.70
    propDY = -0.09
    PltFigText(0.55, propY, '$V_{OC}$', ("%.1f" % (CellInfo['Uoc'] * 1000.)).rjust(5), '$mV$')
    propY += propDY
    PltFigText(0.55, propY, '$I_{SC}$', ("%.3f" % CellInfo['Isc']).rjust(7), '$A$')
    propY += propDY
    PltFigText(0.55, propY, '$J_{SC}$', ("%.2f" % CellInfo['Jsc']).rjust(6), '$mA/cm^{2}$')
    propY += propDY
    PltFigText(0.55, propY, '$FF$', ("%.2f" % CellInfo['FF']).rjust(6), '$\%$')
    propY += propDY
    PltFigText(0.55, propY, '$R_{s}$', ("%.3f" % (CellInfo['Rser']*1.e3)).rjust(7), '$m\Omega$')
    propY += propDY
    PltFigText(0.55, propY, '$R_{p}$', ("%.2f" % CellInfo['Rshunt']).rjust(6), '$\Omega$')
    propY += propDY
    PltFigText(0.55, propY, '$pFF$', ("%.2f" % CellInfo['pFF']).rjust(6), '$\%$')
    propY += propDY * 1.2
    PltFigText(0.55, propY, '$\eta$', ("%.2f" % CellInfo['Eta']).rjust(6), '$\%$')

def TwoDiodeModelJVGraphDraw(JmA, V, JmASolved, VSolved, JmAUnsolved, VUnsolved, CellInfo, fignumber):
    plt.figure(fignumber)
    fig = plt.gcf()
    plt.clf()

    fig.set_size_inches( 14./2.54, 5.5 / 2.54 ) # 14cm x 5.5cm
    plt.axes([0.09, 0.2, 0.39, 0.75])
    plt.xlabel('Voltage [V]')
    plt.ylabel('J [mA/cm$^2$]')
    plt.xlim(0., 0.72)
    plt.ylim(0., 45.)
    plt.xticks(numpy.arange(0.0, 0.75, 0.2))
    plt.yticks(numpy.arange(0., 45.1, 10.))

    plt.scatter(V, JmA, c='#004098', marker='o', lw=0.1, s = 2.5, label='Data') #J-V curve
    plt.vlines([CellInfo['Umpp'], CellInfo['Uoc']], 0., 12./CellInfo['Area']*1.e3, colors="#404040", linestyles='dotted', linewidth=0.4)
    plt.hlines([CellInfo['Impp']/CellInfo['Area']*1.e3, CellInfo['Isc']/CellInfo['Area']*1.e3], 0., 1., colors="#404040", linestyles='dotted', linewidth=0.4)

    plt.plot(VSolved, JmASolved, c='cyan', lw=0.4, label='Two-Diode Model') #Solved
    if len(JmAUnsolved) > 0:
        plt.plot(VUnsolved, JmAUnsolved, c='red', marker='o', lw=0, s=0.25, label='TDM, Out of Bounds') #Unsolved
    plt.legend(loc='lower left', fontsize='x-small', frameon=False)

    plt.figtext(0.285, 0.95, ISODateTimeStr(CellInfo['TestDate'], CellInfo['TestTime']), horizontalalignment='center', verticalalignment='top', size='x-small', bbox=dict(edgecolor='black', facecolor='white'), family='monospace')
    plt.figtext(0.52, 0.98, "%s_%s (%.2f $cm^2$)" % (CellInfo['Title'], CellInfo['ID'].zfill(6), CellInfo['Area']), verticalalignment='top')
    plt.figtext(0.52, 0.87, "%s" % (CellInfo['Comment']), verticalalignment='top', size='small')

    propY = 0.70
    propDY = -0.09
    PltFigText(0.55, propY, '$J_{ph}$', ("%.2f" % CellInfo['Jph']).rjust(6), '$mA/cm^2$')
    propY += propDY
    PltFigText(0.55, propY, '$J_{01}$', ("%.1f" % CellInfo['J01']).rjust(5), '$fA/cm^2$')
    propY += propDY
    PltFigText(0.55, propY, '$J_{02}$', ("%.2f" % CellInfo['J02']).rjust(6), '$nA/cm^2$')
    propY += propDY
    PltFigText(0.55, propY, '$r_{s}$', ("%.2f" % (CellInfo['Rser']*CellInfo['Area'])).rjust(6), '$\Omega\cdot{}cm^{2}$')
    propY += propDY
    PltFigText(0.55, propY, '$r_{p}$', ("%.3f" % (CellInfo['Rshunt']*CellInfo['Area']*1.e-3)).rjust(7), '$k\Omega\cdot{}cm^{2}$')

# Data-loading functions
def CellInfoFromRow(Row):
    result = {}
    for i in range(len(DbFields)):
        if (DbFields[i] in FloatFields) and (Row[i] != None):
            result[DbFields[i]] = float(Row[i])
        else:
            result[DbFields[i]] = Row[i]

    result['ID'] = str(int(result['ID']))
    result['Area'] = result['CellParamArea']*.01

    result['RserType'] = 'RserLf'
    result['RshuntType'] = 'RshuntDf'
    if result['Rshunt'] != result['RshuntDf']:
        result['RshuntType'] = 'RshuntDr'

    return result

def CellInfosFromDb(DBfile):
    result = []

    OldWD = os.getcwd()
    DBfilebasename = os.path.basename(DBfile)
    os.chdir(os.path.dirname(DBfile))

    conn = pyodbc.connect('DRIVER={Microsoft Access Driver (*.mdb)};DBQ=' + DBfilebasename)
    cursor = conn.cursor()
    query = cursor.execute('SELECT ' + ', '.join(DbFields) + ' FROM `Results` ORDER BY `ID`;')
    for row in query:
        result.append(CellInfoFromRow(row))
    cursor.close()
    conn.close()

    os.chdir(OldWD)
    
    return result

def SweepFileName(CellInfo):
    return CellInfo['Title'] + '_' + CellInfo['ID'].zfill(6) + '.prn'

def SweepDataFromFile(Filename):
    result = {}
    numericCol = [u'[mm²] Cell area', u'[°C] T Cell', u'[%/K] Modul-TkI', u'[%/K] Modul TkU', u'[°C] T Monitor Cell', u'[%/K] Monitor Cell TkI', u'[W/m²] corrected to', u'[W/m²/mV] Calibration value', u'Spectral mismatch faktor', u'CellParamTkI_Umpp', u'CellParamTkU_Umpp', u'IEC60891_RsType', u'IEC60891_TkI', u'IEC60891_TkU', u'IEC60891_Kfact', u'[A]Isc', u'[V]Uoc', u'[A]Impp', u'[V]Umpp', u'[W]Pmpp', u'[%]Eta', u'[%]FF', u'[Ohm]Rser', u'[Ohm]Rshunt', u'[W/m²]Irrad.', u'Ufrom[V]', u'Uto[V]', u'Lines', u'Decimals', u'Ifrom[A]', u'Ito[A]', u'[V]Uraw', u'[A]Iraw', u'[mV]Eraw', u'[V]Uflt', u'[A]Iflt', u'[mV]Eflt', u'[V]Ucor', u'[A]Icor', u'[W/m2]Ecor', u'[A]Irev1', u'[V]Urev1', u'[A]Irev2', u'[V]Urev2', u'[A]Irevmax', u'[V]Urefmax']

    with codecs.open(Filename, 'r', encoding='iso-8859-1') as f:
        qIndex = 0
        readerMode = False
        column = []
        colIsFloat = []

        for line in f:
            line = line.strip().rstrip(';')
            if line == "":
                readerMode = False
                column = []
                colIsFloat = []
            else:
                if not (line[0].isdigit() or line[0] == '-'): #line does not start with number
                    if line.startswith('Date;'): #check for the start of a new sweep data
                        qIndex += 1
                        result["Q%d" % qIndex] = {}
                    if qIndex > 0:
                        column = line.split(';') #read in columns
                        for i in range(len(column)):
                            result["Q%d" % qIndex][column[i]] = []
                        readerMode = True
                else: #line starts with number
                    if (qIndex > 0 and readerMode):
                        item = line.split(';')
                        if len(item) == len(column):
                            for i in range(len(item)):
                                itemtext = item[i].strip()
                                if itemtext <> '':
                                    if column[i] in numericCol:
                                        if itemtext <> '-':
                                            result["Q%d" % qIndex][column[i]].append(float(itemtext))
                                    else:
                                        result["Q%d" % qIndex][column[i]].append(itemtext)
    if qIndex <> 4:
        PrintAndLog("RawDataLoadFromFile: Expected 4 data sets but found %d: %s" % (qIndex, argument))

    return result

def InfosAndSweepsFrom(Dir):
    OldWD = os.getcwd()

    infos = {}
    sweepdata = {}
    
    dbfiles = []
    if (os.path.isdir(Dir)):
        dbfiles = glob.glob(os.path.join(Dir, '*.mdb'))
        for dbfile in dbfiles:
            cellinfos = CellInfosFromDb(dbfile)
            title = os.path.basename(os.path.splitext(dbfile)[0])

            datadir = os.path.dirname(dbfile)
            os.chdir(datadir)
            sweeps = {}
            for i in range(len(cellinfos)):
                sweepfn = SweepFileName(cellinfos[i])
                if os.path.isfile(sweepfn):
                    sweeps[cellinfos[i]['ID']] = SweepDataFromFile(sweepfn)

            infos[title] = cellinfos
            sweepdata[title] = sweeps

            cellinfos = {}
            sweeps = {}

    os.chdir(OldWD)

    return infos, sweepdata

# Data analysis functions
def ExpN(V, n):
    return math.expm1(V * 1.60217657e-19 / 1.3806488e-23 / (273.15 + 25.) / n)

def Shift(DarkI, Isc):
    result = []
    for i in range(len(DarkI)):
        result.append(DarkI[i] + Isc)
    return result

def CorrectV(DarkV, DarkI, DarkRs):
    result = []
    for i in range(len(DarkV)):
        result.append(DarkV[i] + DarkI[i] * DarkRs)
    return result

def CellAnalysis(CellInfo, Sweep):
    if 'Q1' in Sweep.keys():
        CellInfo.update(MainValuesFromIV(Sweep['Q1']['[A]Icor'], Sweep['Q1']['[V]Ucor']))
        if 'Q4' in Sweep.keys():
            CellInfo.update(RserValuesFrom(Sweep['Q1']['[A]Icor'], Sweep['Q1']['[V]Ucor'], Sweep['Q4']['[A]Icor'], Sweep['Q4']['[V]Ucor']))
        CellInfo.update(IterateForJ0(CellInfo['Isc'], CellInfo['Uoc'], CellInfo['Impp'], CellInfo['Umpp'], CellInfo['Rser'], CellInfo['Rshunt'], CellInfo['Area'], 500, 0.0001))
    return CellInfo

def JVfromTwoDiodeModel(Jph, J01, J02, Rs, Rp, Vmin, Vmax, Vstep, MaxLoop, Precision):
    J = []
    V = []
    Solved = []
    Loops = []
    
    wNew = 0.75
    wOld = 1. - wNew

    Jval = Jph
    Vval = Vmin

    while Vval <= Vmax:
        counter = 0
        Jvalnew = wOld * Jval + wNew * ( Jph - J01 * ExpN(Vval + Jval * Rs, 1) - J02 * ExpN(Vval + Jval * Rs, 2) - (Vval + Jval * Rs)/Rp )
        Error = abs((Jvalnew - Jval)/Jval)
        while (Error >= Precision) and (counter < MaxLoop):
            counterNew = counter + 1
            if math.floor(counterNew/25.) <> math.floor(counter/25.):
                wNew = 0.5 * wNew
                wOld = 1. - wNew
            counter = counterNew
            Jvalnew = wOld * Jval + wNew * ( Jph - J01 * ExpN(Vval + Jval * Rs, 1) - J02 * ExpN(Vval + Jval * Rs, 2) - (Vval + Jval * Rs)/Rp )
            Error = abs((Jvalnew - Jval)/Jval)
            Jval = Jvalnew

        J.append(Jval)
        V.append(Vval)
        Solved.append(Error < Precision)
        Loops.append(counter)

        Vval += Vstep

    return J, V, Solved, Loops

def IterateForJ0(Isc, Voc, Impp, Vmpp, Rs, Rp, Area, MaxLoop, Precision):
    result = {'DDModelSolved': False}
    wOld = 0.1
    wNew = 1 - wOld

    counter = 0
    error = {'Iph': 1.1 * Precision, 'I01': 1.1 * Precision, 'I02': 1.1 * Precision}

    # 1st calculation
    Iph = Isc * (1. + Rs/Rp)
    I01 = (Iph - Voc/Rp)/ExpN(Voc, 1)
    I02 = (Iph - Impp - I01 * ExpN(Vmpp + Impp * Rs, 1) - (Vmpp + Impp * Rs)/Rp)/ExpN(Vmpp + Impp * Rs, 2)

    # now iterate
    while (counter < MaxLoop) and (not result['DDModelSolved']):
        Iphnew = wOld * Iph + wNew * (Isc + I01 * ExpN(Isc * Rs, 1) + I02 * ExpN(Isc * Rs, 2) + Isc * Rs/Rp)
        I01new = wOld * I01 + wNew * (Iphnew - I02 * ExpN(Voc, 2) - Voc/Rp) / ExpN(Voc, 1)
        I02new = wOld * I02 + wNew * (Iphnew - Impp - I01new * ExpN(Vmpp + Impp * Rs, 1) - (Vmpp + Impp * Rs)/Rp) / ExpN(Vmpp + Impp * Rs, 2)

        error['Iph'] = abs((Iphnew - Iph)/Iph)
        error['I01'] = abs((I01new - I01)/I01)
        error['I02'] = abs((I02new - I02)/I02)
        result['DDModelSolved'] = (error['Iph'] <= Precision) and (error['I01'] <= Precision) and (error['I02'] <= Precision)

        Iph = Iphnew
        I01 = I01new
        I02 = I02new

        counter += 1

    result['Jph'] = Iph/Area * 1.e3
    result['J01'] = I01/Area * 1.e15
    result['J02'] = I02/Area * 1.e9
    return result
    
def RserValuesFrom(LightI, LightV, DarkI, DarkV):
    result = {}
    LightVal = MainValuesFromIV(LightI, LightV, False)
    ShiftedDarkI = Shift(DarkI, LightVal['Isc'])
    DarkVal = MainValuesFromIV(ShiftedDarkI, DarkV, False)

    #Light-forward slope
    if (min(LightI) <= -1.) and (max(LightI) >= 1.):
        sampleV = []
        sampleI = []
        for i in range(len(LightI)):
            if (LightI[i] >= -1.) and (LightI[i] <= 1.):
                sampleI.append(LightI[i])
                sampleV.append(LightV[i])
        if len(sampleI) >= 2:
            coeff = numpy.polyfit(sampleI, sampleV, 1)
            result['RserLf'] = -coeff[0]
            result['Rser'] = result['RserLf']
            result['RserType'] = 'RserLf'

    #Dark-forward slope
    if (min(ShiftedDarkI) <= -1.) and (max(ShiftedDarkI) >= 1.):
        sampleV = []
        sampleI = []
        for i in range(len(ShiftedDarkI)):
            if (ShiftedDarkI[i] >= -1.) and (ShiftedDarkI[i] <= 1.):
                sampleI.append(ShiftedDarkI[i])
                sampleV.append(DarkV[i])
        if len(sampleI) >= 2:
            coeff = numpy.polyfit(sampleI, sampleV, 1)
            result['RserDf'] = -coeff[0]
            result['Rser'] = result['RserDf']
            result['RserType'] = 'RserDf'

    #Light forward Voc and Shifted dark forward Voc comparison
    result['RserLfDf'] = (DarkVal['Uoc']-LightVal['Uoc'])/LightVal['Isc']
    
    #Light forward V at Impp and Dark forward V at Impp comparison
    fVLight = interp1d(LightI, LightV)
    fVDark = interp1d(ShiftedDarkI, DarkV)
    result['RserLfDfIEC'] = (fVDark(LightVal['Impp'])-fVLight(LightVal['Impp']))/(LightVal['Impp'])
    result['Rser'] = result['RserLfDfIEC']
    result['RserType'] = 'RserLfDfIEC'
    
    #Light forward V at Impp and Corrected dark forward V at Impp comparison
    CorrDarkV = CorrectV(DarkV, DarkI, result['RserLfDf'])
    fVCorrDark = interp1d(ShiftedDarkI, CorrDarkV)
    result['RserLfDfIECcorr'] = (fVCorrDark(LightVal['Impp'])-fVLight(LightVal['Impp']))/(LightVal['Impp'])
    result['Rser'] = result['RserLfDfIECcorr']
    result['RserType'] = 'RserLfDfIECcorr'

    pFFvalues = MainValuesFromIV(ShiftedDarkI, CorrDarkV)
    result['pFF'] = pFFvalues['FF']

    return result

def MainValuesFromIV(I, V, FullAnalysis=False):
    Values = {}
    FA = {}
    
    #Isc first guess: I value for smallest abs(V)
    Index = 0
    X = abs(V[Index])
    for i in range(len(V)):
        if abs(V[i]) < X:
            Index = i
            X = abs(V[Index])
    Isc = I[Index]
    FA['Isc'] = {}
    FA['Isc']['Type'] = 0
    FA['Isc']['V'] = [V[Index]]
    FA['Isc']['I'] = [I[Index]]

    #Voc first guess: V value for smallest abs(I)
    Index = 0
    Y = abs(I[Index])
    for i in range(len(I)):
        if abs(I[i]) < Y:
            Index = i
            Y = abs(I[Index])
    Voc = V[Index]
    FA['Voc'] = {}
    FA['Voc']['Type'] = 0
    FA['Voc']['V'] = [V[Index]]
    FA['Voc']['I'] = [I[Index]]

    #Isc guess from cubic spline
    X = []
    Y = []
    XMin = -0.1 * Voc
    XMax = +0.1 * Voc
    while (len(X) < 2) and (XMax < Voc):
        X = []
        Y = []
        for i in range(len(V)):
            if V[i] >= XMin and V[i] <= XMax:
                X.append(V[i])
                Y.append(I[i])
        XMin = 2.0 * XMin
        XMax = 2.0 * XMax
    if len(X) >= 2:
        fIsc = interp1d(X, Y)
        Isc = fIsc(0.)
    FA['Isc']['Type'] = 1
    FA['Isc']['V'] = X
    FA['Isc']['I'] = Y

    #Voc guess from cubic spline
    X = []
    Y = []
    YMin = -0.1 * Isc
    YMax = +0.1 * Isc
    while (len(X) < 2) and (YMax < Isc):
        X = []
        Y = []
        for i in range(len(I)):
            if I[i] >= YMin and I[i] <= YMax:
                X.append(V[i])
                Y.append(I[i])
        YMin = 2.0 * YMin
        YMax = 2.0 * YMax
    if len(X) >= 2:
        fVoc = interp1d(Y, X)
        Voc = fVoc(0.)
    FA['Voc']['Type'] = 1
    FA['Voc']['V'] = X
    FA['Voc']['I'] = Y

    #Impp, Vmpp first guess: largest P
    Index = 0
    Y = I[Index] * V[Index]
    for i in range(len(V)):
        if I[i] * V[i] > Y:
            Index = i
            Y = I[i] * V[i]
    Impp = I[Index]
    Vmpp = V[Index]
    Pmpp = Impp * Vmpp
    FA['MPP'] = {}
    FA['MPP']['Type'] = 0
    FA['MPP']['V'] = [Vmpp]
    FA['MPP']['I'] = [Impp]
    FA['MPP']['P'] = [Pmpp]

    #Impp, Vmpp guess from cubic spline
    YMin = 0.95 * Pmpp
    X = []
    Y = []
    Y2 = []
    for i in range(len(V)):
        if I[i] * V[i] >= YMin:
            X.append(V[i])
            Y.append(I[i] * V[i])
            Y2.append(I[i])
    if len(X) >= 2:
        fPmpp = interp1d(X, Y)
        Vmpp = X[0]
        Pmpp = Y[0]
        for Vval in numpy.arange(min(X), max(X), 0.0001):
            if fPmpp(Vval) > Pmpp:
                Pmpp = fPmpp(Vval)
                Vmpp = Vval
        fImpp = interp1d(X, Y2)
        Impp = fImpp(Vmpp)
    FA['MPP']['Type'] = 1
    FA['MPP']['V'] = X
    FA['MPP']['I'] = Y2
    FA['MPP']['P'] = Y
    
    FF = Pmpp/(Isc * Voc) * 100.

    Values['Isc'] = Isc
    Values['Uoc'] = Voc
    Values['Impp'] = Impp
    Values['Umpp'] = Vmpp
    Values['Pmpp'] = Pmpp
    Values['FF'] = FF

    if FullAnalysis:
        return Values, FA
    else:
        return Values

# Report functions

def OneLine(Cell):
    return "Eta %.2f%% Voc %.1fmV Isc %.3fA Vmpp %.1fmV Impp %.3fA FF %.2f%%" % (Cell['Eta'], Cell['Uoc']*1.e3, Cell['Isc'], Cell['Umpp']*1.e3, Cell['Impp'], Cell['FF'])

def CellReportLine(CellInfo):
    items = []
    for field in ReportFields:
        if field in CellInfo:
            if field in FloatFields:
                if type(CellInfo[field]) != type(None):
                    items.append("%.15g" % CellInfo[field])
                else:
                    items.append('')
            else:
                items.append("%s" % CellInfo[field])
        else:
            items.append('')
    return "\t".join(items)

def ReportSaveToFile(Report, Filename):
    if len(Report) > 1:
        with open(Filename, 'w') as reportFile:
            headerText = "\t".join(ReportFields)
            reportFile.write(headerText + "\n")
            for i in range(0, len(Report)):
                reportFile.write(Report[i] + "\n")

def SweepSaveToFile(Sweep, Filename):
    sweeps = ['Q1', 'Q2', 'Q3', 'Q4']
    fields = ['[V]Uraw', '[A]Iraw', '[mV]Eraw', '[V]Uflt', '[A]Iflt', '[mV]Eflt', '[V]Ucor', '[A]Icor', '[W/m2]Ecor']
    columns = []
    for i in range(len(sweeps)):
        for j in range(len(fields)):
            columns.append("%s-%s" % (sweeps[i], fields[j]))
    maxCount = max(len(Sweep['Q1']['[V]Uraw']), len(Sweep['Q2']['[V]Uraw']), len(Sweep['Q3']['[V]Uraw']), len(Sweep['Q4']['[V]Uraw']))
    
    with open(Filename, 'w') as sweepFile:
        sweepFile.write("\t".join(columns) + "\n")
        for row in range(maxCount):
            items = []
            for i in range(len(sweeps)):
                for j in range(len(fields)):
                    if row < len(Sweep[sweeps[i]][fields[j]]):
                        #items.append("%d_%d_%d" % (i, j, row))
                        items.append("%.9f" % Sweep[sweeps[i]][fields[j]][row])
                    else:
                        items.append("")
            sweepFile.write("\t".join(items) + "\n")

# define main program

def main():
    LogFile = open(os.path.join(os.path.dirname(sys.argv[0]), os.path.splitext(os.path.basename(sys.argv[0]))[0] + '.log'), 'w')
    PrintAndLog("Script: " + os.path.basename(sys.argv[0]), LogFile)
    PrintAndLog("Command-line arguments: " + str(len(sys.argv)-1) + " total", LogFile)

    arguments = []
    for i in range(1, len(sys.argv)):
        PrintAndLog("{0}: {1}".format(i, sys.argv[i]), LogFile)
        arguments.append(sys.argv[i])

    DirIndex = 0
    for argument in arguments:
        DirIndex += 1

        try:
            DataDir = os.path.join(os.path.dirname(argument), Set['DataOutputTo'])
            if not os.path.isdir(DataDir):
                os.makedirs(DataDir)
            GraphDir = os.path.join(os.path.dirname(argument), Set['GraphOutputTo'])
            if not os.path.isdir(GraphDir):
                os.makedirs(GraphDir)

            Infos, Sweeps = InfosAndSweepsFrom(argument)
            if len(arguments) > 1:
                DirCountText = "[d {0}/{1}]".format(DirIndex, len(arguments))
            else:
                DirCountText = ''

            groupcounter = 0
            for Title in Infos.keys():
                groupcounter +=1
                reportData = []
                if len(Infos.keys()) > 1:
                    GroupCountText = "[g {0}/{1}]".format(groupcounter, len(Infos.keys()))
                else:
                    GroupCountText = ''

                cellcounter = 0
                for CellInfo in Infos[Title]:
                    cellcounter += 1
                    PrintAndLog("{0}{1}[{2}/{3}] {4}_{5}".format(DirCountText, GroupCountText, cellcounter, len(Infos[Title]), CellInfo['Title'], CellInfo['ID'].zfill(6)), LogFile)
                    PrintAndLog("- Comment: %s" % CellInfo['Comment'], LogFile)
                    PrintAndLog("- HALM: %s" % OneLine(CellInfo), LogFile)
                    if (Title in Sweeps.keys()) and (CellInfo['ID'] in Sweeps[Title].keys()):
                        CellSweep = Sweeps[Title][CellInfo['ID']]
                        if CellInfo['RshuntType'] == 'RshuntDr':
                            if '[Ohm]Rshunt' in CellSweep['Q1'].keys():
                                if len(CellSweep['Q1']['[Ohm]Rshunt']) > 0:
                                    CellInfo['Rshunt'] = CellSweep['Q1']['[Ohm]Rshunt'][0]
                                    CellInfo['RshuntType'] = 'RshuntLf'
                        CellInfo = CellAnalysis(CellInfo, CellSweep) #Light I-V, Rser, TDM analysis
                        PrintAndLog("- New : %s" % OneLine(CellInfo), LogFile)

                        # Save sweep data in tab-separated text format
                        Dir = os.path.join(os.path.dirname(argument), Set['DataOutputTo'])
                        File = os.path.join(Dir, "%s_%s I-V.txt" % (str(CellInfo['Title']), str(CellInfo['ID']).zfill(6)))
                        SweepSaveToFile(CellSweep, File)
                        PrintAndLog("- Sweep data saved to %s" % File, LogFile)

                        # Light I-V graph
                        P = []
                        for i in range(len(CellSweep['Q1']['[A]Icor'])):
                            P.append(CellSweep['Q1']['[A]Icor'][i] * CellSweep['Q1']['[V]Ucor'][i])

                        Dir = os.path.join(os.path.dirname(argument), Set['GraphOutputTo'])
                        File = os.path.join(Dir, "%s_%s I-V.png" % (str(CellInfo['Title']), str(CellInfo['ID']).zfill(6)))
                        LightIVGraphDraw(CellSweep['Q1']['[A]Icor'], CellSweep['Q1']['[V]Ucor'], P, CellInfo, 0)
                        GraphSaveToFile(0, File)
                        PrintAndLog("- Light I-V graph saved to %s" % File, LogFile)

                        # Resistance analysis graph
                        Dir = os.path.join(os.path.dirname(argument), Set['GraphOutputTo'])
                        File = os.path.join(Dir, "%s_%s R.png" % (str(CellInfo['Title']), str(CellInfo['ID']).zfill(6)))
                        RserGraphDraw(CellSweep['Q1']['[A]Icor'], CellSweep['Q1']['[V]Ucor'], Shift(CellSweep['Q4']['[A]Icor'], CellInfo['Isc']), CellSweep['Q4']['[V]Ucor'], CorrectV(CellSweep['Q4']['[V]Ucor'], CellSweep['Q4']['[A]Icor'], CellInfo['RserLfDf']), CellInfo, 0)
                        GraphSaveToFile(0, File)
                        PrintAndLog("- Rser graph saved to %s" % File, LogFile)

                        #Two-diode model graph
                        JmA = []
                        for i in range(len(CellSweep['Q1']['[A]Icor'])):
                            JmA.append(CellSweep['Q1']['[A]Icor'][i] * 1.e3 / CellInfo['Area'])
                        JmASolved = []
                        VSolved = []
                        JmAUnsolved = []
                        VUnsolved = []
                        Jmodel, Vmodel, Solvedmodel, Loopsmodel = JVfromTwoDiodeModel(CellInfo['Jph']*1.e-3, CellInfo['J01']*1.e-15, CellInfo['J02']*1.e-9, CellInfo['Rser']*CellInfo['Area'], CellInfo['Rshunt']*CellInfo['Area'], 0., CellInfo['Uoc'] + 1.e-3, 1.e-3, 500, 0.0001)
                        for i in range(len(Jmodel)):
                            if Solvedmodel[i]:
                                JmASolved.append(Jmodel[i]*1.e3)
                                VSolved.append(Vmodel[i])
                            else:
                                JmAUnsolved.append(Jmodel[i]*1.e3)
                                VUnsolved.append(Vmodel[i])

                        Dir = os.path.join(os.path.dirname(argument), Set['GraphOutputTo'])
                        File = os.path.join(Dir, "%s_%s TDM.png" % (str(CellInfo['Title']), str(CellInfo['ID']).zfill(6)))
                        TwoDiodeModelJVGraphDraw(JmA, CellSweep['Q1']['[V]Ucor'], JmASolved, VSolved, JmAUnsolved, VUnsolved, CellInfo, 0)
                        GraphSaveToFile(0, File)
                        PrintAndLog("- Two-diode model graph saved to %s" % File, LogFile)

                    reportData.append(CellReportLine(CellInfo))

                Dir = os.path.join(os.path.dirname(argument), Set['DataOutputTo'])
                File = os.path.join(Dir, "%s.txt" % Title)
                ReportSaveToFile(reportData, File)
                PrintAndLog("*Table saved to: %s" % File, LogFile)
        except Exception:
            raise
        
    LogFile.close()


# run main program
try:
    main()
except Exception:
    traceback.print_exc()
    print ""
    raw_input("Press [Enter] to exit")
