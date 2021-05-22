# coding=utf-8

from veusz.plugins import *
from numpy import *
from os import path
import math

class Sentech_ASCII(ImportPlugin):
    """Plugin for importing ASCII files from Sentech Ellipsometer"""

    name = 'Sentech ASCII Plugin'
    description = 'Reads in Sentech ellipsometer map data (ASCII format)'
    promote = 'Sentech'

    file_extensions = set(['.asc'])

    def __init__(self):
        ImportPlugin.__init__(self)

        self.fields = [
            ImportFieldCombo('useheader',
                             items=('0_No_prefix',
                                    '1_Cassette_name',
                                    '2_User-specified'),
                             editable=False,
                             default='0_No_prefix'),
            ImportFieldText('header', descr='User-specified prefix for dataset names', default=''),
            ]

    def HistogramXY(self, binsource, bincount, LSU, ndigits):
        binleft = LSU * math.floor(min(binsource)/LSU)
        binright = LSU * math.ceil(max(binsource)/LSU)
        binwidth = (binright - binleft)/(bincount - 1)
        #binleft -= binwidth
        #binright += binwidth

        bincenters = []
        binvalues = []
        bincenter = binleft
        while (bincenter < binright + 0.5 * binwidth):
            if (ndigits == 0):
                binlabel = "%.0f" % bincenter
            elif (ndigits == 1):
                binlabel = "%.1f" % bincenter
            elif (ndigits == 2):
                binlabel = "%.2f" % bincenter
            elif (ndigits == 3):
                binlabel = "%.3f" % bincenter
            elif (ndigits == 4):
                binlabel = "%.4f" % bincenter
            elif (ndigits == 5):
                binlabel = "%.5f" %  bincenter
                
            bincenters.append(binlabel)
            binvalues.append(0)
            bincenter += binwidth

        for value in binsource:
            binindex = int(round((value - binleft)/binwidth))
            if binindex <0:
                binindex = 0
            elif binindex > len(bincenters) - 1:
                binindex = len(bincenters) - 1
            binvalues[binindex] += 1

        return (bincenters, binvalues)

    def doImport(self, params):
        TextLineTest = ';'
        HeaderLineTest = ';X[mm]'
        IdentifierText = ';cassette_identifier='

        X = []
        Y = []
        D = []
        RI = []
        OL = []
        DError = []
        RIError = []
        OLError = []
        CassetteName = ''
        Prefix = ''
        
        HeaderFound = False

        DataFile = params.openFileWithEncoding()

        # test each line and read in header and data
        for DataLine in DataFile:
            DataLine = DataLine.strip() # remove linebreaks at the end
            LineIsHeader = False
            
            LineIsText = DataLine.startswith(TextLineTest)

            if (LineIsText):
                if (DataLine.startswith(IdentifierText)):
                    CassetteName = DataLine[len(IdentifierText):]

            if ( (not HeaderFound) and LineIsText ): # check if line is header
                LineIsHeader = DataLine.startswith(HeaderLineTest)
                HeaderFound = LineIsHeader

            if (HeaderFound and (not LineIsText)): # read in data values
                Numbers = DataLine.split("\t")
                if (len(Numbers) >= 7):
                    X.append(float(Numbers[0]))
                    Y.append(float(Numbers[1]))
                    D.append(float(Numbers[3]))
                    RI.append(float(Numbers[4]))
                    OL.append(float(Numbers[3]) * float(Numbers[4]))
                if (len(Numbers) >= 9):
                    DError.append(float(Numbers[7]))
                    RIError.append(float(Numbers[8]))
                    OLError.append((DError[len(DError)-1]/D[len(D)-1] + RIError[len(RIError)-1]/RI[len(RI)-1]) * OL[len(OL)-1])

        DataFile.close()

        #form 2D data
        XKeys = {}
        YKeys = {}
        for i in range(len(X)):
            xkey = "%.3f" % X[i]
            if (not (xkey in XKeys)):
                XKeys[xkey] = len(XKeys)
            ykey = "%.3f" % Y[i]
            if (not (ykey in YKeys)):
                YKeys[ykey] = len(YKeys)

        xkeylist = []
        XMin = 0.0
        XMax = 0.0
        for xkey in XKeys:
            xkeylist.append(float(xkey))
        xkeylist.sort()
        if (len(xkeylist) > 1):
            XMin = xkeylist[0] - (xkeylist[1]-xkeylist[0])*0.5
            XMax = xkeylist[len(xkeylist)-1] + (xkeylist[1]-xkeylist[0])*0.5

        ykeylist = []
        YMin = 0.0
        YMax = 0.0
        for ykey in YKeys:
            ykeylist.append(float(ykey))
        ykeylist.sort()
        if (len(ykeylist) > 1):
            YMin = ykeylist[0] - (ykeylist[1]-ykeylist[0])*0.5
            YMax = ykeylist[len(ykeylist)-1] + (ykeylist[1]-ykeylist[0])*0.5

        D2D = []
        RI2D = []
        OL2D = []
        DError2D = []
        RIError2D = []
        OLError2D = []

        for j in range(len(YKeys)):
            D2D.append([])
            RI2D.append([])
            OL2D.append([])
            for i in range(len(XKeys)):
                D2D[j].append(0)
                RI2D[j].append(0)
                OL2D[j].append(0)

        if (len(DError) > 0):
            for j in range(len(YKeys)):
                DError2D.append([])
                RIError2D.append([])
                OLError2D.append([])
                for i in range(len(XKeys)):
                    DError2D[j].append(0)
                    RIError2D[j].append(0)
                    OLError2D[j].append(0)

        #test each data and fill up the matrix
        for i in range(len(D)):
            yindex = len(YKeys) - 1- YKeys["%.3f" % Y[i]]
            xindex = XKeys["%.3f" % X[i]]
            D2D[yindex][xindex] = D[i]
            RI2D[yindex][xindex] = RI[i]
            OL2D[yindex][xindex] = OL[i]

        for i in range(len(D)):
            yindex = len(YKeys) - 1- YKeys["%.3f" % Y[i]]
            xindex = XKeys["%.3f" % X[i]]
            D2D[yindex][xindex] = D[i]
            RI2D[yindex][xindex] = RI[i]
            OL2D[yindex][xindex] = OL[i]
            if (len(DError) > 0):
                DError2D[yindex][xindex] = DError[i]
                RIError2D[yindex][xindex] = RIError[i]
                OLError2D[yindex][xindex] = OLError[i]

        ### histogram ###

        (HistoXD, HistoYD) = self.HistogramXY(D, 7, 0.1, 1)
        (HistoXRI, HistoYRI) = self.HistogramXY(RI, 7, 0.0001, 3)
        (HistoXOL, HistoYOL) = self.HistogramXY(OL, 7, 0.1, 1)

        if (params.field_results['useheader'][:1] == '1'):
            Prefix = CassetteName + '-'
        elif (params.field_results['useheader'][:1] == '2'):
            if (params.field_results['header'] != ''):
                Prefix = params.field_results['header'] + '-'

        result = [ImportDatasetText(Prefix + 'name', [CassetteName]),
                  ImportDataset1D(Prefix + 'x', X),
                  ImportDataset1D(Prefix + 'y', Y),
                  ImportDataset1D(Prefix + 'd', D),
                  ImportDataset1D(Prefix + 'n', RI),
                  ImportDataset1D(Prefix + 'ol', OL),
                  ImportDataset2D(Prefix + 'd2D', D2D, (XMin, XMax), (YMin, YMax)),
                  ImportDataset2D(Prefix + 'n2D', RI2D, (XMin, XMax), (YMin, YMax)),
                  ImportDataset2D(Prefix + 'ol2D', OL2D, (XMin, XMax), (YMin, YMax)),
                  ImportDatasetText(Prefix + 'd-histo_x', HistoXD),
                  ImportDataset1D(Prefix + 'd-histo_y', HistoYD),
                  ImportDatasetText(Prefix + 'n-histo_x', HistoXRI),
                  ImportDataset1D(Prefix + 'n-histo_y', HistoYRI),
                  ImportDatasetText(Prefix + 'ol-histo_x', HistoXOL),
                  ImportDataset1D(Prefix + 'ol-histo_y', HistoYOL),
                  ImportDatasetText(Prefix + 'd-avg', ["%.3f nm" % mean(D)]),
                  ImportDatasetText(Prefix + 'd-stdev', ["%.3f nm" % std(D)]),
                  ImportDatasetText(Prefix + 'd-med', ["%.3f nm" % median(D)]),
                  ImportDatasetText(Prefix + 'd-q1', ["%.3f nm" % percentile(D, 25)]),
                  ImportDatasetText(Prefix + 'd-q3', ["%.3f nm" % percentile(D, 75)]),
                  ImportDatasetText(Prefix + 'n-avg', ["%.4f" % mean(RI)]),
                  ImportDatasetText(Prefix + 'n-stdev', ["%.4f" % std(RI)]),
                  ImportDatasetText(Prefix + 'n-med', ["%.4f" % median(RI)]),
                  ImportDatasetText(Prefix + 'n-q1', ["%.4f" % percentile(RI, 25)]),
                  ImportDatasetText(Prefix + 'n-q3', ["%.4f" % percentile(RI, 75)]),
                  ImportDatasetText(Prefix + 'ol-avg', ["%.3f nm" % mean(OL)]),
                  ImportDatasetText(Prefix + 'ol-stdev', ["%.3f nm" % std(OL)]),
                  ImportDatasetText(Prefix + 'ol-med', ["%.3f nm" % median(OL)]),
                  ImportDatasetText(Prefix + 'ol-q1', ["%.3f nm" % percentile(OL, 25)]),
                  ImportDatasetText(Prefix + 'ol-q3', ["%.3f nm" % percentile(OL, 75)])]
        if ( len(DError) > 0):
            result += [ImportDataset1D(Prefix + 'd_err', DError),
                       ImportDataset1D(Prefix + 'n_err', RIError),
                       ImportDataset1D(Prefix + 'ol_err', OLError),
                       ImportDataset2D(Prefix + 'd_err2D', DError2D, (XMin, XMax), (YMin, YMax)),
                       ImportDataset2D(Prefix + 'n_err2D', RIError2D, (XMin, XMax), (YMin, YMax)),
                       ImportDataset2D(Prefix + 'ol_err2D', OLError2D, (XMin, XMax), (YMin, YMax)),
                       ImportDatasetText(Prefix + 'd_err-avg', ["%.3f nm" % mean(DError)]),
                       ImportDatasetText(Prefix + 'd_err-stdev', ["%.3f nm" % std(DError)]),
                       ImportDatasetText(Prefix + 'd_err-med', ["%.3f nm" % median(DError)]),
                       ImportDatasetText(Prefix + 'd_err-q1', ["%.3f nm" % percentile(DError, 25)]),
                       ImportDatasetText(Prefix + 'd_err-q3', ["%.3f nm" % percentile(DError, 75)]),
                       ImportDatasetText(Prefix + 'n_err-avg', ["%.3f" % mean(RIError)]),
                       ImportDatasetText(Prefix + 'n_err-stdev', ["%.3f" % std(RIError)]),
                       ImportDatasetText(Prefix + 'n_err-med', ["%.3f" % median(RIError)]),
                       ImportDatasetText(Prefix + 'n_err-q1', ["%.3f" % percentile(RIError, 25)]),
                       ImportDatasetText(Prefix + 'n_err-q3', ["%.3f" % percentile(RIError, 75)]),
                       ImportDatasetText(Prefix + 'ol_err-avg', ["%.3f nm" % mean(OLError)]),
                       ImportDatasetText(Prefix + 'ol_err-stdev', ["%.3f nm" % std(OLError)]),
                       ImportDatasetText(Prefix + 'ol_err-med', ["%.3f nm" % median(OLError)]),
                       ImportDatasetText(Prefix + 'ol_err-q1', ["%.3f nm" % percentile(OLError, 25)]),
                       ImportDatasetText(Prefix + 'ol_err-q3', ["%.3f nm" % percentile(OLError, 75)])]
        return result
            
importpluginregistry.append(Sentech_ASCII())
        
