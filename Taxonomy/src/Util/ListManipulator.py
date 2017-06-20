print 'Import ListManipulator.py'
import FASTA as FA

from string import join
#import ListManipulator as LM
import math
import copy
import random
DEBUG = None
location = {}
Verbose = None




def Normalize(inLst, flipOrder = None, newMaxMinRange = (None, None), Verbose = 0):
    ##Raise
    newMax = max(newMaxMinRange)
    newMin = min(newMaxMinRange)
    
    Min = min(inLst)
    if Min < 0:
        RaiseAboveZero_Offset = Min
        adj_inLst = [x-RaiseAboveZero_Offset for x in inLst]
    else:
        adj_inLst = inLst
        RaiseAboveZero_Offset = 0
    
    ##CalculateMaxMin
    Max = max(adj_inLst)
    Min = min(adj_inLst)
    maxminusmin = Max - Min
    ##Normal
    NormalLst = [(x-Min)/maxminusmin for x in adj_inLst]
    ##flip
    if flipOrder:
        NormalLst = [1-x for x in NormalLst]
    ##RestoreRange
    if newMin !=None and newMax!=None and newMax > newMin:
        print 'yo'
        RaiseAboveZero_Offset = 0
        Newmaxminusmin = newMax - newMin
    else:
        newMax = Max
        newMin = Min
        Newmaxminusmin = maxminusmin
        
    outLst = [x*Newmaxminusmin + newMin + RaiseAboveZero_Offset for x in NormalLst]
    if Verbose:
        print "\nNormalize: \nflipOrder %s\nmax %.2f newMax %.2f\nmin %.2f newMin %.2f" %(flipOrder, Max, newMax, Min, newMin)
        print "Normalize Normalized a list of %i elements: "%len(inLst)
        for v in range(min(len(inLst), 50)):
            print "%i\t%.2f\t%.2f" %(v, inLst[v], outLst[v])
        print "\n\n"
    return outLst




def SortLst(inLst, CompareFunctionLst = []):
    cmpObj = MultiCompare(CompareFunctionLst)
    inLst.sort(lambda x,y:cmpObj.compare(x,y))

class MultiCompare:
    ##Usage:
    #cmpObj = MultiCompare([1,2,lambda x:x.seq])
    #probeDB.sort(lambda x,y:cmpObj.cmp(x,y))
    def __init__(self, compLst = []):
        self.compLst = [DataAccessFunctionProcessor(indexer) for indexer in compLst]
    def compare(self, x, y):
        for indexer in self.compLst:
            if indexer(x) > indexer(y):
                return 1
            elif indexer(x) < indexer(y):
                return -1
        return 0


def UseDontUse(instr):
    if not(instr): return None
    if 'str' in str(type(instr)):
        if instr == '0': return None
        if instr == 'None': return
    return 1

def NoNone(Lst):
    def nonone(val):
        if val <> None:
            return 1
    return filter(nonone, Lst)


def genomeChoice(location, locLst,  chrLst = None,  returnChooseIndices = None):
    #assume at least  [genomicLocation] [chromosome positions] chromosme positions are optional assumes lists are  sorted and indexed to each other
    #This program runs Interval chooser it return s singl col of data either the indices or a list of 1,0's
    #Usage: choiceLst = LM.genomeChoice(location, LM.extractCol(indb[1:], 1), chrLst = LM.extractCol(indb[1:], 0), returnChooseIndices = None)
    if not(chrLst): chrLst = ['' for i in range(len(locLst))]
    IntraChromosomalAddVal = 10000
    ContinuousLocLst = ['loc']
    lastChromosome = ''
    lastLoc = 0
    incrementalLoc = 1
    for i in range(len(chrLst)):
        chromosome = chrLst[i]
        loc = int(locLst[i])
        if chromosome != lastChromosome:
            addVal = IntraChromosomalAddVal
            lastChromosome = chromosome
        else:
            addVal = loc - lastLoc
            lastLoc = loc
            incrementalLoc = incrementalLoc + addVal
        ContinuousLocLst.append(incrementalLoc)
    chooseLst = IntervalChoice(ContinuousLocLst[1:],  2400, 5000, location, root = 'genericIntervalChoice')
    if returnChooseIndices:
        return chooseLst
    else:
        outLst = [0 for i in range(len(locLst))]
        for i in chooseLst:
            outLst[i] =1
        return outLst


def ProbeChoice(probeLst, locationAccessFunction, NumProbesToTake, location, HasLL = None, PreSort = 1, root = 'genericIntervalChoice'):
    #Usage: probedb = ProbeChoice(probeLst, locationAccessFunction, NumProbesToTake, location, HasLL = None, root = 'genericIntervalChoice')
    windowSize = 5000

    if NumProbesToTake > len(probeLst): return probeLst
    
    if HasLL: LL = probeLst.pop(0)
    if PreSort: ##sometimes you may not want this function to sort your data ahead of time
        SortLst(probeLst, CompareFunctionLst = [locationAccessFunction])
    locLst = [locationAccessFunction(pobj) for pobj in probeLst]
    
    chosenIndexLst = IntervalChoice(locLst, NumProbesToTake, windowSize, location, root = root,  Report = None)
    outLst = [probeLst[i] for i in chosenIndexLst]
    if HasLL: outLst = [LL] + outLst
    return outLst



def IntervalChoice(locLst,  probesToChoose, windowSize, location, root = 'genericIntervalChoice', Report = None):
    ##Usage: outLst = IntervalChoice(locLst,  probesToChoose, windowSize, root = 'genericIntervalChoice')
    ## where locLst is [['42342341'],['42342390'], . . .] sorted !!!! or ['1232344', '1232390', . . .]
    UD = location['WorkDir']
    if probesToChoose > len(locLst):
        print '/n/n/n/n IntervalChoice: Submitted More requested probes %i than available probes %i \n\n Will Not generate File: %s \n\n\n' %(probesToChoose, len(locLst), 'ChosenProbes_%s_T.txt' %root)
        return [i for i in range(probesToChoose)]
    weightLst, progressiveSum = WeighLocList(locLst, windowSize = windowSize)
    if DEBUG:
        ofn = 'weighedLst_startPoint_T.txt' %root
        OFN = FA.path([UD,ofn],fileType = 'tab')
        OFN.outDB(weightLst)
    
    #probesChosen = 0; ChosenIndices = [0 for i in range(len(locLst))]
    weightLst_LL = ['index', 'Loc_Num', 'startIndex', 'endIndex', 'weight_4', 'extraWeight_5', 'progressiveSum_6', 'presence_7', 'chose_8']
    cycle = 0; nextcycle = cycle + 100
    probesChosen = 0
    trobesTakenInCycle = 100
    while probesChosen  < probesToChoose:
        cycle +=1
        ##Recalculate Interval function
        intervalLst = []
        progressiveSum = 0.0
        for i in range(len(weightLst)):
            x = weightLst[i]
            weight = x[4]
            extraWeight = x[5]
            presence = x[7]
            if presence:
                finalWeight = weight + extraWeight
                if finalWeight > 0:
                    finalWeight = 1.0/finalWeight ##<< this is the only place where inversing takes place
                else:
                    finalWeight = 0
            else:
                finalWeight = 0
            progressiveSum = progressiveSum + finalWeight ## interval = 1/(0 or weight + additionalWeight)*present when chosen, presesent becomes 0
            intervalLst.append(progressiveSum)
            weightLst[i][6] = progressiveSum
        ##Choose randomly
        randNum = random.random()*progressiveSum
        ChosenIndex = binsearch(randNum, 0, len(intervalLst), intervalLst)+1 ##if you don't add this it will always hit the wrong probe
        ##Add probe to chosen list
        weightLst[ChosenIndex][8] = 1
        weightLst[ChosenIndex][7] = 0
        ##Return weights back to weightLst using weight function
        startIndex, endIndex = (weightLst[ChosenIndex][2], weightLst[ChosenIndex][3])
        posMatrix = generatePositionalMatrix(locLst, ChosenIndex, startIndex, endIndex, windowSize)#Returns [[i, weight, Loc]]
        for posobj in posMatrix:
            i, weight, Loc = posobj
            weightLst[i][5] = weightLst[i][5] - weight
        probesChosen = SumArr(weightLst,8)
        if cycle > nextcycle:
            nextcycle = cycle + 100
            print 'cycle %i ProbesChosen %i Probes left to choose %i progressiveSum: %.2f'%(cycle, probesChosen, probesToChoose-probesChosen,progressiveSum)

    weightLst = [weightLst_LL] + weightLst
    if Report:
        ofn = 'ChosenProbes_%s_T.txt' %root
        OFN = FA.path([UD,ofn],fileType = 'tab')
        OFN.outDB(weightLst)
    outLst = []
    for i in range(1,len(weightLst),1):
            if weightLst[i][8] ==1:
                outLst.append(i-1)##weightLst has a label line on top of itself at this point
    return outLst ##This is a list of indices that can be used to access the chosen probes from the original probeLst


def getIndexForWindowSize(locLst, i, windowSize):
    #usage startIndex, endIndex = getIndexForWindowSize(locLst, i, windowSize)
    numElements = len(locLst)
    Loc_Num = locLst[i]
    Loc_start = max(locLst[0], Loc_Num - windowSize)
    Loc_end = min(locLst[-1], Loc_Num + windowSize)
    startIndex = i
    Loc = locLst[startIndex]
    while Loc > Loc_start:
        startIndex = max(startIndex-1, 0)
        Loc = locLst[startIndex]
    if Loc < Loc_start: startIndex +=1
    endIndex = i
    Loc = locLst[endIndex]
    while Loc < Loc_end:
        endIndex = min(endIndex+1, numElements-1)
        Loc = locLst[endIndex]
        if endIndex == numElements-1: break
    if Loc > Loc_end: endIndex = min(endIndex-1, numElements)
    return startIndex, endIndex

def generatePositionalMatrix(locLst, homeIndex, startIndex, endIndex, windowSize):
    ## this generates a set of weights to be applied to numbers near a given position
    posMatrix = []
    if startIndex == endIndex:
        posMatrix = [[homeIndex, 1.0, locLst[homeIndex]]]
    for i in range(startIndex, endIndex+1, 1):##This had to be fixed by adding 1 to endIndex otherwise it consistently missed the end
        Loc = locLst[i]
        delta = float(abs(Loc - locLst[homeIndex]))
        weight = max(1.0 - delta/float(windowSize), 1/float(windowSize))##Zero is not allowed !!
        ## the 1- assures that probes that are closer get more value than probes that are farther away
        posMatrix.append([i, weight, Loc])
    
    return posMatrix


def WeighLocList(locLst, windowSize = 50000):
    #Usage: weighLst = WeighLocList(locLst, windowSize = 50000)
    numElements = len(locLst)
    if 'list' in str(type(locLst[0])):
        for i in range(numElements):
            locLst[i] = int(locLst[i][0])            
    weightLst = []
    t = US.ReportTimer(numElements, ReportInterval = 1000, ProcessName = 'Weight Generator')
    progressiveSum = 0.0
    for elemIndex in range(numElements):
        ret = t(elemIndex)
        Loc_Num = locLst[elemIndex]
        startIndex, endIndex = getIndexForWindowSize(locLst, elemIndex, windowSize)
        calcLst = [];testLst = []
        posMatrix = generatePositionalMatrix(locLst, elemIndex, startIndex, endIndex, windowSize)#Returns [[i, weight, Loc]]
        if DEBUG: print 'startIndex %i endIndex %i Delta %i' %(startIndex, endIndex, abs(posMatrix[0][-1] - posMatrix[-1][-1]))
        weight = max(SumArr(posMatrix,1),1)#doing this avoids a div by zero error 091908
        #print str(posMatrix)
        progressiveSum += 1.0/weight
        #[index, Loc_Num, startIndex, endIndex, weight, extraWeight, progressiveSum, presence]
        weightLst.append([elemIndex, Loc_Num, startIndex, endIndex, weight, 0.0, progressiveSum, 1.0, 0])
    return weightLst, progressiveSum



def binsearch(val, start, end, inlst):
    pivot = int((end+start)/2)
    tval = inlst[pivot]
    #print 'binsearch %.2f tval %.2f %i %i %i' %(val, tval, start, pivot, end)
    if pivot == start: return pivot
    if val >= tval:
        b = binsearch(val, pivot, end, inlst)
    else:
        b = binsearch(val, start, pivot, inlst)
    return b


class Dictionator:
    def __init__(self, queryStr, inLst):
        #Usage:
        # Lstdict = Dictionator('rsID',inLst)
        # object = Lstdict['rs63542']
        if 'int' in str(type(queryStr)):
            execStr = "self.queryFunction = lambda x: x[%i]" %queryStr
        else:
            execStr = "self.queryFunction = lambda x: x.%s" %queryStr
        exec execStr
        self.inLst = inLst
        self.build()
    def build(self):
        print 'Dictionator: Building lists'
        self.captureLstState()
        self.generateDict()
    def captureLstState(self):
        self.first = self.queryFunction(self.inLst[0])
        self.last = self.queryFunction(self.inLst[-1])
        self.n = len(self.inLst)
    def checkLst(self):
        if not(self.first == self.queryFunction(self.inLst[0])) or not(self.last == self.queryFunction(self.inLst[-1])):
            print 'Dictionator: order has changed'
            return None
        if not(self.n == len(self.inLst)):
            print 'Dictionator: elements added'
            return None
        return 1
    def generateDict(self):
        self.vdict = {}
        for i in range(len(self.inLst)):
            self.vdict[self.queryFunction(self.inLst[i])] = i
    def __getitem__(self,val):
        if val in self.vdict:
            i = self.vdict[val]
            return self.inLst[i]
        else:
            self.build()
            if val in self.vdict:
                i = self.vdict[val]
                return self.inLst[i]
            else:
                raise 'Dictionator: %s NOT FOUND' %str(val)
                
Verbose = 0
class DataGrid:
    def __init__(self, emptyCellValue = ''):
        self.emptyCellValue = emptyCellValue
        self.grid = [[]]
        self.MeasureGridSize()
    def reset(self):
        self.grid = [[]]
        self.MeasureGridSize()
    def display(self):
        for r in range(self.NumRows):
            if Verbose: print '%i\t\t->%s' %(r, '\t'.join([str(obj) for obj in self.grid[r]]))
    def Add_ColLst(self, inLst):
        #[[ColData], [ColData]]
        for col in inLst:
            self.AddSingleCol(col)
        self.ResizeGrid()
    def Add_RCLst(self, inLst):
        if Verbose: print 'Add_RCLst:'
        #RCLst = [[rowVal1, rowVal2 . . . ]]
        NumRows, NumCols = self.MeasureGridSize(inLst)##This also Lists free objects
        self.ResizeGrid(NumRows = NumRows, NumCols = None)
        inCols = self.ConvertToCols(inLst)
        for Col in inCols:
            self.AddSingleCol(Col)
        self.ResizeGrid()
    def AddRowsToBottom(self, inLst):
        NumRows, NumCols = self.MeasureGridSize(inLst)##This also Lists free objects
        self.ResizeGrid(NumRows = None, NumCols = NumCols)
        for row in inLst:
            self.grid.append(row)
        self.ResizeGrid()
    def AddSingleCol(self, inCol):
        if Verbose: print 'AddSingleCol:'
        self.ResizeGrid(NumRows = len(inCol))
        for r in range(len(inCol)):
            self.grid[r].append(inCol[r])
    def ResizeGrid(self, NumRows = None, NumCols = None):
        self.ResizeLst(self.grid, NumRows = NumRows, NumCols = NumCols)
        self.MeasureGridSize()
    def ResizeLst(self, inLst, NumRows = None, NumCols = None):
        if Verbose: print 'ResizeLst:'
        if NumRows == None:
            NumRows = self.NumRows
        if NumCols == None:
            NumCols = self.NumCols
        if Verbose: print 'ResizeGrid: from (%i, %i) to (%i, %i)'%(self.NumRows, self.NumCols, NumRows, NumCols) 
        for r in range(NumRows):
            while len(inLst)< r+1:
                if Verbose: print 'Adding row %i to %i' %(len(inLst), NumRows)
                inLst.append([])
            while len(inLst[r]) < NumCols:
                if Verbose: print 'Adding col %i  to %i' %(len(inLst[r]), NumCols)
                inLst[r].append(self.emptyCellValue)
        return inLst
    def ConvertToCols(self, inLst):
        if Verbose: print 'ConvertToCols:'
        NumRows, NumCols = self.MeasureGridSize(inLst)
        inLst = self.ResizeLst(inLst, NumRows, NumCols)
        ColLst = self.Transverse(inLst, NumRows = NumRows, NumCols = NumCols)
        return ColLst
    def Transverse(self, inLst, NumRows = None, NumCols = None):
        if Verbose: print 'Trnasverse:'
        if NumCols == None:
            NumCols = self.CountCols(inLst)
        if NumRows == None:
            NumRows = self.CountRows(inLst)
        OutLst = [[self.emptyCellValue for r in range(NumRows)] for c in range(NumCols)]
        for c in range(NumCols):
            for r in range(NumRows):
                OutLst[c][r] = inLst[r][c]
        return OutLst
    def MeasureGridSize(self, inLst = None):
        if Verbose: print 'MeasureGridSize:'
        if inLst == None:
            self.NumRows, self.NumCols = self.MeasureGridSize(inLst = self.grid)
            return self.NumRows, self.NumCols
        if len(inLst) == 0:
            return 0,0
        for i in range(len(inLst)):
            if not('list' in str(type(inLst[i]))):
                print 'List this: ' + str(inLst[i])
                inLst[i] = [inLst[i]]
        return self.CountRows(inLst), self.CountCols(inLst)
    def CountCols(self, inLst):
        colLst = [len(obj) for obj in inLst]
        return max(colLst)
    def CountRows(self, inLst):
        return len(inLst)
    
def MergeCols(aLst, bLst):
    colsToAdd = len(bLst[0])
    for i in range(len(aLst)):
        aLst[i] = aLst[i] + bLst[i]
    return aLst

def AddCol(Lst, addcol):
    """AddCol(Lst, addcol)|
    adds col to list of rows fills in with empty space||"""
    NumCols = len(Lst[0])
    startLen = len(Lst)
    addcolLen = len(addcol)
    if addcolLen > startLen:
        for i in range(addcolLen-startLen):
            Lst.append([''for i in range(NumCols)])
    j = 0
    for i in range(len(Lst)):
        if j < len(addcol):
            Lst[i].append(addcol[j])
            j+=1
        else:
            Lst[i].append('')
    return Lst

def alt_AddCol(Lst, addcol):
    """AddCol(Lst, addcol)|
    adds col to list of rows fills in with empty space||"""
    NumCols = len(Lst[0])
    startLen = len(Lst)
    addcolLen = len(addcol)
    if addcolLen > startLen:
        for i in range(addcolLen-startLen):
            Lst.append([''for i in range(NumCols)])
    j = 0
    for i in range(len(Lst)):
        if j < len(addcol):
            #print addcol[j]
            addObj = addcol[j]
            if 'list' in str(type(addObj)):
                Lst[i].extend(addObj)
            else:
                Lst[i].append(addObj)
            j+=1
        else:
            Lst[i].append('')
    return Lst

##def extractCol(indb, index):
##    """extractCol(indb, index)|||"""
##    outdb = []
##    for obj in indb:
##        if len(obj) > index:
##            outdb.append(obj[index])
##    return outdb

def extractCol(indb, index):
    """extractCol(indb, index)|||"""
    indexerF = DataAccessFunctionProcessor(index)
    outdb = []
    for obj in indb:
        if len(obj) > index:
            outdb.append(indexerF(obj))
    return outdb

def extract_Cols(indb, start = 0, end = None, exCols = []):
    if not(end):
        end = len(indb)
    outdb = []
    maxCol = max(exCols)
    for i in range(start, end, 1):
        outLine = []
        for col in exCols:
            if col < len(indb[i]):
                outLine.append(indb[i][col])
            else:
                outLine.append('nd')
        outdb.append(outLine)
    return outdb




def NormalizeColumn(inLst, NormalizeTo = 1.0, startRow = 0, colIndex = None):
    """NormalizeColumn(inLst, NormalizeTo = 1.0, startRow = 0, colIndex = None)|
    Normalizes A column to a set number, sum, avg, or median
    colIndex can be int or a number"""
    indb = copy.deepcopy(inLst)
    functSpace = {}
    if colIndex ==None:
        valFunct = lambda x:float(x)
        def access(indb, index, val):
            indb[index] = val
    else:
        cmd = 'valFunct = lambda x: float(x[' + str(colIndex) + '])'
        exec(cmd) in functSpace
        def access(indb, index, val):
            indb[index][colIndex] = val
    if 'str' in str(type(NormalizeTo)):
        if NormalizeTo.lower() =='sum': NormalizeTo = SumArr(inLst[startRow:])
        if NormalizeTo =='avg': NormalizeTo = AvgArr(inLst[startRow:])
        if NormalizeTo =='median': NormalizeTo = MedianArr(inLst[startRow:])
    allvals = []
    for i in range(startRow, len(indb), 1):
        val = valFunct(indb[i])
        allvals.append(val)
    sumOfLine = sum(allvals)
    for i in range(startRow, len(indb), 1):
        if sumOfLine == 0:
            normalizedVal = 0.0
        else:
            normalizedVal = float(valFunct(indb[i]))/float(NormalizeTo)
        access(indb, i, normalizedVal)
    return indb
    
class LLDict:
    """Very uselful little utility for using Label Lines in col data"""
    def __init__(self, inLst):
        self.Ldict = {}
        self.inLst = inLst
        for l in range(len(inLst)):
            label = inLst[l]
            self.Ldict[label] = l
    def __getitem__(self, query):
        return self.Ldict[query]
    def processLine(self, line, queryLst):
        outLst = []
        for query in queryLst:
            if 'int' in str(type(query)):
                outLst.append(line[query])
            elif not(query in self.Ldict):
                outLst.append(str(query)+'_MISSING')
            else:
                outLst.append(line[self.Ldict[query]])
        return outLst
    def getitems(self, queryLst, line):
        outLine = []
        for query in queryLst:
            outLine.append(self.getitem(query, line))
        return outLine
    def getitem(self, query, line):
        return line[self[query]]
    def __call__(self, query, line):
        return line[self[query]]
    def display(self):
        for i in range(len(self.inLst)):
            print '%i\t%s' %(i, self.inLst[i])


def strRepresentList(inLst = None):
    if 'list' in str(type(inLst)):
        ##Convert list to string
        if not('][' in str(inLst)):
            for i in range(len(inLst)):
                inLst[i] = str(inLst[i])
            outStr = ']['.join(inLst)
            return outStr
        else:
            raise 'strRepresentList got a list with ][ in it. Can\'t do anything with it'
            return None
               
    elif 'str' in str(type(inLst)):
        ##convert from str to list
        return inLst.split('][')
    else:
        print 'strRepresentList got neither a list nor a string' 
        return None

def PercentCutoffTallyHo(inList, sortIndex, percCutoff):
    LBU = TallyHo(inList, sortIndex)
    LBU.sort(lambda x,y:cmp(x[-1], y[-1]))
    LBU.reverse()
    sumElementTotal = 0.0
    NumElements = float(SumArr(LBU,1))
    percCutoff = float(percCutoff)
    
    elementsToTake = NumElements * percCutoff/100.0
    objsChosen = []
    elementsTaken = 0
    for i in range(len(LBU)):
        objsChosen.append(LBU[i][0])
        elementsTaken = elementsTaken + LBU[i][1]
        percTaken = int(100.0*elementsTaken/NumElements)
        LBU[i].extend([elementsTaken,percTaken]) 
        if elementsTaken > elementsToTake: break
    percTaken = int(100.0*elementsTaken/NumElements)
    print 'Took %i elements %i groups which is %i percent of population' %(elementsTaken, len(objsChosen), percTaken)
    return objsChosen, LBU



class LabelListObj:
    """This keeps track of Large lists with Label lines at the top"""
    def __init__(self,LLst=None, inDB=[]):
        self.BIGFILE = None
        if 'list' in str(type(inDB)):
            self.inDB = inDB
        else:
            if 'FASTA.path' in str(inDB):
                self.inDBPATH = inDB
                if inDBPATH.getsize() < 100000000:
                    self.inDB = inDBPATH.inDB()
                else:
                    self.BIGFILE = 1
                    self.inDB = inDBPATH.TabChunk(20000)
                    print 'Incomplete file input'
            else:
                self.inDBPATH = None
                
                self.inDB
        if not(LLst):
            LLst = self.inDB.pop(0)
            #if this fails it means you sent to empty objects to this function
        self.LL = LLst
        self.Dict = {}
        self.MakeDict()
    def MakeDict(self):
        self.Dict = {}
        for j in range(len(self.LL)):
            label = self.LL[j]
            if label in self.Dict:
                print 'LLobj: Duplicate Label: ' + label
            self.Dict[label] = j
    def AddElemnts(self,inLst):
        for obj in inLst:
            if not(obj in self.LL):
                self.LL.append(obj)
        self.MakeDict()
    def __getitem__(self,index):
        return self.Dict[index]
    def filterDB(self, indexLst = None, LabelLst = None):
        if indexLst == None and LabelLst == None:
            retObj = [[self.LL]] + self.inDB
            return retObj
        if not(indexLst):
            indexLst = []
            for obj in LabelLst:
                indexLst.append(self[LabelLst])
        outDB = []
        for obj in self.inDB:
            listObj = []
            for i in indexLst:
                listObj.append(obj[i])
            outDB.append(listObj)
        return outDB
    def element(self,inObj,Label):
        index = self[Label]
        return inObj[index]
    
    
        
def printList(inlst, LL= None):
    inLst = copy.deepcopy(inlst)
    if LL and LL.isdigit():
        LL = inLst.pop(LL)
    if 'list' in str(type(LL)):
        labelLine = LL
    else:
        labelLine = []
    inLst = labelLine + inLst
    printStr = ''
    for line in inLst:
        piece = ''
        for obj in line:
            if not(piece):
                piece = obj
            else:
                piece = piece +'\t'+ str(obj)
        piece = piece + '\n'
        printStr = printStr + piece
    print printStr
    
    
def ExcelCI(colstr, INobj = None):
    """ExcelCI: Excell column name to index converter"""
    index = 0
    if len(colstr) >2:
        print colstr + ' is too long'
        return None
    alphabet = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    partA = colstr[0].lower()
    
    if partA in alphabet:
        partA = alphabet.index(partA)
    else:
        print 'BAD partA' + partA
        return None
    if len(colstr) ==2:
        partB= colstr[1].lower()
    else:
        partB = 0
    if partB:
        if partB in alphabet:
            partB = alphabet.index(partB)
        else:
            print 'BAD partB' + partB
            return None
        index = (partA+1) * 26
        index = index + partB
    else:
        index = partA
    if not(INobj):
        return index
    else:
        return INobj[index]


def MakeDictFromList(indb,fromIndex,toindex):
    """MakeDictFromList
     |(indb,fromIndex,toindex)
     |this is self-explanatory
     ||"""
    outDict = {}
    for obj in indb:
        outDict[obj[fromIndex]] = obj[toindex]
    return outDict




def Rank_deprecated(IFNorLst, UD, OFN = None, groupIndex = 0, scoreIndex = 1, reverse = None, RankOutIndex = 2):
     """Rank
     |(IFNorLst,UD,OFN = None,groupIndex = 0,scoreIndex = 1,reverse = None, RankOutIndex = 2)
     |FN or Ranked List
     |It is FN or ListAaware|"""
     
     if str(type(IFNorLst)).find('str'):
             IFN = IFNorLst
             inDB = InArray(IFN,UD,separator ='\t')
     else:
          inDB = IFNorLst
     
     
     inDB_BU = FastListBreakup(inDB,groupIndex)
     outLst = []
     
     for obj in inDB_BU:
         obj.sort(lambda x,y: cmp(float(x[scoreIndex]),float(y[scoreIndex])))
         if reverse:
             obj.reverse()
             
         rankIncr = 1.0/len(obj)
         rank = rankIncr
         
         for item in obj:
             if len(item) <= RankOutIndex:
                 item.append(rank)
             else:
                 item[RankOutIndex] = rank
             rank += rankIncr
             outLst.append(item)
     if OFN:   
         ofn = OutArray(outLst,OFN,UD)
         return ofn
     else:
         return outLst


def Rank(inLst):
    #if a value is None, it gets converted to 0
    #Rank starts at 1
    #inLst [0, .2, None, 3, 4, 5, -0.25]floats and None's
    wLst = []#[[index, value, rank]]
    for v in range(len(inLst)):
        val = inLst[v]
        if val <> None:
            val = float(val)
        wLst.append([v, val, 0])
    wLst.sort(lambda x,y:cmp(x[1], y[1]))
    rank = 0
    NumInsignificantValues = 0
    for v in range(len(wLst)):
        val = wLst[v][1]
        rank +=1#You must increment this each timer otherwise comparison will be totally lost
        if val <> None:
            wLst[v][2] = rank
        else:
            NumInsignificantValues +=1
    insignificantRankValue = 1.0/float(NumInsignificantValues)
    
    wLst.sort(lambda x,y:cmp(x[0], y[0]))#resort inLst to get Ranks out
    for v in range(len(wLst)):##Deal with ties at the bottom
        if wLst[v][2] == 0:
            wLst[v][2] = insignificantRankValue
    outLst = [x[2] for x in wLst]
    print 'Rank: %i  %i insignificant values' %(rank, NumInsignificantValues)
    return outLst


def MakeBinThresh(inLst, numBins):
    wLst = []
    wLstNone = []
    numVals = 0
    for v in range(len(inLst)):
        val = inLst[v]
        if val < 0:
            val = None
        if not(val):
            val = None
        if val:
            numVals +=1
            wLst.append([v,val,0])
        else:
            wLstNone.append([v,val,0])
    intervalSize = int(float(numVals)/float(numBins))
    bin = 1;nexti = intervalSize
    wLst.sort(lambda x,y:cmp(x[1], y[1]))
    for i in range(len(wLst)):
        if i > nexti:
            nexti = i + intervalSize
            bin +=1
        wLst[i][-1] = bin
    oLst = wLst + wLstNone
    oLst.sort(lambda x,y:cmp(x[0], y[0]))
    outLst = [valObj[-1] for valObj in oLst]
    return outLst
 

def QuantBin(SignalLst, NumQuants = 5):
    qbinLst = []
    numElements = len(SignalLst)
    segmentSize = int(numElements/NumQuants)
    positions = [i for i in range(0, numElements+seqmentSize, segmentSize)]
    positions[-1] = numElements
    valLst = []
    for i in range(len(SignalLst)):
        valLst.append([i, float(SignalLst[i]), 0])
    valLst.sort(lambda x,y: cmp(x[1], y[1]))
    compareIndex = 0
    for p in range(0, NumQuants, 1):
        start = positions[p]
        end = positions[p+1]
        for i in range(start, end, 1):
            valLst[i][2] = p
    valLst.sort(lambda x,y: cmp(x[0], y[0]))
    outLst = [x[2] for x in valLst]
    return outLst


def Correl(X,Y):
    """Correl(X,Y)| please fillininfohere|need two paired lists"""
    ##CheckInput
    xLen = len(X)
    yLen =len(Y)
    OK = None
    if yLen == xLen:
        OK = 1
    if yLen > 0 and OK:
        OK = 1
    if not(OK):
        print 'Correl Sez: ERRRRROR Xlen: %i Ylen %i ' %(xLen,yLen)
        return 'ERRROR'
    #Calculate XY
    sumXY = 0
    sumX = 0
    sumXsqrd = 0
    sumY = 0
    sumYsqrd = 0
    N = 0#xLen
    for i in range(xLen):
        x = X[i]
        y = Y[i]
        if x == None or y == None: continue
        N+=1
        x = float(x)
        y = float(y)
        ###print str(x)+ ' and ' + str(y)
        sumX += x
        sumXsqrd += pow(x,2)
        sumY += y
        sumYsqrd += pow(y,2)
        sumXY += x * y
    if N == 0: return None
    rTop = sumXY - (sumX*sumY/N)
    rBottom = math.sqrt((sumXsqrd-pow(sumX,2)/N)*(sumYsqrd-pow(sumY,2)/N))
    if rBottom == 0: return 0
    correl = rTop/rBottom
    return correl

def ThresholdCorrel(X,Y, testFunction = None):
    """Correl(X,Y)| please fillininfohere|need two paired lists"""
    if not(testFunction):
        exec 'testFunction = lambda x:1'
    ##CheckInput
    xLen = len(X)
    yLen =len(Y)
    OK = None
    if yLen == xLen:
        OK = 1
    if yLen > 0 and OK:
        OK = 1
    if not(OK):
        print 'Correl Sez: ERRRRROR Xlen: %i Ylen %i ' %(xLen,yLen)
        return 'ERRROR'
    #Apply thresholding

    X = copy.deepcopy(X)
    Y = copy.deepcopy(Y)
    acceptedNums = 0
    for i in range(xLen):
        x = float(X[i])
        y = float(Y[i])
        if testFunction(x) and testFunction(y):
            acceptedNums +=1
            X[i] = x
            Y[i] = y
        else:
            X[i] = ''
            Y[i] = ''
    print 'ThresholdCorrel accepted %i of %i numbers' %(acceptedNums, xLen)
    #Calculate XY
    sumXY = 0
    sumX = 0
    sumXsqrd = 0
    sumY = 0
    sumYsqrd = 0
    N = xLen
    for i in range(xLen):
        x = X[i]
        y = Y[i]
        if not(x):
            N -=1
            continue
        if not(y):
            N -=1
            continue
        x = float(x)
        y = float(y)
        ###print str(x)+ ' and ' + str(y)
        sumX += x
        sumXsqrd += pow(x,2)
        sumY += y
        sumYsqrd += pow(y,2)
        sumXY += x * y
    rTop = sumXY - (sumX*sumY/N)
    rBottom = math.sqrt((sumXsqrd-pow(sumX,2)/N)*(sumYsqrd-pow(sumY,2)/N))
    if rBottom == 0: return 0
    correl = rTop/rBottom
    return correl

def ThresholdBinaryCorrel(X,Y, testFunction = None):
    """Correl(X,Y)| please fillininfohere|need two paired lists"""
    if not(testFunction):
        exec 'testFunction = lambda x:1'
    ##CheckInput
    xLen = len(X)
    yLen =len(Y)
    OK = None
    if yLen == xLen:
        OK = 1
    if yLen > 0 and OK:
        OK = 1
    if not(OK):
        print 'Correl Sez: ERRRRROR Xlen: %i Ylen %i ' %(xLen,yLen)
        return 'ERRROR'
    #Apply thresholding

    X = copy.deepcopy(X)
    Y = copy.deepcopy(Y)
    acceptedNums = 0
    for i in range(xLen):
        Xval = X[i]
        Yval = Y[i]
        if (Xval =='') or (Yval ==''): continue
        x = float(Xval)
        y = float(Yval)
        if testFunction(x) and testFunction(y):
            acceptedNums +=1
            X[i] = x
            Y[i] = y
        else:
            X[i] = ''
            Y[i] = ''
    #print 'ThresholdBinaryCorrel accepted %i of %i numbers' %(acceptedNums, xLen)

    return float(acceptedNums)
        
def ThresholdBinaryScore(X,Y, threshold = None, greaterThan = 1):
    """Correl(X,Y)| please fillininfohere|need two paired lists"""
    
    ##CheckInput
    xLen = len(X)
    yLen =len(Y)
    OK = None
    if yLen == xLen:
        OK = 1
    if yLen > 0 and OK:
        OK = 1
    if not(OK):
        print 'Correl Sez: ERRRRROR Xlen: %i Ylen %i ' %(xLen,yLen)
        return 'ERRROR'
    #Apply thresholding
    if threshold:
        X = copy.deepcopy(X)
        Y = copy.deepcopy(Y)
        for i in range(xLen):
            x = float(X[i])
            y = float(Y[i])
            if greaterThan:
                if x > threshold:
                    X[i] = x
                else:
                    X[i] = ''
                if y > threshold:
                    Y[i] = y
                else:
                    Y[i] = ''
            else:
                if x < threshold:
                    X[i] = x
                else:
                    X[i] = ''
                if y < threshold:
                    Y[i] = y
                else:
                    Y[i] = ''
    score = 0.0
    for i in range(xLen):
        if X[i] and Y[i]:
            xy = X[i]/Y[i]
            yx = Y[i]/X[i]
            ratio = min(yx, xy) ##This covers the interaction of two imperfect probes
                                ## The most you an add at each probe is one probe
            score += ratio
    
    return score


def AnnotateDirectRepeats(IFNorLst,UD,sortCol = 'seq',OFN = None):
    """AnnotateDirectRepeats
    |(IFNorLst,sortCol = 'seq',OFN = None)
    |Output: filename or outDB
    |writes clusters of duplications to each probe or object name returns data in same
    order as input
    ||"""
    if str(type(IFNorLst)).find('list') <>-1:
        inDB = IFNorLst
    else:
        IFN_T = IFNorLst
        inDB = InArray(IFN_T,UD,separator = '\t',Verbose = 1)
    ##Prep for unsort
    for j in range(len(inDB)):
        inDB[j] = [j] + inDB[j]
    ##decide sortCol
    if str(type(sortCol)).find('str'):
        labelLine = inDB.pop(0)
        sortCol = labelLine.index(sortCol)
    else:
        sortCol = int(sortCol) + 1
    inDB_BU = ListBreakup(inDB,sortCol)
    labelLine.extend(['clust','#members'])
    outDB = []
    outDB.append(labelLine)
    clust = 1
    for group in inDB_BU:
        groupLen = len(group)
        if groupLen == 1:
            groupLen = 'uniq'
        for obj in group:
            label = obj[2]
            obj.extend([clust,groupLen])
            outDB.append(obj)
        clust +=1
    outDB.sort()
    for j in range(len(outDB)):
        outDB[j] = outDB[j][1:]
    if OFN:
        clOFN_T = 'directRepeatLabelled_'+ IFN_T
        clOFN_T = OutArray(outDB,clOFN_T,UD)
        return clOFN_T
    else:
        return outDB

def TransverseArray(inDB, NothingChar = '', Verbose = 1):
    """TransverseArray
    |(inDB, NothingChar = '')
    |Output: outDB
    |Transverses an ixj array where i <> j
    ||"""
    inDBLen = len(inDB)
    lengths = []
    for obj in inDB:
        lengths.append(len(obj))
    maxLen = max(lengths)
    
    outDB = [[NothingChar for h in range(inDBLen) ] for i in range(maxLen)]
    nextc = 1000
    for c in range(inDBLen):
        if c > nextc: #Verbose:
            nextc = c + 1000
            print 'TransverseArray: is doing column %i of %i.' %(c,inDBLen)
        for r in range(maxLen):
            procList = inDB[c]
            LenProcList = len(procList)
            if r >= LenProcList: break
            outDB[r][c] = procList[r]
    return outDB



def TransverseDataSet(inDB):
    """TransverseDataSet(inDB)| Only takes symmtrical matrices|"""
    ##converts a list of lists into Excel-exportable text
    cols = len(inDB)
    #All this builds up the empty set
    row = []
    for j in range(cols):
        row.append('')
    rowLenList = []
    for r in inDB:
        rowLenList.append(len(r))
    maxRow = max(rowLenList)
    print 'TransverseDataSet: is transversing an array with %i rows and %i cols' %(maxRow,cols)
    outDB = []
    for obj in range(maxRow):
        outDB.append(row[:])
    #now we are ready to populate the dataset
    for i in range(cols):
        obj = inDB[i]
        for j in range(len(obj)):
            outDB[j][i] = obj[j]
    return outDB


def IndexFileLineNumbers(IFN):
    inFP = IFN.FP(RW='r')
    FilePosIndexLst = []
    line = ' '
    lineNum = 0
    nextL = 0
    while line:
        FilePos = inFP.tell()
        line = inFP.readline()
        lineNum +=1
        FilePosIndexLst.append(FilePos)
        if lineNum > nextL:
            nextL = lineNum + 100
            print 'indexed %i lines last line was: %s' %(lineNum, line[:15])
    inFP.close()
    return FilePosIndexLst


def TransverseArray_ME(IFN, RowsToRemoveFromTop = 0):
    BlockToProcessSize = 100
    
    DEBUG = 0
    
    inFP = IFN.FP(RW='r')
    FilePosIndexLst = IndexFileLineNumbers(IFN)
    RowSize = len(FilePosIndexLst)
    
    fpos = FilePosIndexLst[0]
    inFP.seek(fpos)
    lineLst = inFP.readline().replace('\n','').split('\t')
    ColSize = len(lineLst)
    fpos = FilePosIndexLst[-2]
    inFP.seek(fpos)
    lineLst = inFP.readline().replace('\n','').split('\t')
    EndColSize = len(lineLst)
    print 'FirstRowLen %i LastRowLen %i' %(ColSize, EndColSize)
    #ret = raw_input('Continue ? (any key):')

    ofn = 'Transversed_' + IFN.fn
    OFN = FA.path([IFN.d, ofn],fileType = 'tab')
    outFP = OFN.FP('w')
    
    
    IntervalsToProcess = range(0,ColSize,BlockToProcessSize)
    i = 0
    MatrixSize = 0
    for interval in IntervalsToProcess:
        procMatrix = []
        IntervalStart = interval
        IntervalEnd = interval + BlockToProcessSize
##        if IntervalEnd >= MatrixSize:
##                IntervalEnd = MatrixSize
        print 'Processing %i of %i intervals cols %i to %i' %(i, len(IntervalsToProcess),IntervalStart,IntervalEnd)
        l = 0
        nextl = 1000
        for linePos in range(RowSize):
            fpos = FilePosIndexLst[linePos]
            inFP.seek(fpos)
            lineLst = inFP.readline().replace('\n','').split('\t')
            if len(lineLst) ==1: continue ##010908_this avoids padding the rows with garbageThis may cause problems later
            appendLst = lineLst[IntervalStart:IntervalEnd]
            procMatrix.append(appendLst)
            l +=1
            if l > nextl:
                showString = str(appendLst)[:25]
                print 'LinePos %i fpos: %s   %s' %(linePos,str(fpos),showString)
                nextl = l+ 100
        TV_procMatrix = TransverseArray(procMatrix, Verbose = 1)

        ###########Just this once!!!
        if RowsToRemoveFromTop:
            for p in range(RowsToRemoveFromTop):
                yo = TV_procMatrix.pop(0)
                print 'Popping ' + str(yo)[:20]
            RowsToRemoveFromTop = None
            
        yo = OutArray(TV_procMatrix, outFP, 'll')
        TV_procMatrix = []
        i +=1
    outFP.close()
    inFP.close()
    return OFN



def Histo(inDB,index=None,numBins = 10, minNum = None, maxNum = None):
    """Histo(inDB,index=None,multiplier = 1)| Exclusive property of IM|"""
    inDATA = []
    Indexfunct = DataAccessFunctionProcessor(index)
##    if index:
##        for obj in inDB:
##            if str(obj[index]) <> 'None':
##                inDATA.append(float(obj[index]))
##    else:
    for obj in inDB:
        inDATA.append(float(Indexfunct(obj)))
    
    if not(maxNum): maxNum  = max(inDATA)
    if not(minNum): minNum = min(inDATA)
    ##Dont Normalize !!!!! I do it but I dont use it
    Normalize = None
    if Normalize:
        NormalinDATA = []
        for n in range(len(inDATA)):
            num = inDATA[n]
            num = ((num-minNum)/(maxNum-minNum))
            NormalinDATA.append(num)
        inDATA = NormalinDATA
    
    intervals = (maxNum-minNum)/numBins
    ##print str(intervals)
    histo = []
    binList = []
    val = minNum
    for k in range(numBins+2):
        histo.append([val,0])
        binList.append(val)
        val = val + intervals
    inDATA.sort(lambda x,y: cmp(float(x),float(y)))
    k = 0
    i = 0
    while i < len(inDATA) and k < numBins+2:
        above = float(binList[k])
        val = inDATA[i]
        ###print 'i %i %f between %f and %f' %(i,val,below,above)
        if float(val) <= above:
            histo[k][1] +=1
            i+=1
        else:
            k +=1
                
    return histo

def PrintHistoBars(histo, breakup = None):
    maxHisto = MaxArr(histo,1)
    minHisto = MinArr(histo,1)
    printObj = 'val\t\telems\n'
    for obj in histo:
        numElems = obj[1]
        numBars = int(numElems/maxHisto * 40) + 1
        binNum = str(obj[0])[:5]
        prio = binNum + '\t\t' +str(numElems) + '\t' + ''.join(['|' for i in range(numBars)]) + '\n'
        printObj = printObj + prio
    print printObj
    if breakup:
        tot = SumArr(histo,1)
        stoppoint = int(tot * breakup)
        s = 0
        i = 0
        while s <stoppoint:
            i +=1
            s = s + histo[i][1]
            
        print 'bin %s represents %s point s = %s total = %s' %(str(histo[i][0]), str(breakup),str(s),str(tot))
        return histo[i][0]
    return printObj


def stdev(inList):
    """stdev(inList)| Exclusive property of IM|"""
    avg = Average(inList)
    i = 0
    sumSqrDiff = 0.0
    for x in inList:
        if x == None: continue ##Added 050709
        x = float(x)
        i += 1
        sumSqrDiff += math.pow(avg - x,2)
    return math.sqrt(sumSqrDiff/i)


def stdevNonBiased(inList):
    
    """stdevNonBiased| Excel used this it is diviede by n-1|"""
    avg = Average(inList)
    i = 0
    sumSqrDiff = 0.0
    for x in inList:
        x = float(x)
        i += 1
        sumSqrDiff += math.pow(avg - x,2)
    i = max(1, i-1)
    return math.sqrt(sumSqrDiff/i)

def IndexedZscore(inList,startIndexOrLabel = 1, DataIndex = None):
    """Zscore(inList,startIndexOrLabel = 1)| Exclusive property of IM|"""
    ##Adds the z-score to your data This function is stupid!
    if DataIndex==None:
        zscore = Zscore(inList,startIndexOrLabel)
        return zscore

    zlist = []
    for i in range(startIndexOrLabel,len(inList),1):
        zlist.append(float(inList[i][DataIndex]))
    zscore = Zscore(zlist,startIndexOrLabel)
    for i in range(startIndexOrLabel,len(inList),1):
        zlist_Index = i - startIndexOrLabel
        inList[i].append(zscore[zlist_Index])
    
    return inList



def Zscore(inList,startIndexOrLabel = 1):
    """Zscore(inList,startIndexOrLabel = 1)| Exclusive property of IM|"""
    #This function is also stupid it trie sto do too much it attempt to append the z
    if 'str' in str(type(startIndexOrLabel)):
        label = startIndexOrLabel
        startIndex = 0
    elif startIndexOrLabel == 0:
        startIndex = 0
        outList = []
    else:
        label = 'Z-'+ inList[0]
        startIndex = startIndexOrLabel
        outList = [label]
    avg = Average(inList[startIndex:])
    standdev = stdev(inList[startIndex:])
    
    for x in inList[startIndex:]:
        if x == None:##Added 050709
            outList.append(None)
            continue
        x = float(x)
        if standdev == 0:
            val = 0
        else:
            val = (x-avg)/standdev
        outList.append(val)
    return outList




def justZscore(inList,indexer = None):
    """Zscore(inList,startIndexOrLabel = 0)| Exclusive property of IM|"""
    DataAccF = DataAccessFunctionProcessor(indexer)
    
    outList = []
    
    avg, stdev = AvgStdevArr(inList,function = DataAccF)
    
    for x in inList:
        x = float(DataAccF(x))
        outList.append((x-avg)/stdev)
    return outList



def StreamReader(Size,infile,Demarkator = '//',Verbose = 0):
    """StreamReader
    |(Size,infile,Demarkator,Verbose = 0)
    |Output: big string
    |This returns Size (bytes) of data from an open file already named infile
    in calling function or in module namespace It also reads a little extra
    to get to the next demarkator
    ||"""
    
    p = ' '
    OutStr = infile.read(Size)
    p = infile.readline()
    OutStr = OutStr + p
    if OutStr[-5:].find(Demarkator)<>-1:
        print 'easy return'
        return OutStr     
    else:
        atEnd = None
        while not(atEnd) and p <> '':
            p = infile.readline()
            OutStr = OutStr + p
            if p.find(Demarkator)<>-1:
                atEnd = 1
    if Verbose:
        readlen = len(OutStr)
        print ' StreamReader just read chars length of: ' + str(readlen)
    return OutStr

def FindNth(Instr, char,n):
        """FindNth
        |(Instr, char,n)
        |Output: int
        |Finds the Nth occurrence of a string in another string.||"""
        x = 0
        b = 0
        for y in range(n):
            x = Instr.find(char, x + 1)
        return x
    
def DumpDuplicates(inLst, cmpFunction = None, Verbose = 0):
    
    inLst_BU = FastListBreakup(inLst, cmpFunction, Verbose = Verbose)
    outLst = []
    for obj in inLst_BU:
        outLst.append(obj[0])
    return outLst

def FastListBreakup(inList,indexer = None,Verbose = None):
    """FastListBreakup
    |(inList,index)
    |Output: LIst of Lists
    |Sorts at the begining, gains a lot of speed that way!!||"""
    
    if not(inList):
        return []
    Indexfunct = DataAccessFunctionProcessor(indexer)
    
##    if indexer ==None:
##        indexer = lambda x:x
##    if 'int' in str(type(indexer)):
##        if len(inList[0]) <= indexer:
##            print 'Size of index is Wrong!! %s, it should be lessthanequalTo %i' %(str(indexer), len(inList[0]))
##            return []
##    if 'function' in str(type(indexer)):
##        Indexfunct = indexer
##    elif 'int' in str(type(indexer)):
##        Indexfunct = lambda x :x[indexer]
##    else:
##        print 'Indexer: %s is wrong please correct this' %(str(indexer))
##        return []
    inList.sort(lambda x,y:cmp(Indexfunct(x),Indexfunct(y)))
    previous = Indexfunct(inList[0])
    Out_BU = []
    miniList = []
    i = 0
    clusterSize = 0
    ObjectsToProcess = len(inList)
    objectsProcessed = 0
    for obj in inList:
        objectsProcessed +=1
        if previous == Indexfunct(obj):
            miniList.append(obj)
        else:
            if len(miniList) <> 0:
                Out_BU.append(miniList)
                clutserSize = len(miniList)
            miniList = [obj]
            previous = Indexfunct(obj)
            i +=1
            if Verbose:
                print 'Cluster %i \t has %i members. \t Processed %i of %i' %(i-1,clutserSize,objectsProcessed,ObjectsToProcess)
    if len(miniList) >0:
        Out_BU.append(miniList)
        if Verbose:
            clusterSize = len(miniList)
            print 'Cluster %i \t has %i members' %(i,clusterSize)
    return Out_BU

def FastDictIndex(inList,Sortindex,TakeIndex, returnDict = 1,verbose = None):
    """FastListIndex
    |(inList,index)
    |Output: LIst of Lists
    |Sorts at the begining, gains a lot of speed that way!!||"""
    Out_BU = FastListBreakup(inList,Sortindex,Verbose = verbose)
    if returnDict:
        outObj = {}
    else:
        outObj = []
    for obj in Out_BU:
        item = obj[0][Sortindex]
        BU_Obj = [item,[]]
        for buObj in obj:
            ##CAUTION CAUTION CAUTION This clears redundancies
            if not(buObj[TakeIndex] in BU_Obj[1]):
                BU_Obj[1].append(buObj[TakeIndex])
        if returnDict:
            outObj[BU_Obj[0]] = BU_Obj[1]
        else:
            outObj.append(BU_Obj)
    return outObj


def LambdaListBreakup(inList,inFunct,Verbose = None):
    """LambdaFastListBreakup
    |(inList,inFunct,Verbose = None)
    |Output: LIst of Lists
    |Sorts at the begining, gains a lot of speed that way!!||"""

    inList.sort(lambda x,y:cmp(inFunct(x),inFunct(y)))
    previous = inFunct(inList[0])
    Out_BU = []
    miniList = []
    i = 0
    clusterSize = 0
    ObjectsToProcess = len(inList)
    objectsProcessed = 0
    for obj in inList:
        objectsProcessed +=1
        if previous == inFunct(obj):
            miniList.append(obj)
        else:
            if len(miniList) <> 0:
                Out_BU.append(miniList)
                clutserSize = len(miniList)
            miniList = [obj]
            previous = inFunct(obj)
            i +=1
            if Verbose:
                print 'Cluster %i \t has %i members. \t Processed %i of %i' %(i-1,clutserSize,objectsProcessed,ObjectsToProcess)
    if len(miniList) >0:
        Out_BU.append(miniList)
        if Verbose:
            clusterSize = len(miniList)
            print 'Cluster %i \t has %i members' %(i,clusterSize)
    return Out_BU
def ListBreakup(inList,index):
    """ListBreakup
    |(inList,index)
    |Output: LIst of Lists
    |Allows the user to breakup a list into groups having the same values for
    one element in the list.||"""
    group = []
    outList = []
    for miniList in inList:
        if len(miniList) < index:
            print 'listBreakup: BAD index sent' + str(index) + str(miniList)
            return None
        if not(miniList[index] in group):
            group.append(miniList[index])
    for element in group:
        subList = []
        for miniList in inList:
            if element == miniList[index]:
                subList.append(miniList)
        outList.append(subList)
    return outList


def nTallyHo(inDB, indexer = None, Verbose = None):
    """TallyHo
    |(inDB, refIndex, averageIndex = None, Verbose = None)
    |Output: List of lists
    |||"""
    outDB = []
    if not(inDB):
        return []
    Indexfunct = DataAccessFunctionProcessor(indexer)
    processedinDB = FastListBreakup(inDB,Indexfunct)
    
    
    for node in processedinDB:
        outDB.append([Indexfunct(node[0]),len(node)])
    if Verbose:
        print 'tallyHo returns %i groups from %i members' %(len(outDB),len(inDB))
        outDB.sort(lambda x,y:cmp(x[-1],y[-1]))
        outDB.reverse()
        for obj in outDB:
            print 'Cluster:\t %s \t\t\t members %i' %(str(obj[0]),obj[1])
    return outDB


def TallyHo(inDB, refIndex = None, averageIndex = None, Verbose = None):
    """TallyHo
    |(inDB, refIndex, averageIndex = None, Verbose = None)
    |Output: List of lists
    |||"""
    outDB = []

    processedinDB = FastListBreakup(inDB,refIndex)
    
    if averageIndex == None:
        for node in processedinDB:
            if 'function' in str(type(refIndex)):
                outDB.append([refIndex(node[0]),len(node)])
            else:
                outDB.append([node[0][refIndex],len(node)])
    else:
        avgLst = []
        for node in processedinDB:
            label = node[0][refIndex]
            avgLst = []
            for piece in node:
                avgLst.append(piece[averageIndex])
            avg = Average(avgLst,returned = True)
            outDB.append([label,avg[0],avg[1]])
    if Verbose:
        print 'TallyHo returns %i groups from %i members' %(len(outDB),len(inDB))
        outDB.sort(lambda x,y:cmp(x[-1],y[-1]))
        outDB.reverse()
        for obj in outDB:
            print 'Cluster:\t %s \t\t\t members %i' %(str(obj[0]),obj[1])
    return outDB



def oldTallyHo(inDB, refIndex = None, averageIndex = None, Verbose = None):
    """TallyHo
    |(inDB, refIndex, averageIndex = None, Verbose = None)
    |Output: List of lists
    |||"""
    outDB = []
    if refIndex == None:
        for j in range(len(inDB)):
            inDB[j] = [inDB[j]]
        refIndex = 0
    processedinDB = ListBreakup(inDB,refIndex)
    
    if averageIndex == None:
        for node in processedinDB:
            outDB.append([node[0][refIndex],len(node)])
    else:
        avgLst = []
        for node in processedinDB:
            label = node[0][refIndex]
            avgLst = []
            for piece in node:
                avgLst.append(piece[averageIndex])
            avg = Average(avgLst,returned = True)
            outDB.append([label,avg[0],avg[1]])
    if Verbose:
        print 'tallyHo returns %i groups from %i members' %(len(outDB),len(inDB))
        outDB.sort(lambda x,y:cmp(x[-1],y[-1]))
        outDB.reverse()
        for obj in outDB:
            print 'Cluster:\t %s \t\t\t members %i' %(str(obj[0]),obj[1])
    return outDB

    
def Countif(inArr, ifFunction = None, colIndex = None, rowIndex = None):
    arr = []
    if not(ifFunction):
        ifFunction = lambda x: (1 if x > 0 else 0)
    if colIndex:
        for obj in inArr:
            arr.append(obj[colIndex])
    else:
        arr = copy.deepcopy(inArr)
    if rowIndex:
        arr = arr[rowIndex:]
    for i in range(len(arr)):
        arr[i] = ifFunction(arr[i])
    return SumArr(arr)

def SumArr(Arr,index=None):
    """SumArr
    |(Arr,index=None)
    |Arrsum
    |can find sum of list of nums or of indexed numbers in list of lists
    | |"""
    Indexfunct = DataAccessFunctionProcessor(index)
    Arrsum = 0
    for obj in Arr:
        val = Indexfunct(obj)
        if val <> None:
            Arrsum += float(val)
    return Arrsum

def BinSumArr(Arr,index=None):
    """BinSumArr
    |(Arr,index=None)
    |Arrsum
    |can find sum of list of nums or of indexed numbers in list of lists
    | |"""
    Arrsum = 0
    for obj in Arr:
        if index <> None:
            if obj[index]:
                Arrsum += 1
        else:
            if obj:
                Arrsum += 1
    return Arrsum

def Average(inLst, returned = None):
    """Average
    |(inLst, returned = False)
    |Output: List or number
    |Average of list OR [Avg, # of elements in List] (returned = T / F)||"""
    length = len(inLst)
    i = 0
    Sum = 0
    if length>1:
        for elem in inLst:
            if elem == None: continue ##I sure hope this doesn't fuck anything up 051107
            if 'str' in str(type(elem)):
                if elem.isdigit():
                    elem = float(elem)
                else:
                    return None
            i += 1
            Sum += float(elem)
        if i ==0:
            avg = 'nd'
        else:
            avg = Sum/i
        
    elif length == 1:
        avg = inLst[0]
    else:
        avg = 'nd'
    if returned:
        return [avg,length]
    else:
        return avg
    

def AvgArr(Arr,index = None):
    """AvgArr
    |(Arr,index)
    |Output: number
    |Returns the average of an array at a particular index.||"""
    i = 0.0
    ssum = 0.0
    for obj in Arr:
        if obj =='': continue
        if index == None:
            ssum = ssum + float(obj)
        else:
            ssum = ssum + float(obj[index])
        i += 1
    if i >0:
        avg = ssum/i
    else:
        avg = None
    return avg

##def isOdd(val):
##    halfValFloat = float(val)/2.0
##    halfValInt = int(val)/2
##    if int(abs(halfValFloat - float(halfValInt))*10) ==5:
##        return 1
##    else:
##        return None



def Median(signalLst, indexer = None):
    IndexFunct = DataAccessFunctionProcessor(indexer)
    tLst = []
    for val in signalLst:
        tLst.append(float(IndexFunct(val)))
    v = len(tLst)
    tLst.sort()
    if isOdd(len(tLst)):
        median = tLst[int(v/2)+1]
    else:
        median = (tLst[v/2] + tLst[v/2+1])/2
    return median
                

def MedianAbsoluteDeviation(signalLst, indexer = None):
    Indexfunct = DataAccessFunctionProcessor(indexer)          
    median = Median(signalLst, indexer = indexer)
    dLst = []
    for v in signalLst:
        dLst.append(abs(float(Indexfunct(v))-median))
    return Median(dLst)

def MADZscore(inList,indexer = None):
    #Median Absolute Ddeviation
    #To z-score non-Normal data
    """MADZscore(inList,startIndexOrLabel = 1)| Exclusive property of IM|"""
    Indexfunct = DataAccessFunctionProcessor(indexer)
    
    median = Median(inList, indexer = indexer)
    MAD = MedianAbsoluteDeviation(inList)
    outLst = []
    for x in inList:
        x = float(Indexfunct(x))
        outLst.append((x-median)/MAD)
    return outLst

def NormalScore(signalLst):
    #ratio should be ~ 1.4826 for Normal dist data
    stdev = LM.Stdev(signalLst)
    MAD = MedianAbsoluteDeviation(signalLst)
    ratio = stdev/MAD
    return ratio/1.4826

def isEven(n):
   """Return true if n is even."""
   return n%2==0

def isOdd(n):
   """Return true if n is odd."""   
   return not isEven(n)
    
def MedianArr(Arr,index = None):
    """MedianArr
    |(Arr,index)
    |Output: number
    |Returns the average of an array at a particular index.||"""
    Indexfunct = DataAccessFunctionProcessor(index)
    valLst = []
    for obj in Arr:
        valLst.append(float(Indexfunct(obj)))
    numVals = len(valLst)
    if numVals == 0:
        return None
    valLst.sort()
    if isOdd(numVals):
        #return Central number
        return valLst[numVals/2]
    else:
        #isEven average the two mid vals
        halfPoint = numVals/2
        return Average([valLst[halfPoint-1], valLst[halfPoint]])

def AvgStdevNArr_outlierremover(Arr,indexer = None):
    ##Can dump up to 10% of samples to decrese stdev scatter
    """INCOMPLETE
    |(Arr,index)
    |Output: number
    |Returns the average of an array at a particular index.||"""
    Indexfunct = DataAccessFunctionProcessor(indexer)
##    if not(function) and not(index):
##        function = lambda x: x
##    elif not(index == None):
##        exec 'function = lambda x:x[%i]' %index
    Array = [Indexfunct(x) for x in Arr]
    array = []
    while Array:
        a = Array.pop(0)
        if str(a) == 'None': continue
        if str(a).isalpha():continue
        array.append(a)
    array.sort()
    n = len(array)
    outlierD =(1+int(n/10))
    avg, stdev = justAvgStdevArr(array)
    testLst = [[n, avg, stdev]]
    testRangeSize = n -outlierD
    for i in range(1, outlierD, 1):
        runlst =array[i:testRangeSize+i]
        n = len(runlst)
        avg, stdev = justAvgStdevArr(runlst)
        n = testRangeSize
        print '%i from %i to %i avg: %.2f stdev is %.2f' %(i, i, testRangeSize+i, avg, stdev)
        testLst.append([n, avg, stdev])
    testLst.sort()
    w = testLst[0]
    return w

    
def justAvgStdevArr(array):
    avg = AvgArr(array)
    XminusXmean = 0.0
    n = 0
    for n in range(len(array)):
        num = array[n]
        x = float(num)
        XminusXmean += math.pow(x-avg,2)
    if n ==0:
        return avg, 0.0
    stdev = XminusXmean/(n+1)
    stdev = math.sqrt(stdev)
    return avg, stdev

def AvgStdevArr(Arr,index = None, function = None, startIndex = 0):
    """INCOMPLETE
    |(Arr,index)
    |Output: number
    |Returns the average of an array at a particular index.||"""
    Indexfunct = DataAccessFunctionProcessor(function)
##    if not(function) and not(index):
##        function = lambda x: x
##    elif not(index == None):
##        exec 'function = lambda x:x[%i]' %index
    Array = [Indexfunct(x) for x in Arr[startIndex:]]
    array = []
    while Array:
        a = Array.pop(0)
        if str(a) == 'None': continue
        if str(a).isalpha():continue
        array.append(a)
            
    index = None
    #print str(array[:10])
    avg = AvgArr(array)
    XminusXmean = 0.0
    n = 0
    for n in range(len(array)):
        num = array[n]
        x = float(num)
        XminusXmean += math.pow(x-avg,2)
    if n ==0:
        return avg, 0.0
    stdev = XminusXmean/(n+1)#n is zero-based therfore we add to it in order to get N
    #N-1 (or juts plain n) is correct for a sample stdev.
    stdev = math.sqrt(stdev)
    return avg, stdev

def AvgStdevNArr_NminusOne(Arr,index = None, function = None, startIndex = 0):
    """INCOMPLETE
    |(Arr,index)
    |Output: number
    |Returns the average of an array at a particular index.||"""
    Indexfunct = DataAccessFunctionProcessor(function)
##    if not(function) and not(index):
##        function = lambda x: x
##    elif not(index == None):
##        exec 'function = lambda x:x[%i]' %index
    Array = [Indexfunct(x) for x in Arr[startIndex:]]
    array = []
    while Array:
        a = Array.pop(0)
        if str(a) == 'None': continue
        if str(a).isalpha():continue
        array.append(a)
            
    index = None
    avg = AvgArr(array, index)
    XminusXmean = 0.0
    n = 0
    for n in range(len(array)):
        num = array[n]
        x = float(num)
        XminusXmean += math.pow(x-avg,2)
    if n ==0:
        return avg, 0.0, 0
    stdev = XminusXmean/n##The n-1 is to the left The original is on the right sinze we are zero-based (n+1)
    stdev = math.sqrt(stdev)
    return avg, stdev, n

##Deprecated 060708 this is a mess
##def MaxArr(Arr,index = None,retIndex = None):
##    """MaxArr
##    |(Arr,index,retIndex = None)
##    |return maxNum
##    |returns max number of index of array
##    | |"""
##    if not(Arr):
##        return 0
##    if index:
##        maxNum = -999999999.0
##        i = 0
##        for obj in Arr:
##            num = float(obj[index])##Here you will ned a exception catcher for int or float casting
##            if  num > maxNum:
##                maxNum = num
##                maxIndex = i
##            i += 1
##    else:
##        maxNum = -999999999.0
##        for obj in Arr:
##            if not(index):
##                num = float(obj)
##            else:
##                num = float(obj[index])
##            if  num > maxNum:
##                maxNum = num
##                maxIndex = i
##    if retIndex:
##        return maxIndex
##    else:
##        return maxNum
##def MinArr(Arr,index = None, retIndex = None):
##    """MinArr
##    |(Arr,index,retIndex = None)
##    |return minNum
##    |returns min number af array in index
##    | |"""
##    if not(index):
##        minNum = 999999999
##        for obj in Arr:
##            if not(index):
##                num = int(obj)
##            else:
##                num = int(obj[index])
##            if  num < minNum:
##                minNum = num
##    else:
##        minNum = 999999999
##        i = 0
##        for obj in Arr:
##            num = int(obj[index])
##            if  num < minNum:
##                if retIndex:
##                    minNum = i
##                else:
##                    minNum = num
##            i += 1
##    return minNum

def MaxArr(Arr,index = None,retIndex = None):
    """MaxArr
    |(Arr,index,retIndex = None)
    |return maxNum
    |returns max number of index of array
    | |"""
    if not(Arr):
        return None
    maxNum = -9e100
    if index <> None:
        for i in range(len(Arr)):
            num = float(Arr[i][index])##Here you will ned a exception catcher for int or float casting
            if  num > maxNum:
                maxNum = num
                maxIndex = i
    else:
        for i in range(len(Arr)):
            num = float(Arr[i])
            if  num > maxNum:
                maxNum = num
                maxIndex = i
    if retIndex:
        return maxIndex
    else:
        return maxNum
    
def MinArr(Arr,index = None, retIndex = None):
    """MinArr
    |(Arr,index,retIndex = None)
    |return minNum
    |returns min number af array in index
    | |"""
    if not(Arr):
        return None
    minNum = 9e100
    if index <> None:
        for i in range(len(Arr)):
            num = float(Arr[i][index])##Here you will ned a exception catcher for int or float casting
            if  num < minNum:
                minNum = num
                minIndex = i
    else:
        for i in range(len(Arr)):
            num = float(Arr[i])
            if  num < minNum:
                minNum = num
                minIndex = i
    if retIndex:
        return minIndex
    else:
        return minNum


def MaMiRec(inDB,Index=None,MaMi = 'Min'):
    """MaMiRec
    |(inDB,index,MaMi = 'Min')
    |inDB[mami]
    |Max or Min of an array object on specified index it calls Max or Min only Min is specified and is default
    | |"""
    if MaMi == 'Min':
        mami = MinArr(inDB,index=Index,retIndex = 1)
    else:
        mami = MaxArr(inDB,index=Index,retIndex = 1)
    return inDB[mami]

def pInArray(path):
    return InArray(path.fn , path.d, separator = '\t',Verbose = 0)

def InArray(IFN , userDir, separator = '\t',Verbose = 0):
    """InArray
    |(infilename , userDir, separator = '\t',Verbose = 0)
    |Output: List of FASTAs
    |This is for inputting tab-delimited arrays||"""
    if str(IFN).find('FASTA.path instance') <> -1:
        infilename = IFN.fn
        userDir = IFN.d
    else:
        infilename = IFN
    outFasta = []
    infile = open(userDir + infilename, 'r')
    s = infile.read()
    infile.close()
    Li = s.split('\n')
    for line in Li:
        if (line <> 0 and line <> ''):
            fastRecord = line.split(separator)
            for o in range(len(fastRecord)):
                if fastRecord[o] == 'None':
                    fastRecord[o] = None
            outFasta.append(fastRecord)
    if Verbose == 1:
        print 'InArray: pulled in ' + str(len(outFasta)) + ' Objects'
    return outFasta


def OutArray(dataArray, outFileNameOrFilePointer, userDir, DEBUG = None,outType = '\t', Verbose = 0):
    ##Modified 10-08-08
    """OutArray
    |(dataArray, outFileNameOrFilePointer, userDir, outType = '\t', Verbose = 0)
    |Output: File
    |Exports an array of data to file using either tabs between pieces
    This program now allows streamwriting if you send it a filepointer instead of a filename
    ||"""
    if str(type(outFileNameOrFilePointer)).find('file')<>-1:
        outFP = outFileNameOrFilePointer
        outfilename = outFP.tell()
        CloseFile = None
    else:
        outfilename = outFileNameOrFilePointer
        if outType == ',':
            if outfilename.find('txt')<>-1:
                outfilename = outfilename.replace('.txt','.csv')
            else:
                outfilename = outfilename + '.csv'
        outFP = open(userDir + outfilename,'w')
        CloseFile = 1
    for line in dataArray:
        if str(type(line)).find('list')<>-1:
            if not(DEBUG):
                wLine = ''
                for element in line:
                    wLine = wLine + str(element) + outType
                outFP.writelines(wLine[:-1])
            else:
                for element in line:
                    outFP.write(str(element) + outType)
        else:
            outFP.write(str(line))
        outFP.writelines('\n')
    if CloseFile:
        outFP.close()
    if Verbose == 1:
        print 'YoYoYo' + str(outFP.tell())
        print 'OutArray: Wrote ' + str(len(dataArray)) + ' Objects to file: ' + str(outFP)
        
    return outfilename


def ExcelOutArray(inDB, ofn, UD, breakupSize = 65000, LL =-1):
    if LL == -1:
        LL = inDB.pop(0)
    if LL<>-1:
        outDBgenerator = lambda :[LL]
    else:
        outDBgenerator = lambda :[]
    outDB = []
    pt = 0
    outObj = US.FilePathLstObj()
    outObj.type = 'tab'
    for s in range(0,len(inDB), breakupSize):
        pt +=1
        e = s + breakupSize
        outDB = outDBgenerator() + inDB[s:e]
        ptofn = 'pt' + str(pt) + '_' + ofn
        OutArray(outDB, ptofn, UD, Verbose = 1)
        outObj.add([UD + ptofn])
    return outObj



def TopPart(inDB, part, refIndex, Verbose = 1):
    """TopPart
    |(inDB,part,refIndex,Verbose = 1)
    |Output: 
    |This function pulls the top part (real fractional number <1, >0
    from list of normalized objects||"""
    maxVal = -100000
    minVal = 999999
    normalIndex = []
    i = 0
    for obj in inDB:
        val = obj[refIndex]
        normalIndex.append([i,val])
        i += 1
        if val >maxVal:
            maxVal = val
        if val <minVal:
            minVal = val
    normalIndex.sort
    normalIndex.reverse()
    inDB.sort(lambda x,y,refIndex: cmp(x[refIndex],y[refIndex]))
    inDB.reverse()
    for i in range(len(normalIndex)):
        normNum = (float(normalIndex[i][1])-minVal)/(maxVal-minVal)
        if normNum < part:
            return inDB[:i]
def deleteProbeDBRepeats(inDB, sortCol = None, NR = None, labelCol = None):
    if sortCol == None:
        sortCol = 1
    #NR is if you want to get rid of a probe if it has other copies
    if labelCol <> None:#deal with Labels
        LL = inDB.pop(labelCol)
        
    for i in range(len(inDB)):#For Unsort
        inDB[i].append(i)
        
    Indexfunct = lambda x :x[sortCol].upper().strip()
    indb_BU = FastListBreakup(inDB, Indexfunct)
    indbLen = len(inDB)
    duplicationClusterLen = len(indb_BU)
    if indbLen <> duplicationClusterLen:
        print 'deletduplictes: %i probeObjects in %i duplication Clusters' %(indbLen, duplicationClusterLen)
        outDB = []
        for obj in indb_BU:
            if not(NR):
                outDB.append(obj[0])
    else:
        outDB = inDB
        
    outDB.sort(lambda x,y: cmp(int(x[-1]),int(y[-1])))
    for i in range(len(outDB)):#For Unsort
        outDB[i][0] = i+1
        p = outDB[i].pop(-1)
    if labelCol <> None:#deal with Labels
        outDB = [LL] + outDB
    return outDB

    
def DeleteRepeats(IFN, UD, sortCol, OFNorOutType, INtype = 'fasta', Verbose = 1):
    """DeleteRepeats
    |(IFN, UD, sortCol, OFNorOutType, INtype = 'fasta', Verbose = 1)
    |Output: List or File
    |||Removes repeated sequences in fastaTB file formats. It can accept a list
    (send a list and a dummy for UD) or an in filename. It can either
    return an outfilename or it can return a list. If you want a list return
    specify 'list' for OFNorOutType, and send a list instead of IFN. If you
    want to send a filename, it will default to tab database and detect the
    filename by looking for a '.' within the OutType."""
    verbose = Verbose
    if type(IFN) == type([]):
        INArray = IFN
    else:
        if INtype == 'fasta':
            INArray = FA.InFASTANL(IFN , UD, Verbose = verbose)
        else:
            INArray = InArray(IFN , UD, separator = '\t', Verbose = verbose)
    startLen = len(INArray)
    outArray = DelRepeats(INArray, sortCol, Verbose = verbose)[0]

    if OFNorOutType == 'list':
        return outArray
    else:
        if '.' in OFNorOutType:
            outF = OutArray(outArray, OFNorOutType,UD,Verbose = Verbose)
        elif OFNorOutType == 'fasta':
            OFN = 'nr_' + IFN
            outF = FA.OutFASTA(outArray, OFN, UD, outType = '\n', Verbose = Verbose)
            return outF
        else:
            return 'somethings fuckedup with deleteRepeats'

def ImprovedDeleteRepeats(IFNorLst, UD, sortCol, OFNorOutType, tallyCol = None, #
                          INtype = 'fasta',Verbose = 1):
    """ImprovedDeleteRepeats
    |(IFN, UD, sortCol, OFNorOutType, tallyCol = None, INtype = 'fasta',
    Verbose = 1)
    |Output: List or File
    |||Removes repeated sequences in fastaTB file formats. It can accept a list
    (send a list and a dummy for UD) or an in filename. It can either
    return an outfilename or it can return a list. If you want a list return
    specify 'list' for OFNorOutType, and send a list instead of IFN. If you
    want to send a filename, it will default to tab database and detect the
    filename by looking for a '.' within the OutType. Like DeleteRepeats, slower
    but has more options for tallying the removed repeats"""
    if str(type(IFNorLst)).find('list')<>-1:
        INArray = IFNorLst
    else:
        if INtype == 'fasta':
            INArray = FA.inFASTANL(IFNorLst , UD, Verbose = verbose)
        else:
            INArray = inArray(IFNorLst , UD, separator = '\t',Verbose = verbose)
    startLen = len(INArray)
    subLst = FastListBreakup(INArray,sortCol)
    outArray = []
    tally = None
    for obj in subLst:
        obj.sort(lambda x,y: cmp(x[sortCol],y[sortCol]))
        if tallyCol:
            tally = TallyHo(obj,tallyCol)
            obj.append(tally)
            outArray.append([obj[0],tally])
        else:
            outArray.append(obj[0])
    endLen = len(outArray)
    if Verbose:
        print 'deleteRepeats started with ' + str(startLen) + ' elements and ended up with elements: ' + str(endLen)
    if OFNorOutType == 'list':
        return outArray
    else:
        if '.' in OFNorOutType:
            outF = MB.outDATA(outArray, OFNorOutType,UD,Verbose = Verbose)
        elif OFNorOutType == 'fasta':
            OFN = 'nr_' + IFNorLst
            outF = MB.outFASTA(outArray, OFN, UD, outType = '\n', Verbose = Verbose)
            return outF
        else:
            return 'somethings fuckedup with deleteRepeats'


def DelRepeats(inDB, iColNum,Verbose = 1):
    """DelRepeats
    |(inDB, iColNum)
    |Output: List
    |||Removes repeated sequences in fastaList or other list type. """
    CleanDB = []
    RepDB = []
    NameL = []
    for i in range(len(inDB)):
        #print i
        if not(inDB[i][iColNum] in NameL):
            #print inDB[i][iColNum]
            NameL.append(inDB[i][iColNum])
            CleanDB.append(inDB[i])
        else:
            RepDB.append(inDB[i])
    if Verbose:
        print 'From ' + str(len(inDB)) + ' reduced to ' + str(len(CleanDB)) + '.'
        print len(NameL)
    return [CleanDB , RepDB]


def FastDelRepeats(inDB, iColNum = None,Verbose = 1):
    """FastDelRepeats
    |(inDB, iColNum)
    |Output: List
    |||Removes repeated sequences in fastaList or other list type. """
    CleanDB = []
    RepDB = []
    if iColNum ==None:
        inDB.sort()
    else:
        inDB.sort(lambda x,y:cmp(str(x[iColNum]),str(y[iColNum])))
    if iColNum ==None:
        before = inDB[0]
    else:
        before = inDB[0][iColNum]
    CleanDB.append(inDB[0])
    for i in range(len(inDB)):
        if iColNum ==None:
            piece = inDB[i]
        else:
            piece = inDB[i][iColNum]
        if piece <> before:
            CleanDB.append(inDB[i])
            if iColNum == None:
                before = inDB[i]
            else:
                before = inDB[i][iColNum]
        else:
            RepDB.append(inDB[i]) 
    if Verbose:
        print 'From ' + str(len(inDB)) + ' reduced to ' + str(len(CleanDB)) + '.'

    return [CleanDB , RepDB]


def PaintDB(dbOld, dbAddData, UD, iKeyOldCol, iPosOldCol, iKeyAddCol, iPosAddCol, ColLabel = 'label', startRow = 1):
    """PaintDB
    |(dataSet, dbINF, UD, KeyDBCol, DBCol, KeyDataSetCol, WriteDataSetCol,
    ColLabel = 'label', startRow = 1)
    |Output: List
    |This function takes a dataset writeToDB and a dbIN as input, it must
    have a column with unique Keys, it writes columns from dbIN onto the
    dataSet and returns the new dataset outDB||"""
    indexLst = []
    if type(dbAddData) == type([]):
        dbIN = dbAddData
        dbAddData = 'DataStruct with elements: %s' % (len(dbIN))
    else:
        dbIN = InArray(dbAddData,UD)
    #Build indexer for DB from which we take data
    for obj in dbIN:
        indexLst.append(obj[iKeyOldCol])

    #Write Label
    if startRow > 0:
        diff = iPosOldCol - len(dbOld[0]) + 1
        for z in range(diff):
            dbOld[0].append(' ')
        dbOld[0][iPosOldCol] = ColLabel

    for r in range(1, len(dbOld)):
        obj = dbOld[r]
        if obj[iKeyOldCol] in indexLst:
            i = indexLst.index(obj[iKeyOldCol])
            diff = iPosOldCol - len(obj) +1 
            for z in range(diff):
                obj.append(' ')
            dbOld[r][iPosOldCol] =  dbIN[i][iPosAddCol]
        else:
            print 'PaintDB could not find: %s' % (obj[iKeyOldCol])
    
    print 'PaintDB processed %i lines from database with name: %s' % (len(dbOld), dbAddData)
    return dbCanvas


def DataAccessFunctionProcessor(accessFunct):
    #Usage: accessFunct = AccessFunction(accessFunct)
    if 'function' in str(type(accessFunct)):
        return accessFunct
    if accessFunct == None:
        return lambda x:x
    if 'int' in str(type(accessFunct)):
        exec 'outFunct = lambda x: x[%i]' %accessFunct
        return outFunct
    raise 'error incorrect function sent to AccessFunction %s' %str(accessFunct)





def RankingFunction(inLst_x, inLst_y, accessFunct = None):
    #Usage: Lst_x, Lst_y = RankingFunction(inLst_x, inLst_y, accessFunct = None)
    accessFunct = DataAccessFunctionProcessor(accessFunct)
    Lst_x = []
    Lst_y = []
    Nx = len(inLst_x)
    Ny = len(inLst_y)
    
    if Nx <> Ny:
        raise 'RankingFunction: got two lists that are of different sizes %i and %i' %(Nx, Ny)
    for i in range(Nx):
        Lst_x.append([i, accessFunct(inLst_x[i])])
        Lst_y.append([i, accessFunct(inLst_y[i])])
    Lst_x.sort(lambda x,y: cmp(x[1], y[1]))
    Lst_y.sort(lambda x,y: cmp(x[1], y[1]))
    for i in range(Nx):
        Lst_x[i].append(i)
        Lst_y[i].append(i)
    Lst_x.sort(lambda x,y: cmp(x[0], y[0]))
    Lst_y.sort(lambda x,y: cmp(x[0], y[0]))
    for i in range(Nx):
        Lst_x[i] = Lst_x[i][-1]
        Lst_y[i] = Lst_y[i][-1]
    return Lst_x, Lst_y
        
        
def SpearmanDistance_INPUTRANK(Lst_x, Lst_y):
    #before you package this, add the ranking function here
    if len(Lst_x) <2 or len(Lst_y) <2: return None
    Sumxy = 0.0
    Sumx = 0.0
    Sumy = 0.0
    SumxSqrd  = 0.0
    SumySqrd = 0.0
    if len(Lst_x) == len(Lst_y):
        N = len(Lst_x)
        valsum = 0.0
        for i in range(N):
            Sumxy += Lst_x[i]*Lst_y[i]
            Sumx += Lst_x[i]
            Sumy += Lst_y[i]
            SumxSqrd  += math.pow(Lst_x[i],2)
            SumySqrd += math.pow(Lst_y[i],2)
        SC = (N*Sumxy - Sumx*Sumy)/(math.sqrt(N*SumxSqrd - math.pow(Sumx,2))*math.sqrt(N*SumySqrd - math.pow(Sumy,2)))
    else:
        raise 'Woopsie !! SpearmanRankCorrelation wrong exported number of elements'
    return SC

def SpearmanDistance(inLst_x, inLst_y, accessFunct = None):
    #Usage: SC = SpearmanDistance(Lst_x, Lst_y, accessFunct = lambda x: x.Signal)
    if len(inLst_x) <2 or len(inLst_y) <2: return None
    Lst_x, Lst_y = RankingFunction(inLst_x, inLst_y, accessFunct = accessFunct)
    SC = SpearmanDistance_INPUTRANK(Lst_x, Lst_y)
    return SC





def UnitTest():
    print Average
    print AvgArr
    print DelRepeats
    print DeleteRepeats
    print FindNth
    print ImprovedDeleteRepeats
    print ListBreakup
    print MaMiRec
    print MaxArr
    print MinArr
    print OutArray
    print PaintDB
    print SumArr
    print TallyHo
    print TopPart


if __name__ == '__main__':
    import random
    inlst = [random.random()*10 for i in range(600)]
    inlst.sort()
    for i in [int(random.random()*10) for i in range(20)]:
        b = binsearch(i, 0, len(inlst), inlst)
        print 'BinarySearch: %.2f  result is %i which begets %.2f which is between %.2f and %.2f' %(i, b, inlst[b], inlst[b-1], inlst[b+1])
    Lst_x = [0,20,28,27,50,29, 7, 17, 6, 12]#http://en.wikipedia.org/wiki/Spearman_rank_correlation
    Lst_y = [86, 97, 99, 100, 101, 103, 106, 110, 112, 113]
    SC = SpearmanDistance(Lst_x, Lst_y, accessFunct = None)
    print 'AccessFunction: None SpearmanDistance: %0.4f' %SC
    Lst_x = [[4, Lst_x[i]] for i in range(len(Lst_x))]
    Lst_y = [[4,Lst_y[i]] for i in range(len(Lst_x))]

    SC = SpearmanDistance(Lst_x, Lst_y, accessFunct = 1)
    print 'integer 1: None SpearmanDistance: %0.4f' %SC
    SC = SpearmanDistance(Lst_x, Lst_y, accessFunct = lambda x:x[1])
    print 'lambda: x:x[1]  SpearmanDistance: %0.4f' %SC
    class bullshit:
        def __init__(self, val):
            self.val = val[1]
    Lst_x = [bullshit(Lst_x[i]) for i in range(len(Lst_x))]
    Lst_y = [bullshit(Lst_y[i]) for i in range(len(Lst_x))]
    SC = SpearmanDistance(Lst_x, Lst_y, accessFunct = lambda x:x.val)
    print 'lambda: x:x.val  SpearmanDistance: %0.4f' %SC
    print 'Spearman should be -0.175758'
    
