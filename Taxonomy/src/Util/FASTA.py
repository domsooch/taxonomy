print 'Import Fasta.py'
import sys
import os
import fileinput
import linecache
import pickle, cPickle
#import anydbm
#import FASTA as FA
#import MolBio as MB
#import BLAST as BL
print 'Yo'
#import ProbeDesignFunctions as PD
#import ListManipulator as LM
#import Useability as US
#import nCluster as CL
import copy
import re
import random
location = {}
Verbose = None
DEBUG = None
import pyCore.pyLib.ListManipulator as LM



IncludedProbes_LL = ['index', 
                'captureprobe', 
                'name', 
                'include', 
                'replicates', 
                'sourcetarget', 
                'notes']
IncludedTargets_LL = ['index', 
                'seq_accession', 
                'seq_custom', 
                'targetname', 
                'probes', 
                'replicates', 
                'notes']



def InLstOrPath(inLstorPATH = None, FileQuery = 'Get Primer File', defaultfn = '*.*', defaulttype = 'fasta'):
    #Certain Functions take either a PATH or List as input this allowes you to send any old thing
    #And then it outputs just the list
    #Perhaps we will write a file iterator in the future
    if (inLstorPATH):
        if 'list' in str(type(inLstorPATH)):
            return inLstorPATH
        elif 'FASTA.path' in str(inLstorPATH):
            IFN = inLstorPATH
            seqLst = None
        else:
            raise 'inLstOrPath(inLstorPATH): inLstorPATH must be either PathObject or List'
    else:
        seqLst = None
        IFN = None
    if not(seqLst):
        if not(IFN):
            if not(location):
                UD = os.getcwd()+'\\'
            else:
                UD = location['WorkDir']
            IFN =  path([UD, ''])
        if not(IFN.pathexists()):
            IFN.EZAsk(title = FileQuery, defaultfn = defaultfn, defaulttype = defaulttype)
        seqLst = IFN.inDB()
    return seqLst


class seqBreakupObj:
    #Usage BUObj = seqBreakupObj(IFN, breakupSize = 10000000, Verbose = 1)
    #BUObj.Breakup()
    #FilePathLstObj = BUObj.returnFilePathLstObj()
    def __init__(self, IFN, breakupSize = 10000000, Verbose = 1):
        self.Verbose = Verbose
        self.breakupSize = breakupSize
        self.IFN = IFN
        self.UD = IFN.d
        self.fileType = IFN.type
        self.ifn = IFN.fn
        self.root = IFN.root
        self.inFP = IFN.FP(RW='r')
        self.outFP = None
        self.FilePATHLst = []
        self.CurrentOFN = None
        self.fileIndex = 0
    def generateNewFile(self):
        if self.CurrentOFN:
            self.outFP.close()
            self.FilePATHLst.append(self.CurrentOFN)
            print 'Closed File %s Size: %skb' %(self.CurrentOFN.fn, str(self.CurrentOFN.getsize())),
        self.fileIndex +=1
        fn = 'part%i_%s' %(self.fileIndex,self.ifn)
        self.CurrentOFN = path([self.UD, fn], fileType = self.fileType)
        self.outFP = self.CurrentOFN.FP(RW='w')
        print '  Opened File %s ' %(self.CurrentOFN.fn)
            
    def Breakup(self, breakupSize = None):
        if breakupSize:
            self.breakupSize = breakupSize
        self.generateNewFile()
        inbuff = ''
        while 1:
            readbuff = self.inFP.read(self.breakupSize)
            if not(readbuff): break
            inbuff = inbuff + readbuff
            if self.Verbose:
                print 'Read %i from %s at pos %s' %(len(inbuff), self.CurrentOFN.fn, str(self.inFP.tell()))
            
            pos = self.FindLastCaret(inbuff)
            if pos == -1:
                self.outFP.write(inbuff)
                inbuff = ''
            else:
                if len(inbuff) < self.breakupSize:
                    self.outFP.write(inbuff)
                    break
                else:
                    self.outFP.write(inbuff[:pos])
                    self.generateNewFile()
                    self.outFP.write(inbuff[pos:])
                    inbuff = ''
        
        self.inFP.close()
        self.outFP.close()
        self.FilePATHLst.append(self.CurrentOFN)
        self.CurrentOFN = None
        self.display()
    def returnFilePathLstObj(self):
        return US.FilePathLstObj(path_Lst =self.FilePATHLst, ud = None, fnLst = [], fileType = self.fileType)
    def FindLastCaret(self, buff):
        if self.Verbose: print 'FindLastCaret:  ',
        pos = -1
        start = self.FindStartSearchPoint(buff)-100
        caretsFound = 0
        
        while 1:
            caretPoint = buff.find('>', start)
            if  caretPoint <> -1:
                start = caretPoint + 1
                pos = caretPoint
                caretsFound +=1
            else:
                break
        if self.Verbose:
            print 'Found %i Carets. Last one at position %i' %(caretsFound, pos)
        return pos
    def FindStartSearchPoint(self, buff):
        for i in range(len(buff), 0, -100):
            testarea = buff[-i:]
            if '>' in testarea:
                print 'found ending caret at %i' %i
                break
        startpos = len(buff) -i
        return startpos
    def display(self):
        for F in self.FilePATHLst:
            print '%s\t\t%i'%(F.fn, F.getsize())
            
 


def MergeFiles(ud, fnLst = None, ofn = None):
    if not(fnLst):
        fnLst = US.ChooseInputFiles(ud)[1]
    if not(ofn):
        ofn = 'merged_.txt'
    OFN = path([ud,ofn],fileType = 'fasta')
    outFP = OFN.FP(RW='w')
    
    for ifn in fnLst:
        print 'Merging ifn: %s '  %(ifn)
        IFN = path([ud,ifn],fileType = 'fasta')
        inFP = IFN.FP(RW='r')
        instr = ' '
        while instr:
            instr = inFP.read(100000000)
            if instr:
                outFP.write(instr)
        inFP.close()
    outFP.close()
    return OFN

def MergeFilesPath(PATHLst, ofn = None, ud = None):
    if not(PATHLst):
        PATHLst = US.FilePathLstObj(fileType = 'fasta')
        PATHLst.askUser()
    if not(ofn):
        ofn = 'merged_.txt'
    if not(ud):
        ud = PATHLst[0].d
    OFN = path([ud,ofn],fileType = 'fasta')
    outFP = OFN.FP(RW='w')
    
    for IFN in PATHLst:
        print 'Merging ifn: %s of size %s'  %(IFN.fn, str(IFN.getsize()))
        inFP = IFN.FP(RW='r')
        instr = ' '
        while instr:
            instr = inFP.read(100000000) + inFP.readline()
            if instr:
                outFP.write(instr)
        inFP.close()
    outFP.close()
    return OFN

def IndexFastaFiles(inPathObjLst = None):
    """MapFastaFile
    |inPathObjLst = None
    |Takes in list of file paths and then return list of every Fasta in file and its location
    """
    ## The reason this sucks is because it uses the flawed lineinput module
    if not(inPathObjLst):
        inPathObjLst = US.FilePathLstObj()
        inPathObjLst.askUser()
    outDB = []
    for inPathObj in inPathObjLst: 
        inp = inPathObj.p
        print 'Processing: ' + inp
        inf = fileinput.FileInput(files = [inp])
        for line in inf:
            i = inf.filelineno()
            if '>' in line[:10]:
                print str(i) + '  ' + line[:-1]
                outDB.append([inp,line[1:-1], i,0])
            if len(outDB) > 0:
                outDB[-1][3] = i
    return outDB


def getFASTAObjects(inList = None, returnAsFastaLstObject = 1):
    """getFASTAObjects
    |inList = None, returnAsFastaLstObject = 1
    |Searches list of fastafiles for specific Fasta Objects
    """
    if not(inList):
        print 'Please choose a set of Fasta Files to choose from'
        tInList = IndexFastaFiles()
        chooseList = []
        indexList = []
        for obj in tInList:
            chooseList.append(obj[1])
            indexList.append(obj[1])
        chosen = US.ChooseFromList(chooseList, Title = 'Choose Fasta Objects')
        inList = []
        for obj in chosen:
            inList.append(tInList[indexList.index(obj)])
        print 'You have chosen %i fasta objects.' %(len(inList))
    outDB = []
    for obj in inList:
        inp, label, fromLine, toLine = obj
        fromLine = int(fromLine)
        toLine = int(toLine)
        print 'Extracting: ' + str(obj[1])
        inf = fileinput.FileInput(files = [inp])
        seq = ''
        for lineNum in range(fromLine+1,toLine+1,1):
            line = linecache.getline(inp,lineNum)[:-1]
            seq = seq + line
        outDB.append([label,seq,len(seq)])
        linecache.clearcache()
    if returnAsFastaLstObject:
        outObj = FASTALstObj(inList = outDB)
        return outObj
    else:
        return outDB



def DumpSuffix(Instr):
    lenStr = len(Instr)
    found = None
    i = 0
    for i in range(lenStr-1,0,-1):
        #print str(i) + '  ' + Instr[i]
        if Instr[i] =='.':
            found = 1
            break
        if not(found): i = lenStr
    retStr = Instr[:i]
    return retStr
  
class path:
    """path"""
    def __init__(self,locationLst=None,fileType = 'tab',Label = None, verify = None):
        self.Label = Label
        self.type = fileType
        self.fn = ''
        self.d = ''
        self.p = ''
        self.root = ''
        self.fdbroot = ''
        self.fp =None
        self.memReadSize = 10000000
        self.outLineBuff = []
        self.SuffixlessPath = ''
        self.outLineBuff = ''
        self.iterfunct = None
        self.timer = None
        locLstType = str(type(locationLst))
        if Verbose: print 'FA.path: inpath: %s  of type %s' %(str(locationLst), locLstType)
        if 'str' in locLstType or 'unicode' in locLstType:
            self.inpath(str(locationLst))
        elif not(locationLst):
            if verify:
                ud, fn = self.askUser()
            else:
                ud, fn = '', ''
            self.allocate(ud,fn)
        elif str(type(locationLst)).find('list')<>-1:
            ud, fn = locationLst
            if verify:
                ud = self.verifyUD(ud)
                ud, fn = self.verifyFN(ud,fn)
            self.allocate(ud,fn)
        elif 'FASTA.path' in str(locationLst):
            ud = locationLst.d
            fn = locationLst.fn
            self.type = locationLst.type
            self.Label = locationLst.Label
            self.allocate(ud,fn)
    def __str__(self):
        return 'FA.path FASTA.path instance :' + self.p
    def inpath(self, pathstr):
        ud = os.path.dirname(pathstr)+'\\'
        fn = os.path.basename(pathstr)
        self.allocate(ud,fn)
    def askUser(self, queryUD = None):
        if not(queryUD):
            queryUD = os.getcwd()+'\\'
        locationLst = US.ChooseInputFile(queryUD,DefaultValue = None)
        self.askFileType()
        ud, fn = locationLst
        self.allocate(ud,fn)
        return ud, fn ##maybe you should return self here
    def askFileType(self):
        retType = raw_input('What type of file is this? (fasta, tab, csv, csa, raw, default is %s):' %self.type)
        if retType:
            self.type = retType
    def EZAsk(self, title = 'title', defaultfn = '', defaulttype = 'tab', Interactive = None):
        self.type = defaulttype
        import ds_easygui as EZ
        inp = EZ.fileopenbox(title=title, default = defaultfn)
        fn = os.path.basename(inp)
        ud = os.path.dirname(inp) + '\\'
        if Interactive:
            self.askFileType()
        self.allocate(ud,fn)
        return ud, fn ##maybe you should return self here
    def FP(self,RW = 'r'): 
        self.fp = open(self.p,RW)
        return self.fp
    def allocate(self,d =None, fn = None):
        if d:
            self.d = d
        if fn:
            self.fn = fn
        if self.d == None: self.d = ''
        if self.fn == None: self.fn = ''
        self.p = self.d + self.fn
        self.root = DumpSuffix(self.fn)
        self.fdbroot = DumpSuffix(self.fn)
        self.SuffixlessPath = DumpSuffix(self.p)
    def changeSuffix(self, newSuffix):
        newfn = self.fdbroot + '.' + newSuffix
        self.allocate(fn = newfn)
    def rename(self, newfn):
        oldpath = os.path.abspath(self.p)
        newpath = self.d + newfn
        newpath = os.path.abspath(newpath)
        if os.path.exists(newpath):
            print 'Newpath exists: %s we ware erasing it' %newpath
            os.remove(newpath)
        os.rename(oldpath, newpath)
        print 'path: renamed %s to %s' %(oldpath, newpath)
        self.fn = newfn
        self.allocate()
    def ERASE(self):
        if self.fp:
            self.fp.close()
        os.remove(self.p)
    def Chunk(self, size = None):
        infp = self.FP(RW = 'r')
        
        if not(size):
            instr = infp.read()
        else:
            instr = infp.read(size) + infp.readline()
        infp.close()
        return instr
        
    def TabChunk(self, size):
        instr = self.Chunk(size)
        inLst = instr.split('\n')
        outDB = []
        for obj in inLst:
            outDB.append(obj.split('\t'))
        return outDB
    def getsize(self):
        return os.path.getsize(self.p)
    def display(self):
        print '\t' + self.p
    def verifyUD(self,ud):
        if os.path.exists(ud):
            return ud
        return os.getcwd()
    def erase(self):
        if self.pathexists():
            print 'Erasing ' + self.p
            os.remove(self.p)
    def pathexists(self):
        if self.type == 'csa':
            return self.csapathexists()
        if self.type =='fdb':
            return self.fdbpathexists()
        if self.fn and os.path.exists(self.p):
            return 1
        else:
            return None
    def csaDict(self):
        csaDict = {}
        if self.type =='csa':
            if self.fn[-4:] != '.csa':
                fn = self.fn + '.csa'
            else:
                fn = self.fn
            indb = LM.InArray(fn,self.d)
        for ob in indb[1:]:
            csaDict[int(ob[0])] =  ob[1].replace('>','')
        return csaDict
    def csaLst(self):
        csaLst = []
        if self.type =='csa':
            indb = LM.InArray(self.fn,self.d)
        for ob in indb:
            csaLst.append(ob[1].replace('>',''))
        return csaLst
    def fdbpathexists(self):
        dbfilenames = []
        for ending in ['.nhr', '.nin', '.nsq']:
            fname = self.fn + ending
            dbfilenames.append(fname)
        for fn in dbfilenames:
            if not(os.path.exists(self.d + fn)):
                print 'FormatDB verification: Did not find %s in Directory %s' %(fn, self.d)
                return None
            else:
                print 'FormatDB Found: %s in Directory %s' %(fn, self.d)
        return 1
    def csapathexists(self):
        if self.type <> 'csa':
            return None
        p = self.p
        if p[-4:].lower() <> '.csa': p = p + '.csa'
        if self.fn and os.path.exists(p):
            return 1
        else:
            return None
    def direxists(self):
        if os.path.exists(self.d):
            return 1
        else:
            return None
    def verifyFN(self,ud,fn):
        if os.path.exists(ud+fn) and fn:
            return ud,fn
        print 'Invalid Path given: ' + ud + fn
        return US.ChooseInputFile(ud,DefaultValue = self.Label)
    def CountLines(self):
        indb = self.inDB()
        LineLen = len(indb)
        print 'FileName: %s has \t\t%i\tLines' %(self.fn, LineLen)
        return LineLen
    def EstimateNumLines(self):
        if not(self.outLineBuff):
            yo = self.__iter__()
            self.ReloadLines()
        filesize = self.getsize()
        bufferReadSize = max(min(filesize, self.memReadSize),1)
        numLines = len(self.outLineBuff)
        estLineLen = int(float(filesize)/bufferReadSize * numLines)
        print 'File %s is %s long and has ~ %s lines' %(self.fn, str(self.getsize()), str(estLineLen))
        return estLineLen
    def selfTimer(self, procName = 'FileReader'):
        ##Part of Iterator machinery
        if self.timer == None:
            totalProcessElements = self.getsize()
            self.timer = US.ReportTimer(totalProcessElements,
                                        ReportInterval = self.memReadSize,
                                        ProcessName = procName,
                                        msec = 1,
                                        maxIntervalsToKeep = 10)
        else:
            filePosition = self.tellfp()
            fileSize = self.getsize()
            if filePosition == 'closed':
                filePosition = fileSize
            self.timer(filePosition)
    ##Iteration group initiator
    def __iter__(self, memReadSize = None): ##That's the function you need in your object 1
        ##Part of Iterator machinery
        ##copy of Create_IterLineObj
        ##Can't figure out how to get looper to take StopIteration error
        ##In reality, the __iter__ returns an iter wrapper object
        ##not waht we are sending with the return self
        ##Fix this later!!
        print '__iter__(self, memReadSize = None)'
        #self.allocate()
        if memReadSize:
            self.memReadSize = memReadSize
        self.iteropen()
        return self##Once you open the fp you pass object to iterator with next and __iter__, it does the rest

    def next(self):##This is the function you need to have in your function 2
        ##Part of Iterator machinery
        if DEBUG and (len(self.outLineBuff) ==0):
            print 'next(self): len of outLineBuff: ' + str(len(self.outLineBuff))
            print 'fploc %s fp: %s' %(str(self.tellfp()), str(self.fp))
        if not(self.outLineBuff):
            self.ReloadLines()
        if self.outLineBuff:
            return self.outLineBuff.pop(0)
        else:
            raise StopIteration ##That's the error you need to raise 3
    def iterclose(self):
        ##Part of Iterator machinery
        if DEBUG: print 'iterclose:'
        if self.fp:
            self.fp.close()
            #self.fp = None if this line were included you would constantly be reopening the same file
    def iteropen(self):
        ##Part of Iterator machinery
        if DEBUG: print 'iteropen:'
        if self.fp:
            if self.fp.closed:
                return None
            else:
                return 1
        elif self.fp == None: #if self.fp is None it means it the fiel was never opened and it will open it otherwise the file would keep getting open
            self.openFile()
            return 1
    def PrepareToReOpenFile(self):
        #This erases the memory that a file has been parsed before and allowes reparsing of the same file
        self.fp = None
        
    def openFile(self):
        ##Part of Iterator machinery
        ##Call this if you explicitly wnat to reopen a file you may have parsed before
        if DEBUG: print 'openFile:'
        self.fp = self.FP(RW = 'r')
    def fillOutLineBuffer(self):
        ##Part of Iterator machinery
        outStr = self.fp.read(self.memReadSize) + self.fp.readline()
        if not(outStr):
            self.iterclose()
            return None
        self.outLineBuff = outStr.split('\n')
        #reporting
        filePosition = self.tellfp()
        fileSize = self.getsize()
        if filePosition == 'closed':
            filePosition = fileSize
        if DEBUG: print 'Reloading %s at position %s' %(self.fn,filePosition)
        return len(self.outLineBuff)
    
    def ReloadLines(self):
        ##Part of Iterator machinery
        #1. if file has not been opene already, open it
        #2. if file is closed return None
        #3. read lines
        #4. if no lines read, then close file  and return None
        if DEBUG: print 'ReloadLines: ' + str(self.fp)
        if self.iteropen():
            self.fillOutLineBuffer()
        self.selfTimer(procName = 'FileLineIterator')
        
    def IterLineObj(self, Kill = None):
        #Alternative Iter ation method: pls deprecate some day
        if Kill:
            self.iterclose()
            linecache.clearcache()
            return 0
        inf = fileinput.FileInput(files = [self.p])
        retFunction = lambda x:linecache.getline(self.p,x)[:-1]
        return retFunction
    def Create_IterLineObj(self, Kill = None, memReadSize = 10000000):
        #Alternative Iter ation method: pls deprecate some day
        self.allocate()
        self.memReadSize = memReadSize
        if Kill:
            self.iterclose()
            self.outLineBuff = None
            return 0
        if not(self.fp):## or (self.fp.closed):
            self.iteropen()
        return self.IterLineFunction
    def IterLineFunction(self):
        #Alternative Iter ation method: pls deprecate some day
        if self.outLineBuff:
            ret = self.outLineBuff.pop(0)
            if ret:
                return ret
        self.ReloadLines()
        if self.outLineBuff:
            return self.outLineBuff.pop(0)
        else:
            return None
    def makeGenBankSeqDict(self):
        #>gi|227809833|gb|FJ966084.1| Influenza A virus (A/California/04/2009(H1N1)) segment 6 neuraminidase (NA) gene, complete cds
        if not(self.type =='fasta'):
            raise 'makeGBSeqDict: Wrong kind of file for this function'
        indb = self.inDB()
        seqDict = {}
        for faObj in indb:
            label = faObj[0]
            labelLst = label.split('|')
            acc = labelLst[3].split('.')[0]
            seqDict[acc] = faObj[1]
        return seqDict
    def inDict(self, keyIndex = 0, valIndex = None):
        if self.type =='fasta':
            keyIndex = 0
            valIndex = 1
        indict = {}
        indb = self.inDB()
        for obj in indb:
            label = obj[keyIndex].strip()
            if valIndex <> None:
                val = obj[valIndex]
            else:
                val = obj
            indict[label] = val
        return indict
    
    def tellfp(self, num = None):
        if num: return self.tellfp_num()
        if not(self.fp):
            return 'neveropened'
        if not(self.fp.closed):
            return self.fp.tell()
        else:
            return 'closed'
    def tellfp_num(self):
        if not(self.fp):
            return 0
        if not(self.fp.closed):
            return self.fp.tell()
        else:
            return self.getsize()
    def PositionInFile(self):
        loc = self.tellfp()
        loctype = str(type(loc))
        if ('long' in loctype) or ('int' in loctype):
            loc = int(loc/1000)
        size = int(self.getsize()/1000)
        p = '%s we are at pos %s of %s' %(self.fn, loc, size)
        return p
    def inDB(self):
        self.allocate()
        print "type is: " + self.type
        if self.type == 'tab':
            inDB = LM.InArray(self.fn,self.d)
            print 'inDB has imported %i lines from file %s of filetype %s' %(len(inDB), self.fn, self.type)
            return inDB
        elif self.type =='fasta':
            inDB = InFASTANL(self.fn,self.d)
            print 'inDB has imported %i lines from file %s of filetype %s' %(len(inDB), self.fn, self.type)
        elif self.type == 'csv':
            inDB = LM.InArray(self.fn,self.d, separator = ',')
            return inDB
        elif self.type == 'csa':
            return self.csaDict()
        elif self.type == 'raw':
            infp = open(self.p, 'r')
            inDB = infp.read()
            infp.close()
            print 'inDB has imported %i chars from file %s of filetype %s' %(len(inDB), self.fn, self.type)
        elif self.type == 'lines':
            infp = open(self.p, 'r')
            inDB = infp.readlines()
            infp.close()
            print 'inDB has imported %i lines from file %s of filetype %s' %(len(inDB), self.fn, self.type)
        elif self.type =='pickle':
            infp = self.FP(RW='r')
            inDB = pickle.load(infp)
            infp.close()
            print 'inDB has imported %i lines from file %s of filetype %s' %(len(inDB), self.fn, self.type)
            return inDB
        elif self.type == 'anydbm':
            inDB = Serializer(self)
            inDB.openRead()##I hope this is ok022207
        return inDB
    def readChunk(self, readSize):
        self.allocate()
        infp = open(self.p, 'r')
        inDB = infp.read(readSize)
        infp.close()
        print 'readChunk has imported %i chars from file %s of filetype %s' %(len(inDB), self.fn, self.type)
        return inDB
    def outDB(self,outDB, prefix = None, ofn = None,suffix = None):
        if prefix:
            self.fn = prefix + self.fn
            
        if ofn:
            self.fn = ofn
            
        if suffix:
            self.fn = self.fdbroot + suffix
        self.allocate()
        if self.type == 'tab':
            writeObj = self.fn
            if self.fp:
                if not(self.fp.closed):
                    writeObj = self.fp
            ofn = LM.OutArray(outDB,writeObj,self.d)
            print 'outDB has exported %i lines to file %s of filetype %s\n in directory: %s' %(len(outDB), self.fn, self.type, self.d)
        elif self.type == 'FPtab':
            if not(self.fp):
                self.fp =FP(self,RW = 'w')
            inDB = LM.OutArray(outDB, self.fp,'l')
            print 'inDB has written %i lines to open file %s of filetype %s at fp %s' %(len(outDB), self.fn, self.type, str(self.tellfp()))
        elif self.type == 'csv':
            ofn = LM.OutArray(outDB,self.fn,self.d,outType = ',')
            print 'outDB has exported %i lines to file %s of filetype %s\n in directory: %s' %(len(outDB), self.fn, self.type, self.d)
        elif self.type =='fasta':
            writeObj = self.fn
            if self.fp:
                if not(self.fp.closed):
                    writeObj = self.fp
            
            ofn = OutFASTA(outDB, writeObj,self.d, outType = '\n')
            print 'outDB has exported %i lines to file %s of filetype %s' %(len(outDB), self.fn, self.type)
        elif self.type == 'csa':
            ofn = LM.OutArray(outDB,self.fn,self.d)
            print 'outDB has exported %i lines to file %s of filetype %s' %(len(outDB), self.fn, self.type)
        elif self.type == 'raw' or self.type == 'lines':
            outfp = open(self.p, 'w')
            t = outfp.write(outDB)
            ofn = self.fn
            outfp.close()
            print 'outDB has exported %i chars to file %s of filetype %s' %(len(outDB), self.fn, self.type)
        elif self.type =='pickle':
            outfp = self.FP(RW='w')
            yo = pickle.dump(outDB,outfp)
            outfp.close()
            print 'inDB has exported %i lines from file %s of filetype %s' %(len(outDB), self.fn, self.type)
        self.allocate()
##    def __getitem__(self, i):
##        #Alternative Iter ation method: pls deprecated
##        if 'slice' in  str(type(i)):
##            pass
##        ##you need to come up with something here using the slice object slice.start, slice.stop, slice.step
##        else:
##            return self.getsingleitem(i)
##    def getsingleitem(self,i):
##        #Alternative Iter ation method: pls deprecated
##        count = 0
##        realcount = 0
##        self.close();self.outLineBuff=[]
##        for line in self:
##            #print 'line: %s i %i realcount %i' %(line, count, realcount)
##            if count == i: break
##            count +=1
##            if line:realcount +=1
##            #Self keeps returning None's until self.fp gets reset and outlinebuffer gets refilled
##            #if you don't do this all hell will break loose as your system will start looping back to the beginning of file
##            
##        if count == realcount:
##            self.close();self.outLineBuff=[]
##            return line
##        else:
##            raise 'Item index %i requested higher than number of lines in file %i' %(i, realcount)
##
####Usage for Iteration function
##BOIter = BO_PATH.Create_IterLineObj()
##
##l = 0; nextl = 1000
##line = ' '
##while line:
##    line = BOIter()
##    l+=1
##    if l > nextl:
##        nextl = l+ 5000
##        p = BO_PATH.PositionInFile()
##        print 'CBO_Parsing line %i %s' %(l, p)
##    line = line.split('\t')
##BO_PATH.Create_IterLineObj(Kill = 1)


class Tabpath(path):
    ##This is designed to return a set of items sharing the same value at lineindex
    #It is designed for Tab-delimited datafiles
    ##It has no implemented disc caching since I can't guarantee that the file will be sorted
    ##perhapps soem ME sorting function can be written one day
    
    def __init__(self,locationLst=None,fileType = 'tab',Label = None, verify = None, lineindex = 0):
        path.__init__(self, locationLst, fileType, Label, verify)
        self.lineIndex = lineindex
        self.AllowReload = 1
    def sort(self):
        if self.getsize() < 100000000:
            inList = self.inDB()
            inList.sort(lambda x,y:cmp(x[self.LineIndex], y[self.lineIndex]))
            self.outDB(inList)
        else:
            raise "File is too big to sort in place write a better function"
            #write a sorting function one day that will sort a large file in place
            #breakup file into 1-10 pieces, then sort each one, then using buffered line reads,
            #place into a final sorted file each next piece from each stream
    def bufferedReloadLines(self):
        pass
        ##write a function that will only read one breakup object at a time from a file
        #a breakup object is an object in which each lineIndex value is equal
    def __iter__(self): ##That's the function you need in your object 1
        self.AllowReload = 1
        self.outLineBuff = []
        return self
    def ReloadLines(self):
        self.AllowReload = None
        inList = self.inDB()
        self.outLineBuff = LM.FastListBreakup(inList, self.lineIndex)
        print 'Tabpath loaded %s There are %i iterable Units in this file.' %(self.fn, len(self.outLineBuff))
    def next(self):##This is the function you need to have in your function 2
        ##Part of Iterator machinery
        if not(self.outLineBuff) and self.AllowReload:
            self.ReloadLines()
        if self.outLineBuff:
            return self.outLineBuff.pop(0)
        else:
            raise StopIteration ##That's the error you need to raise 3
    def __len__(self):
        return len(self.outLineBuff)
    




class workdirObj:
    UD = ''
    installdir = os.getcwd()+'\\'
    
    IFN = path([installdir, 'WorkDir.txt'], fileType = 'raw')
    def __init__(self, UD = ''):
        installdir = self.installdir
        if not(os.path.exists(self.installdir)):
            os.mkdir(self.installdir)
        fnlst = os.listdir(self.installdir)
        if 'WorkDir.txt' in fnlst:
            self.UD = self.get()
        else:
            self.save()
    def save(self, UD = ''):
        self.UD = self.assignUD(UD)
        writestring = self.UD
        self.IFN.outDB(writestring)
        return UD
    def get(self):
        ud = self.IFN.inDB().replace('\n','')
        self.UD = self.assignUD(ud)
        return self.UD
    def assignUD(self, UD = ''):
        if not(UD):
            UD = self.UD
        if UD and os.path.exists(UD):
            self.UD = UD
            return UD
        else:
            print 'defaulting to install directory.'
            UD = self.installdir
            self.UD = UD
        return UD


class folder:
    """folder"""
    def __init__(self,ud=None,Label = None):

        self.d = ud 
    def askUser(self):
        if not(self.d):
            self.d  = US.ChooseInputDirectory(os.getcwd()+'\\',DefaultValue = None)
        else:
            self.d  = US.ChooseInputDirectory(self.d,DefaultValue = None)
    def verifyUD(self):
        if os.path.exists(self.d):
            return self.d
        return os.getcwd()+'\\'
    def display(self):
        print 'folderObject is set to %s' %(str(self.d))

class Serializer:
    def __init__(self,OUT_PATH = None, Overwrite = None, Verbose = 0):
        if not('FASTA.path' in str(OUT_PATH)):
            if 'str' in str(type(OUT_PATH)):
                if '.' in OUT_PATH:
                    ofn = OUT_PATH[:OUT_PATH.find('.')] + '_anydbm.txt'
                else:
                    ofn = OUT_PATH + '_anydbm.txt'
            else:
                ofn =  'Serialized_generic_anydbm.txt'
            UD = os.getcwd()+'\\'
            print 'Making new anydbm: ' + UD + ofn
            OUT_PATH = path([UD,ofn], fileType = 'anydbm')
            if Overwrite:
                OUT_PATH.erase()
        self.Verbose = Verbose
        self.path = OUT_PATH
        self.keyLst = []
        self.n = 0
        self.anydbm = None
        self.anydbmKeys = []
    def __setitem__(self,k,ob):
        self.write(k,ob)
    def __len__(self):
        self.n = len(self.keyLst)
        return self.n
    def write(self, key, obj):
        key = str(key)###This is an experiemnet
        if  self.anydbm == None:
            self.openNew()
        if key in self.keyLst:
            if 0: print 'You are overwriting record %s' %(key)
        else:
            self.keyLst.append(key)
        if self.Verbose :
            print 'Serializer:Copying Object'
        Cobj = copy.deepcopy(obj)
        if self.Verbose :
            print 'Serializer:writing Object'
        self.anydbm[key] = cPickle.dumps(Cobj,1)
        self.n = len(self.keyLst)
    def openNew(self):
        self.anydbm = anydbm.open(self.path.p,'n')
        self.keyLst = self.keys()
    def openRead(self):
        self.anydbm = anydbm.open(self.path.p,'r')
        self.keyLst = self.keys()
        self.n = len(self.keyLst)
    def openReadWrite(self):
        self.anydbm = anydbm.open(self.path.p,'c')
        self.keyLst = self.keys()
        self.n = len(self.keyLst)
    def close(self):
        if self.anydbm:
            self.anydbm.close()
            self.anydbm = None
    def keys(self):
        self.anydbmKeys = self.anydbm.keys()
        return self.anydbmKeys
    def __getitem__(self, key):
        if 'int' in str(type(key)):
            if not(self.anydbmKeys):
                self.keys()
            retObj = cPickle.loads(self.anydbm[self.anydbmKeys[key]])
        else:
            retObj = cPickle.loads(self.anydbm[key])
        return retObj


class newpath:
    """path"""
    def __init__(self,locationLst=None,fileType = None,Label = None):
        self.type = fileType
        if not(locationLst):
            locationLst = US.ChooseInputDirectory(os.getcwd()+'\\')
            self.fn = raw_input('Give me a filename:')
            self.type = raw_input('What type of file is this? (fasta,tab):')
        else:
            ud = locationLst[0].replace('\\','\\')
            fn = locationLst[1]
        if not(os.path.exists(ud)):
            print 'class path: Directory ' + ud + ' does not yet exist.'
        self.d = os.path.abspath(ud)+'\\'
        self.fn = fn
        self.p = self.d + self.fn
        self.root = self.fn[:self.fn.find('.')]
    def askUser(self):
        locationLst = US.ChooseInputFile(os.getcwd()+'\\',DefaultValue = 'FilePath')
        self.type = raw_input('What type of file is this? (fasta,tab):')
        ud, fn = locationLst
        self.allocate(ud,fn)
        return ud, fn
    def allocate(self,d,fn):
        self.d = d
        self.fn = fn
        self.p = self.d + self.fn
        self.root = self.fn[:self.fn.find('.')]




class PDRun:
    ##USAGE
    ##pdrun = PDRun(NumProbes = 12000)
    ##pdrun.doMe(runNum, targetifn)
    probeDB = []
    def __init__(self, DD_IFN = None, PL_IFN = None, Target_IFN = None, NumProbes = 12000, Interactive = None):
        self.Interactive = Interactive
        self.NumProbes = NumProbes
        self.DD_IFN = DD_IFN
        self.PL_IFN = PL_IFN
        self.probeDB = []
        self.Target_IFN = Target_IFN
        self.UD = os.getcwd() + '\\'
        self.ImportData()
        self.NotDesignedList = [] ##['Target Label', 'Sequence', 'Probes Found', Probes Requested']
        self.dumpedduplicates = []
        self.NumProbesDesigned = 0
    def doMe(self, runNum, targetifn):
        ##INPUT PREP
        if not(runNum):
            runNum = raw_input('Give me the 4-digit RunNumber for the PDrun you want to analyze:')
        designDesc_FN = 'Design_Information_File_%s.txt' %runNum
        self.DD_IFN = path([self.UD,designDesc_FN], fileType = 'raw')
        #self.MakeNotDesignedReport()
        probeList_FN = 'Probe_List_%s.csv' %runNum
        self.PL_IFN = path([self.UD,probeList_FN], fileType = 'csv')
        if not(self.PL_IFN.pathexists()):
            probeList_FN = 'Probe_List_%s.txt' %runNum
            self.PL_IFN = path([self.UD,probeList_FN], fileType = 'tab')
        if self.PL_IFN.pathexists():
            self.repairStupidCombiCSVBUG()
            self.PL_IFN.type = 'csv'
        else:
            probeList_FN = 'Probe_List_%s.txt' %runNum
            self.PL_IFN = path([self.UD,probeList_FN], fileType = 'tab')
            self.repairStupidCombiCSVBUG()
            self.PL_IFN.type = 'csv'
        self.Target_IFN = path([self.UD,targetifn], fileType = 'csv')
        if not(self.Target_IFN.pathexists()):
            print '\n\n\nInput Included Target List: \n\n\n'
            self.Target_IFN.askUser()

        ## IMPORT DATA
        self.ImportData()


        if self.Interactive:
            ret = raw_input('Do you want to add more probes? (1 = yes)')
        else:
            ret = None
        if ret ==1:
            addIFN = path([UD,''], fileType = 'tab')
            print '\n\n\nInput Probes you wish to add: \n\n\n'
            addIFN.askUser()
            if addIFN.pathexists():
                self.AddProbes(addIFN,filetype ='tab')
        
        self.DumpDuplicates()
        self.TallyProbes()
        self.ExportNotDesignedList()
        self.MakeIncludedProbes(self.NumProbes)
        print '\n\n\nBLURB:\n'
        print 'For Design id %s' %str(runNum)
        print 'TargetFile %s was submitted.' %self.Target_IFN.fn
        print self.TargetBLURB
        print self.NDReportBLURB
        print self.dupBLURB
        print self.NotDesignedBLURB
        print self.FinalProbeLst_BLURB
        print self.DuplicateFile_BLURB
        print '\n\n\n'

        
        #clusterList = pdrun.ClusterTargetDB()
    def repairStupidCombiCSVBUG(self):
        ##repair retarded Combi BUG
        self.PL_IFN.type = 'raw'
        rawDB = self.PL_IFN.inDB()
        rawDB = rawDB.replace('"','')
        rawDB = rawDB.replace('\t',',')
        self.PL_IFN.outDB(rawDB)
        
        
    def ImportData(self):
        if self.PL_IFN:
            self.probeDB = self.PL_IFN.inDB()
            self.probeLL = self.probeDB.pop(0)
        self.NumProbesDesigned = len(self.probeDB)
        if self.Target_IFN:
            self.TargetDB = self.MakeTargetDB()
            self.UD = self.Target_IFN.d
        
        
        
    def TallyProbes(self):
        for TargObj in self.TargetDB.DataLst:
            TargObj.probesDesigned = 0
        for probeObj in self.probeDB:
            TargetName = probeObj[5]
            self.TargetDB[TargetName].probesDesigned +=1
    def MakeTargetDB(self):
        TargetDB = TargetListDictSeqObj(self.Target_IFN)
        TargetDB.ImportIncludedTargets()
        
        TargetsSubmitted = len(TargetDB)
        ProbesRequested = 0
        TargetsWithProbeRequests = 0
        for TargObj in TargetDB.DataLst:
            ProbesRequested += TargObj.probes
        self.TargetBLURB = '%i Targets were submitted %i Probes were requested %i probes were designed' %(TargetsSubmitted, ProbesRequested, self.NumProbesDesigned)
        return TargetDB
    
    def MakeNotDesignedReport(self):
        ##Extract Data from and parse Description file
        if self.DD_IFN.pathexists():
            DDesc = self.DD_IFN.inDB()
        else:
            self.NDReportBLURB = 'No DesignReportFile Imported'
            return 0
        DDprototype = re.compile(r'Target : \'(.*?)\' found (.*?) of (.*?) requested probes',re.DOTALL)
        TargFindList = DDprototype.findall(DDesc)
        #print str(TargFindList)
        ##('gi|37196787|dbj|AB107208.1| Pisum sativum gene for ITS2', '0', '1')
        self.NotDesignedList = []
        for obj in TargFindList:
            targetname = obj[0]
            targetObject = self.TargetDB[targetname]
            label, acc, seq,  seqlen = targetObject.exportFaLst()
            #print str(obj),
            if obj[1] <> '0': ## Because I'll get it when i export the probelist
                outObj = [label, acc, seq, obj[1], obj[2]]
                self.NotDesignedList.append(outObj) 
        self.NDReportBLURB = 'ProbeDesign could not design all requested probes for %i targets' %len(self.NotDesignedList)
        

    def ExportNotDesignedList(self):
        self.NotDesignedList = []
        self.MakeNotDesignedReport()
        LL = ['Target Label', 'AccNum', 'Sequence', 'Probes Found',
              'Probes Requested']
        
        targsNotDesigned = 0 
        for targObj in self.TargetDB.DataLst:
            if targObj.probesDesigned == 0 and targObj.probes > 0:
                targsNotDesigned +=1
                label, acc, seq, length = targObj.exportFaLst()
                self.NotDesignedList.append([label,acc, seq, 0, targObj.probes])
        self.NotDesignedList.sort(lambda x,y:cmp(x[3], y[3]))
        self.NotDesignedList = [LL] + self.NotDesignedList
        
        ofn = 'Report_' + self.PL_IFN.root + ''
        OFN = path([self.PL_IFN.d, ofn], fileType = 'csv')
        OFN.outDB(self.NotDesignedList)
        self.NotDesignedBLURB = 'Did not design any probes for %i Targets.\n File with full Target Probe acoounting is %s' %(targsNotDesigned, ofn)
    def DumpDuplicates(self):
        lbu = LM.FastListBreakup(self.probeDB,0)
        numDups = len(self.probeDB) - len(lbu)
        #self.MakeTargetsDesigned()
        probeDB = []
        ofn = ''
        self.dupBLURB = 'No duplicates.'
        self.dupBLURB = 'We have %i duplicate clusters in a probe list of %i probes' %(numDups, len(self.probeDB))
        print self.dupBLURB
        
        dumpedduplicates = []
        for obj in lbu:
            probeDB.append(obj[0])
            dumpedduplicates.extend(obj[1:])
        
        ofn = 'dumpedduplicates_' + self.PL_IFN.root + '_T.txt'
        yo = LM.OutArray(dumpedduplicates, ofn,self.UD)
        self.dumpedduplicates = dumpedduplicates
        self.probeDB = probeDB
        self.DuplicateFile_BLURB = 'File with duplicte probes that were removed is: %s' %ofn
    def AddProbes(self, IFN,filetype ='probe'):
        indbObj = ListDictSeqObj(IFN)##This class may be deprecated or non-existant
        if filetype == 'probe':
            indbObj.ImportIncludedProbes()
        else: ##fasta, tab etc that part is handled by the path object
            ##labelIndex = 0, seqIndex = 1, startRow = 0
            indbObj.ImportFASTALst()
        for obj in indbObj.FALst:
            outObj = [obj[1], obj[0],'TRUE',1] + ['' for i in range(17)]
            self.probeDB.append(outObj)
    def MakeIncludedProbes(self, availableElectrodes):
        probeDB = self.probeDB
        repeatsToAdd = availableElectrodes - len(probeDB)
        probeDBlen =  len(probeDB)
        if probeDBlen:
            repeatMultiplier = int(repeatsToAdd/len(probeDB))+1
        else:
            repeatMultiplier = 1
        repeats = max(repeatMultiplier,1)
        numProbes = len(probeDB)
        for i in range(len(probeDB)):
            probeDB[i].append(i)
            probeDB[i][2] = 'TRUE'
            probeDB[i][3] = repeats
        probesChosen = LM.SumArr(probeDB, index = 3)
        repeatsToAdd = availableElectrodes - probesChosen
        print 'repeatMultiplier: %i We started with %i repeats. We have %i probes, \n%i repeats chosen, \n%i repeast to Add' %(repeatMultiplier, repeats, numProbes, probesChosen, repeatsToAdd)
        t = 1;nextt = t
        if numProbes:
            while repeatsToAdd>0:
                if t > 195000:break
                t +=1
                index = random.randrange(numProbes)
                if probeDB[index][3] == repeats:
                    probeDB[index][3] +=1
                    repeatsToAdd -=1
                if t>nextt:
                    nextt = t + 10000
                    print 'Random Cycle %i still have %i repeats to Add' %(t,repeatsToAdd)
        probeDB = [self.probeLL] + probeDB
        outud = self.PL_IFN.d
        ofn = 'Final_' + self.PL_IFN.fn
        self.out_FN = path([outud, ofn], fileType = 'csv')
        self.out_FN.outDB(probeDB)
        self.FinalProbeLst_BLURB = 'Final probe list is in file %s' %self.out_FN.fn
    def ClusterTargetDB(self):
        ##This takes a long time
        TargetList = self.TargetDB.exportFALst()
        faofn = self.Target_IFN.root + '_txt'
        FA_OFN = path([self.Target_IFN.d,faofn], fileType = 'fasta')
        FA_OFN.outDB(TargetList)
        SM, pairDB, fastaLstObj = CL.Seq_GenerateSM(FA_OFN)
        self.clObj = CL.ClusterRun(SM, fastaLstObj, percent=None, unitdistance = 80.0)
        self.clObj.FakeHeirarchicalClustering(targetNumClusters = None, ClusterPercent=None, ClusterUnitdistance = 90.0)
        
        clusterReport = self.clObj.exportClusterReport()
        clofn = 'clusterReport_' + FA_OFN.root + '_T.txt'
        ofn = LM.OutArray(clusterReport,clofn,self.UD)
        
        outFASTALst, clusterArray = self.clObj.AttachLabelToCluster(fastaLstObj)
        clofn = 'clusterArray_' + FA_OFN.root + '_T.txt'
        ofn = LM.OutArray(clusterArray,clofn,self.UD)
        return clusterArray
        


class LabelResolver:
    labelIndexDict = {}
    def __init__(self, LL, LabelLst):
        self.Labels = LabelLst ##Sometimes you only want a subset of the labels
        ##LabelLst = ['index', 'seq_accession', 'seq_custom', 'targetname', 'probes', 'replicates', 'notes']
        for lab in self.Labels:
            if lab in LL:
                self.labelIndexDict[lab] = LL.index(lab)
            else:
                self.labelIndexDict[lab] = None
    def __getitem__(self, item):
        return self.labelIndexDict[item]

class TargetObj:
    probesDesigned = 0
    def __init__(self, inLst, LabelObj, LL):
        self.LR = LabelObj
        self.inLst = inLst
        #LL = ['index', 'seq_accession', 'seq_custom', 'targetname', 'probes', 'replicates', 'notes']
        ##This needs to get rid of " ' and spaces from the labels of both targets and probes. Ideally you want this done before ant after the run
        if len(inLst) > 6:
            for l in LL:
                ##print 'self.%s = inLst[self.LR[\'%s\']]' %(l,l)
                exec 'self.%s = inLst[self.LR[\'%s\']]' %(l,l)
        elif not('index' in LL):
            ##This is how I handle having an included targets file without the index column!!
            for l in LL:
                exec 'self.%s = ""' %(l)
            self.index = 0
            self.seq_custom = inLst[1]
            self.targetname = inLst[2]
            if len(inLst)> 3:
                self.probes = int(inLst[3])
            else:
                self.probes = 0
                
        else:
            #It must mean that the line is not standard size. Like it might be missing probes and notes
            ##This was a bug!! fixed 11-22-07
            for l in LL:
                exec 'self.%s = ""' %(l)
            
            self.seq_custom = inLst[2]
            self.targetname = inLst[3]
            if len(inLst)> 4:
                self.probes = int(inLst[4])
            else:
                self.probes = 0
        
        self.probes = int(self.probes)
        self.probesDesigned = 0
    def __getitem__(self, item):
        #print 'retobj = self.%s' %item
        exec 'retobj = self.%s' %item
        return retobj
    def exportLst(self):
        outLst = []
        for l in ['index', 'seq_accession', 'seq_custom', 'targetname', 'probes', 'replicates', 'notes']:
            outLst.append(self[l])
        return outLst
    def exportFaLst(self): 
        return [self['targetname'], self['seq_accession'], self['seq_custom'], len(self['seq_custom'])]
        
class TargetListDictSeqObj:
    ##Usage
    #TargetDB = ListDictSeqObj(self.Target_IFN)
    #TargetDB.ImportIncludedTargets() or ImportFASTALst
    #or use Make with your own specs
    def __init__(self,IFN):
        self.n = 0
        
        self.LabelDict = {}
        self.DataLst = []
        self.IFN = IFN
        self.indb = self.IFN.inDB()
    def __iter__(self):
        return self.DataLst.__iter__()
        
    
    def __getitem__(self,index):
        if 'int' in str(type(index)):
            return self.DataLst[index]
        elif 'str' in str(type(index)):
            ind = self.LabelDict[index]
            return self.DataLst[ind]
    def __len__(self):
        return len(self.DataLst)
    def exportFALst(self):
        outLst = []
        getmeLst = []
        for tobj in self.DataLst:
            label, acc, seq, slen = tobj.exportFaLst()
            if seq:
                outLst.append([label,seq,slen])
            else:
                getmeLst.append(acc)
        if getmeLst:
            print 'Nee dto get these:'
            for l in getmeLst:
                print l
        return outLst
    def ImportIncludedTargets(self, IFN = None):
        if not(IFN):
            IFN = self.IFN
        inDB = IFN.inDB()
        LL = inDB.pop(0)
        self.LabelResolver = LabelResolver(LL, ['index', 'seq_accession', 'seq_custom', 'targetname', 'probes', 'replicates', 'notes'])
        for obj in inDB:
            addobj = TargetObj(obj, self.LabelResolver, LL)
            self.DataLst.append(addobj)
            self.n = len(self.DataLst)
            self.LabelDict[addobj.targetname] = self.n-1






class FASTAObj:
    """FASTAObj
    |__init__(self,label,seq)
    |def ValidateSeq(self)
    """
    def __init__(self,label,seq,extraData =None):
        self.extra = None
        if extraData:
            self.extra = extraData
        self.index = None
        self.templabel = None
        self.reallabel = label
        self.label = label
        self.accessNum = None
        self.locus = None
        self.slice = None
        self.seq = seq
        self.n = len(seq)
        self.probs = 0
        self.badchars = ''
    def __len__(self):
        return self.n
    def ValidateSeq(self):
        for b in self.seq:
            outSeq = ''
            if b in 'gatcGATC':
                outSeq = outSeq + b
            else:
                self.probs +=1
                self.badchars = self.badchars + b
            #self.seq = outSeq
        return self.probs
    def List(self):
        return [self.reallabel, self.seq, self.n]
    def export(self):
        return [self.reallabel, self.seq, self.n] + self.extra


class BetterFASTALstObj:
    """BetterFASTALstObj
    |__init__(self,inList = None, IFNpath = None)
    |def Add_T(self,path),def Add_FA(self,path),def AddLst(self,inLst),def AddFA(self, inFALst)
    |def AddLst(self,inLst),def AddFA(self, inFALst),def Export_FA(self,outPath)
    |def Export_T(self,outPath),def Lst(self),def Validate(self)
    |
    """
    ##Usage: SubjectDict = BetterFASTALstObj(IFNpath = BGDB_IFN)
    def __init__(self,IFNpath = None, inList = None):
        self.FALst = []#FASTA Objects
        self.IFN_PATH = IFNpath
        if self.IFN_PATH:
            self.ud = IFNpath.d
            self.fn = IFNpath.fn
        else:
            self.ud = ''#location['WorkDir']
            self.fn = ''
        self.mappath = ''
        self.IndexFlag = None
        self.labelDict = {}
        
        if inList:
            self.AddLst(inList)
        if IFNpath:
            inList = IFNpath.inDB()
            self.AddLst(inList)
        self.n = len(self.FALst)
        self.CreateIndexList()
    def __getitem__(self,index):
        if 'int' in str(type(index)):
            return self.FALst[index]
        elif index.isdigit():
            return self.FALst[int(index)]
        else:
            if len(self.labelDict) <> self.n:
                self.remakeLabelDict()
                print 'FastaLstObject Had to remake labelDict'
            ind = self.labelDict[index]
            return self.FALst[ind]
    def __len__(self):
        return len(self.FALst)
    def Add_PATH(self,inpath=None):
        if not(inpath): inpath = path()
        inList = inpath.inDB()
        self.AddLst(inList)
        self.n = len(self.FALst)
    def remakeLabelDict(self):
        self.labelDict = {}
        for j in range(len(self.FALst)):
            faObj = self.FALst[j]
            self.labelDict[faObj.label] = j
        self.n = len(self.labelDict)
    def AddLst(self,inLst, labelIndex = 0, seqIndex = 1):
        for obj in inLst:
            if len(obj) < 3: continue
            label = obj[labelIndex]
            seq = obj[seqIndex]
            addObj = FASTAObj(label,seq, extraData =obj[3:])
            self.FALst.append(addObj)
            index = len(self.FALst)-1
            self.labelDict[label] = index
            self.FALst[-1].index = index
        self.n = len(self.FALst)
        #self.CreateIndexList(self, outpath = None)
    
    def AddFAobj(self, inFALst):
        for obj in inFALst: 
            self.FALst.append(obj)
            self.labelDict[obj[0]] = len(self.FALst)-1
        self.n = len(self.FALst)
    def Export_FA(self,outpath=None):
        if not(outpath):outpath = newpath()
        outdb = self.Lst()
        self.IFN_path.type = 'fasta'
        self.IFN_path.outDB(outdb, prefix = 'Indexed_')
    def Export_T(self,outpath=None):
        if not(outpath):outpath = newpath()
        outdb = self.Lst()
        self.IFN_path.type = 'tab'
        self.IFN_path.outDB(outdb, prefix = 'Indexed_')
    def Lst(self):
        outdb =[]
        for obj in self.FALst:
            outdb.append([obj.label,obj.seq,obj.n])
        return outdb
    def Validate(self):
        for obj in self.FALst:
            ret =obj.ValidateSeq()
            if ret:
                print '%s has %i bad chars: %s' %(obj.label ,ret,obj.badchars)
    def CreateIndexList(self, outpath = None):
        outDB = []
        ofn = 'LabelIndex_' + self.fn + '_T.txt'
        if not(outpath):
            outpath = path([self.ud,ofn],fileType = 'tab')
        else:
            outpath = path(outpath)
        for j in range(len(self.FALst)):
            obj = self.FALst[j]
            if obj.label <> j:#prevents losing the label
                obj.templabel = obj.label
            obj.label = j
            outDB.append([j,obj.templabel,obj.seq,obj.n])
        outpath.outDB(outDB)
        ##ofn = LM.OutArray(outDB,ofn,self.ud,Verbose = 1)
        self.mappath = outpath
        return self.mappath
    def DeIndexList(self,MapPath= None):
        if MapPath:
            inLst = LM.InArray(inpath.fn,inpath.d,separator = '\t')
            self.FALst = []
            self.AddLst(inLst, labelIndex = 1, seqIndex = 2)
        else:
            for obj in self.FALst:
                x = obj.label 
                obj.label = obj.templabel
                obj.templabel = x
        

    
class FASTALstObj:
    """FASTALstObj
    |__init__(self,inList = None, Path_T = None, Path_FA = None)
    |def Add_T(self,path),def Add_FA(self,path),def AddLst(self,inLst),def AddFA(self, inFALst)
    |def AddLst(self,inLst),def AddFA(self, inFALst),def Export_FA(self,outPath)
    |def Export_T(self,outPath),def Lst(self),def Validate(self)
    |
    """
    def __init__(self,inList = None, path_T = None, path_FA = None,extraData = None):
        self.FALst = []#FASTA Objects
        self.fn = ''
        self.ud = ''
        self.mappath = ''
        self.IndexFlag = None
        self.labelDict = {}
        #print str(len(inList))
        if inList:
            if not(extraData):
                self.AddLst(inList)
            else:
                self.AddLstExtraData(inList)
        elif path_T:
            self.fn = os.path.basename(path_T.p)
            self.ud = os.path.dirname(path_T.p)+'\\'
            if not(extraData):
                self.Add_T(path_T)
            else:
                inList = LM.InArray(self.fn,self.ud)
                self.AddLstExtraData(inList)
        elif path_FA:
            self.Add_FA(path_FA)
            self.fn = os.path.basename(path_FA.p)
            self.ud = os.path.dirname(path_FA.p)+'\\'
        self.n = len(self.FALst)
    def __getitem__(self,index):
        if 'int' in str(type(index)):
            return self.FALst[index]
        elif index.isdigit():
            return self.FALst[int(index)]
        else:
            if len(self.labelDict) <> self.n:
                self.remakeLabelDict(self)
                print 'FastaLstObject Had to remake labelDict'
            ind = self.labelDict[index]
            return self.FALst[ind]
    def __len__(self):
        return len(self.FALst)
    def Add_T(self,inpath=None):
        if not(inpath): inpath = path()
        inList = LM.InArray(inpath.fn,inpath.d,separator = '\t')
        self.AddLst(inList)
        self.n = len(self.FALst)
    def Add_FA(self,inpath):
        if not(inpath): inpath = path() 
        inList = InFASTANL(inpath.fn,inpath.d)
        self.AddLst(inList)
        self.n = len(self.FALst)
    def remakeLabelDict(self):
        self.labelDict = {}
        for j in tange(len(self.FALst)):
            faObj = self.FALst[j]
            self.labelDict[faObj.label] = j
    def AddLst(self,inLst, labelIndex = 0, seqIndex = 1):
        for obj in inLst:
            label = obj[labelIndex]
            seq = obj[seqIndex]
            addObj = FASTAObj(label,seq)
            self.FALst.append(addObj)
            self.labelDict[label] = len(self.FALst)-1
        self.n = len(self.FALst)
    def AddLstExtraData(self,inLst, labelIndex = 0, seqIndex = 1):
        for obj in inLst: 
            self.FALst.append(FASTAObj(obj[labelIndex],obj[seqIndex],extraData =obj))
            self.labelDict[label] = len(self.FALst)-1
        self.n = len(self.FALst)
    def AddFAobj(self, inFALst):
        for obj in inFALst: 
            self.FALst.append(obj)
            self.labelDict[obj[0]] = len(self.FALst)-1
        self.n = len(self.FALst)
    def Export_FA(self,outpath=None):
        if not(outpath):outpath = newpath()
        outdb = self.Lst()
        ofn = OutFASTA(outdb,outpath.fn,outpath.d,outType = '\n')
        return ofn
    def Export_T(self,outpath=None):
        if not(outpath):outpath = newpath()
        outdb = self.Lst()
        ofn = LM.OutArray(outdb,outpath.fn,outpath.d,outType = '\t')
        return ofn
    def Lst(self):
        outdb =[]
        for obj in self.FALst:
            outdb.append([obj.label,obj.seq,obj.n])
        return outdb
    def Validate(self):
        for obj in self.FALst:
            ret =obj.ValidateSeq()
            if ret:
                print '%s has %i bad chars: %s' %(obj.label ,ret,obj.badchars)
    def CreateIndexList(self, outpath = None):
        outDB = []
        ofn = 'LabelIndex_' + self.fn + '_T.txt'
        if not(outpath):outpath = newpath([self.ud,ofn],fileType = 'tab')
        
        for j in range(len(self.FALst)):
            obj = self.FALst[j]
            obj.templabel = obj.label
            obj.label = j
            outDB.append([j,obj.templabel,obj.seq,obj.n])
        
        ofn = LM.OutArray(outDB,ofn,self.ud,Verbose = 1)
        self.mappath = path([self.ud,ofn])
        return self.mappath
    def DeIndexList(self,MapPath= None):
        if MapPath:
            inLst = LM.InArray(inpath.fn,inpath.d,separator = '\t')
            self.FALst = []
            self.AddLst(inLst, labelIndex = 1, seqIndex = 2)
        else:
            for obj in self.FALst:
                x = obj.label 
                obj.label = obj.templabel
                obj.templabel = x
        








def renameFASTARecords(genome_IFN=None,map_IFNPATH=None,designRun = None):
    """renameFASTARecords
    |(genome_IFN=None,map_IFNPATH=None)
    |Output: OFNPATH
    |Renames all record in fastafile using the mapfile which is  labelmapfile Itis memory Lean using fileinput
    ||"""
    import fileinput
    
    if designRun:
        prefix = 'Dez_'
    else:
        prefix = 'RN_'
    
    print 'Input Genome Files'

    if not(genome_IFN):
        gpath = 'FA'
    else:
        gpath = genome_IFN.p
        
    if not(os.path.exists(gpath)):
        print 'This Path doesnt exist give me another:' + gpath
        print 'Input Genome Files'
        genome_IFN = US.FilePathLstObj() 
        genome_IFN.askUser()
        genpathLst = []
        basename = genome_IFN[0].fn
        userDir = genome_IFN.d
        for pth in genome_IFN:
            genpathLst.append(pth.p)
    else:
        genpathLst = [genome_IFN.p]
        basename = genome_IFN.fn
        userDir = genome_IFN.d
    

    

    if not(map_IFNPATH):
        mpath = 'FA'
    else:
        mpath = map_IFNPATH.p
    if not(os.path.exists(mpath)):
        print 'Input Map Names File:'
        map_IFNPATH = path()
        map_IFNPATH.askUser()
    
    mapDB = LM.InArray(map_IFNPATH.fn,map_IFNPATH.d)
    mapDict = {}
    for obj in mapDB:
        mapDict[obj[0].strip()] = obj[1].strip()
    if 'FilePathLstObj' in str(genome_IFN):
        rootFNPATH = genome_IFN[0]
    else:
        rootFNPATH = genome_IFN
    print str(rootFNPATH)
    ofnpath = rootFNPATH.d + prefix + rootFNPATH.fn

    
    outFP = open(ofnpath ,'w')
    inf = fileinput.FileInput(files = genpathLst)
    write = None
    if 'RN_' in rootFNPATH.fn:
        accIndex = 1
    else:
        accIndex = 3
    for line in inf:
        if '>' in line:
            i = inf.filelineno()
            print str(i) + '  ' + line
            lLst = line.split('|')
            accNum = lLst[accIndex]
            if '.' in accNum:
                accNum = accNum[:accNum.find('.')]
            if accNum in mapDict:
                newLabel = mapDict[accNum]
            else:
                print 'Found no New Label for ' + accNum
                newLabel = line.replace('\n','')
                newLabel = newLabel.replace('>','')
            if ('Delete' in newLabel and prefix == 'Dez_') or 'D3L3T3' in newLabel: ##D3L3T3 Trumps all other considerations
                write = None
            else:
                write = 1
            line = '>'+ newLabel + '\n'
        if write:
            outFP.write(line)
    outFP.close()
    return path(locationLst=ofnpath)



def FASTACMD(INForLst,DBName,UserDir,subRoutineDir,InParams = {},OFNorFASTAObject = 'FAobject'):
    """FASTACMD
    |(INForLst,DBName,UserDir,subRoutineDir,InParams = {},OFNorFASTAObject = 'FAobject')
    |Output: OFN
    |Very generic Batch File Runner calls createBat
    |fasacmdParams = {'d':'','p':'F','i':'','o':'','l':'10000000','L':'0,0'}|"""
    fastacmdParams = {'d':'','p':'F','i':'','o':'','l':'10000000','L':'0,0'}
    fastacmdParams.update(InParams)
    fastacmdParams['d'] = DBName
    if str(type(INForLst)).find('list')<>-1:
        if len(INForLst) > 1:
            INForLst = LM.OutArray(INForLst,'fastacmdIN.txt',UserDir,outType = '')
            fastacmdParams['i'] = INForLst
        else:
            fastacmdParams.pop('i')
            fastacmdParams['s'] = INForLst[0]
    elif str(type(INForLst)).find('str')<>-1:
        fastacmdParams['i'] = UserDir + INForLst
    else:
        print 'FASTACMD got bad input here are the first 100 chars: \n' + str(INForLst)[0:100]
        return 0
    ## The output will be spit out to the UserDir no matter what since os.chdir runs  
    if OFNorFASTAObject == 'FAobject':
        fastacmdParams['o'] = UserDir+OFNorFASTAObject + '.txt'
    else:
        fastacmdParams['o'] = UserDir+OFNorFASTAObject

    FAOFN = BL.RunBat('fastacmd',UserDir,subRoutineDir,Params=fastacmdParams)
    if OFNorFASTAObject == 'FAobject':
        FAObject = InFASTANL(FAOFN,'')
    else:
        FAObject = FAOFN
    return FAObject

def ConvertFASTA(infilename , userDir, Verbose = 0, inType = '\n'):
    """ConvertFASTA
    |(infilename , userDir, Verbose = 0, inType = '\n')
    |Output: filename of newly converted Fasta file converts to other type
    |Accepting either tab or nl FASTA Files, returns a list of FASTA objects.
    |Input: File, directory, type of FASTA input file.|"""
    response = 'inFASTA:ERRR'
    outfilename = infilename
    if inType == '\n':
        inDB = InFASTANL(infilename, userDir, Verbose)
        if infilename.find('_FA.')<>-1: #outfasta adds _T by itself
            outfilename = outfilename.replace('_FA.','.')
        response = OutFASTA(inDB, outfilename, userDir, outType = '\t', Verbose = 0)
    elif inType == '\t':
        inDB = InFASTATab(infilename , userDir, Verbose)
        if outfilename.find('_T.')<>-1:
            outfilename = outfilename.replace('_T.','_FA.')
            print outfilename + '  from' + infilename
        else:
            index = outfilename.index('.')
            outfilename = outfilename[:index] + '_FA' + outfilename[index:]
        response = OutFASTA(inDB, outfilename, userDir, outType = '\n', Verbose = 0)        
    return response

def StreamConvertFASTANL(infilename , userDir, OFN = None, MemoryReadSize = 100000, Verbose = 0):
    """StreamConvertFASTA
    |(infilename , userDir, OFN = None, MemoryReadSize = 100000, Verbose = 0)
    |Output: filename of newly converted Fasta file converts to other type
    |Accepts only NL type fasta's as inputs stream converts them to tab delimitted files
    |Input: File, directory|"""
    UD = userDir
    if OFN:
        outfilename = OFN
    else:
        outfilename = infilename[:infilename.find('.')] + '_T.txt'
    outFP = open(userDir + outfilename, 'w')
    inFP = open(userDir + infilename, 'r')
    s = ' '
    previousPart = ' ' 
    j = 0
    numRecsProcessed = 0
    while s:
        s = inFP.read(MemoryReadSize)
        s = s + inFP.readline()
        s = s.replace('\r','')
        if '>' in s:
            processString = previousPart + s
            fasta , previousPart = ProcessStringToFasta(processString)
            yo = OutFASTA(fasta, outFP, UD,outType = '\t', Verbose = 0)
            numRecsProcessed = numRecsProcessed + len(fasta)
            fasta = []
        else:
            previousPart = previousPart + s
        if Verbose == 1:
            j+=1
            if j == 0 or j == 10:
                inFileLoc = str(inFP.tell())
                outFileLoc = str(outFP.tell())
                ppLen = str(len(previousPart))
                print '#FASTA\'s Procd: ' + str(numRecsProcessed) + ' InfilePos: ' + inFileLoc + ' OutFilePos: ' + outFileLoc + 'PiggyTailsize:' + ppLen
                j = 0
    if '>' in previousPart:
            fasta, nothing = ProcessStringToFasta(previousPart, End = 1)
            yo = OutFASTA(fasta, outFP, UD,outType = '\t', Verbose = 0)
    inFP.close()
    outFP.close()
    return outfilename


def streamFASTAConvertToDB(genome_IFN=None,DB_PATH=None,LabelProcFunction = None, Verbose = 1):
    """streamFASTAConvertToDB
    |(genome_IFN=None,DB_PATH=None,LabelProcFunction = None, Verbose = 1)
    |Output: PATH of database
    |Accepts only NL type fasta's as inputs stream converts them to database
    |Input: |"""
    MemoryReadSize = 1000000
    if not(LabelProcFunction):
        def procLabel(line):
            accNum  = line.split('|')[3]
            accNum = DumpSuffix(accNum)
            return accNum
        LabelProcFunction = procLabel
    
    if not(genome_IFN):
        genome_IFN = path([os.getcwd()+'\\',''],fileType = 'fasta')
        genome_IFN.askUser()
    if not(DB_PATH):
        dbofn = 'anyDBM_' + genome_IFN.root + '.anydbm'
        DB_PATH = path([genome_IFN.d ,dbofn],fileType = 'anydbm')
        
    DB_writeObj = Serializer(DB_PATH)
    inFP = open(genome_IFN.p,'r')
    
    s = '  '
    previousPart = inFP.read(MemoryReadSize)
    j = 0
    numRecsProcessed = 0
    while s:
        s = inFP.read(MemoryReadSize)
        s = s + inFP.readline()
        s = s.replace('\r','')
        if '>' in s:
            processString = '\n' + previousPart + s
            fasta , previousPart = ProcessStringToFasta(processString)
            for fastaObj in fasta:
                AccNum = LabelProcFunction(fastaObj[0])
                DB_writeObj.write(AccNum,fastaObj)
            numRecsProcessed = numRecsProcessed + len(fasta)
            fasta = []
        else:
            previousPart = previousPart + s
        if Verbose == 1:
            j+=1
            if j == 0 or j == 10:
                inFileLoc = str(inFP.tell())
                recordsWritten = DB_writeObj.n
                ppLen = str(len(previousPart))
                print '#FASTA\'s Procd: ' + str(numRecsProcessed) + ' InfilePos: ' + inFileLoc + ' recordsWritten: ' + str(recordsWritten) + '  PiggyTailsize:' + ppLen
                j = 0
    if '>' in previousPart:
        previousPart = '\n' + previousPart
        fasta, nothing = ProcessStringToFasta(previousPart, End = 1)
        for fastaObj in fasta:
            AccNum = LabelProcFunction(fastaObj[0])
            DB_writeObj.write(AccNum,fastaObj)
    inFP.close()
    DB_writeObj.close()
          
    return DB_writeObj.path



def ProcessStringToFasta(inStr, End = None):
    """ProcessStringToFasta
    |(inStr, previousEnd = None)
    |Careful how you use this, if the first char is not '>' of '\n' then trouble
    |Accepts only NL type fasta's as inputs stream converts them to tab delimitted files
    |Input: File, directory|"""

    inStr = inStr.replace('\n>','\nNEWFasTa_>')
    Li = inStr.split('NEWFasTa_')
    outDB = []
    straggler = ''
    if End:
        Li.append('End')
    for piece in Li[:-1]:
        NLpos = piece.find('\n')
        label = piece[:NLpos]
        sequence = piece[NLpos:]
        sequence = sequence.replace('\n','')
        length = len(sequence)
        ##print label
        if label.find('>') <> -1:
            label = label.replace('>','')
            outDB.append([label,sequence,length])
        else:
            straggler = Li[-1]   
            
    return outDB, straggler

def InFASTA(infilename , userDir, Verbose = 0, inType = '\n'):
    """InFASTA
    |(infilename , userDir, Verbose = 0, inType = '\n')
    |Output: List of FASTA objects.
    |Accepting either tab or nl FASTA Files, returns a list of FASTA objects.
    |Input: File, directory, type of FASTA input file.|"""
    response = 'inFASTA:ERRR'
    if inType == '\t':
        response = InFASTATab(infilename, userDir, Verbose)
    elif inType == '\n':
        response = InFASTANL(infilename , userDir, Verbose)
    return response

def pInFASTANL(path):
    return InFASTANL(path.fn , path.d, Verbose = 0)

def InFASTANL(IFN , userDir, Verbose = 0):
    """InFASTANL
    |Input:(infilename , userDir = None, Verbose = 0)
    |Output: List of FASTA Objects.
    |Either imports a FASTA multi-sequence file or a list of FASTA
    strings and converts to standard working format: [title, sequence,
    seqlength] If UD is sent in as 'lst' then it will operate on the infilename
    as if it was a file input stream, bypassing file opening. If None is sent as userDir
    Then it wil assume that a full path has been speced by the infilename
    |Input: File, Dir|"""
    splitChar = '\n'
    outFasta = []
    titleList = [['',0,0]]
    counter = 0
    
    if userDir.find('string') <> -1:
        s = IFN
    else:
        if str(IFN).find('FASTA.path instance') <> -1:
            infilename = IFN.fn
            userDir = IFN.d
        else:
            infilename = IFN
        infile = open(userDir + infilename, 'r')
        print 'Reading FastaFile'
        s = infile.read()
        infile.close()
    print 'splitting FastaFile'
    Li = s.split(splitChar)
    sequence = ''
    #titleList is [title,startseq, endseq]
    for piece in Li:
    #Kludgy way to deal with biog fasta files assumes '>' is true delimitter
        counter = counter + 1
        if '>' in piece:
            titleList[-1][2] = counter - 1
            titleList.append([piece[1:],counter,0])
    titleList[-1][2] = counter
    for title in titleList[1:]:
        for i in range(title[1],title[2],1):
            sequence = sequence + Li[i]
            sequence = sequence.replace('\r','')
        title[0] = title[0].replace('\r','')
        outFasta.append([title[0],sequence,len(sequence)])
        sequence = ''
    if Verbose == 1:
        if userDir == 'lst': infilename = infilename[0:20]
        print 'inFASTA: pulled ' + str(len(outFasta)) + ' objects From: ' + str(infilename)
    return outFasta

def InFASTANLME(IFN , userDir, Verbose = 0):
    """InFASTANLME
    |Input:(infilename , userDir = None, Verbose = 0)
    |Output: List of FASTA Objects.
    |Either imports a FASTA multi-sequence file or a list of FASTA
    strings and converts to standard working format: [title, sequence,
    seqlength] If UD is sent in as 'lst' then it will operate on the infilename
    as if it was a file input stream, bypassing file opening. If None is sent as userDir
    Then it wil assume that a full path has been speced by the infilename
    |Input: File, Dir|"""
    splitChar = '\n'
    outFasta = []
    titleList = [['',0,0]]
    
    if 'int' in str(type(Verbose)):
        reportInterval = Verbose
    else:
        reportInterval = 0
    NextReportInterval = 0
    if userDir.find('string') <> -1:
        s = IFN
    else:
        if str(IFN).find('FASTA.path instance') <> -1:
            infilename = IFN.fn
            userDir = IFN.d
        else:
            infilename = IFN
    #infile = open(userDir + infilename, 'r')
    inf = fileinput.FileInput(files = [userDir + infilename])
    Li = []
    counter = 0
    
    for line in inf:
        piece = line.replace('\n','')
        counter = counter + 1
        if '>' in piece:
            titleList[-1][2] = counter - 1
            titleList.append([piece[1:],counter,0])
        Li.append(line)
        if Verbose and counter > NextReportInterval:
            NextReportInterval += reportInterval
            print 'Parsed %s to Line %i' %(str(IFN),inf.filelineno())
    if Verbose:
        print 'Putting seqs Together.'
    titleList[-1][2] = counter
    for title in titleList[1:]:
        sequence = ''
        for i in range(title[1],title[2],1):
            sequence = sequence + Li[i]
            sequence = sequence.replace('\r','')
        outFasta.append([title[0],sequence,len(sequence)])
    if Verbose:
        if userDir == 'lst': infilename = infilename[0:20]
        print 'inFASTA: pulled ' + str(len(outFasta)) + ' objects From: ' + str(infilename)
    return outFasta


def InFASTATab(infilename , userDir, Verbose = 0):
    """InFASTATab
    |Input: (infilename , userDir = None, Verbose = 0)
    |Output: List of FASTA Objects.
    |Imports a Tab-delimited FASTA multi-sequence file into standard working
    format: [title, sequence, seqlength] This is for inputting tab-delimited
    FASTA's that are output by me.
    Not getting a userDir assumes explicit path in infilename.
    ||"""
    outFasta = []
    infile = open(userDir + infilename, 'r')
    s = infile.read()
    infile.close()
    Li = s.split('\n')
    for line in Li:
        if (line <> 0 and line <> ''):
            fastRecord = line.split('\t')
            outFasta.append(fastRecord)
    if Verbose == 1:
        print 'inFASTATab: pulled in ' + str(len(outFasta)) + ' fastaFiles'
    return outFasta

def InConcatFASTA(userDir, fileSuffix = 'seq', Verbose = 0):
    """InConcatFASTA
    |(userDir, fileSuffix = 'seq',Verbose = 0)
    |Output: List of FASTA objects.
    |This function concatenates a whole directory of fasta files
    |Input: Directory of FASTA files, suffix to choose for. You must specify
    a suffix.|"""
    filelist = []
    filelist = os.listdir(userDir)
    outFasta = []
    counter = 0
    for infilename in filelist:
        if infilename.find(fileSuffix)<>-1:
            infile = open(userDir + infilename, 'r')
            s = infile.read()
            Li = s.split('\n')
            title = '>' + infilename[0:infilename.find(".")]
            sequence = ""
            for piece in Li:
                    sequence = sequence + piece      
            outFasta.append([title,sequence,len(sequence)])
            infile.close()
            counter = counter + 1
            print 'pulled file: ' + infilename
    if Verbose == 1:
        print 'inConcatFASTA: pulled' + str(counter) + ' files'
    return outFasta

def ConcatFASTA(fileList, userDir, OFN = None, outType = '\n', Verbose = 0):
    """ConcatFASTA
    |(fileList, userDir, OFN = None, outType = '\n', Verbose = 0)
    |Output: Either a list of FASTA objects or a CSV File.
    |This concatenator makes one big file out of a host of little FASTAFiles
    and outputs to a list at the prompt, or a CSV file.
    |Input: File list is [[filname,type],[]], directory, Optional CSV output
    file, Type is '\t' or '\n'.|"""
    outLst = []
    print fileList
    for fileObj in fileList:
        filename = fileObj[0]
        print 'Getting: ' + filename
        sepType = fileObj[1]
        if sepType == '':
            sepType = '\n'
        outPut = InFASTA(filename, userDir, Verbose, sepType)
        outLst.extend(outPut)
    if OFN <> None:
        if OFN.find('_FA') == -1 and outType == '\n':
            OFN = OFN.replace('.','_FA.')
        response = OutFASTA(outLst, OFN, userDir, outType, Verbose)
    else:
        response = outLst
    if Verbose:
        print 'ConcatFASTA has just put together ' + str(len(outLst)) + ' sequences.'
    return response

def IndexInFASTA(infilename , userDir, Verbose = 0):
    """IndexInFASTA
    |(infilename , userDir, Verbose = 0)
    |Output: A list of Lists of the input file items capable of being indexed.
    |Imports a FASTA multi-sequence file into indexable format:
    [[title],[sequence], [seqlength]]
    |Input: file, directory|"""
    outFasta = [[],[],[]] #note distinct structure
    titleList = [['',0,0]]
    counter = 0
    infile = open(userDir + infilename, 'r')
    s = infile.read()
    infile.close()
    Li = s.split('\n')
    sequence = ""
    #titleList is [title, startseq, endseq]
    #Kludgy way to deal with biog fasta files assumes '>' is true delimitter
    for piece in Li:
        counter = counter + 1
        if '>' in piece:
            titleList[-1][2] = counter -1
            titleList.append([piece[1:],counter,0])
    titleList[-1][2] = counter
    for title in titleList[1:]:
        for i in range(title[1],title[2],1):
            sequence = sequence + Li[i]
        outFasta[0].append(title[0])
        outFasta[1].append(sequence)
        outFasta[2].append(len(sequence))
        sequence = ''
    if Verbose == 1:
        print 'IndexInFASTA: pulled ' + str(len(outFasta[0])) + ' files From:' + str(infilename)
    return outFasta

def seqBU(seq, lineLen = 1000, punct = '\n'):
    outStr = ''
    for i in range(0,len(seq),1000):
        outStr = outStr + seq[i:i+1000] + punct
    return outStr

def OutFASTA(dataArray, OFNorFP, userDir, outType = '\t', Verbose = 0):
    """OutFASTA
    |(dataArray, outfilename, userDir, outType = '\t', Verbose = 0)
    |Output: A file of FASTA objects.
    |Exports an array of data to file using either tabs or newline chars
    between pieces.
    |Input: List, File name, Directory, delimit type.|"""
    CloseOFN = 1   
    # OutFasta can accept a filepointer instead of a filename
    if str(type(OFNorFP)).find('file') <> -1:
        if not(OFNorFP.closed):
            outFP = OFNorFP
            CloseOFN = None
        else:
            print 'OutFASTA failed execution, sent bad file pointer: ' + str(outfilename)
            return OFNorFP
    else:
        if outType == '\t':
            outfilename = OFNorFP[0:OFNorFP.find('.')] + '_T' + OFNorFP[OFNorFP.find('.'):]
            OFNorFP = outfilename
        else:
            outfilename = OFNorFP
        outFP = open(userDir + outfilename,'w')
        
    for line in dataArray:
        if outType == '\n':
            outFP.write('>')
            outFP.writelines(str(line[0]) + outType)
            seq = line[1]
            writeseq = seqBU(seq, lineLen = 1000, punct = '\n')
            if len(writeseq) >10:
                pass#print writeseq
            outFP.writelines(writeseq)
        else:
            Parts = len(line)
            for i in range(Parts):
                outFP.writelines(str(line[i]) + outType)
        outFP.writelines('\n')
    if CloseOFN:
        outFP.close()
    if Verbose:
        print 'outFASTA: Wrote ' + str(len(dataArray)) + ' fasta\'s to file: ' + str(OFNorFP)
    return OFNorFP

def SubsetFASTA(subsetListOrFN, inFasta, Verbose = 0):
    """SubsetFASTA
    |(subsetListOrFN, inFasta, Verbose = 0)
    |Output: List of FASTA objects.
    |This is for extracting a subset of FASTA files from a bigger
    database. A list or a listfile is pulled in depending on whether list is
    string or actual list||"""
    outFasta = []
    inFN = 'List Object'
    if type(subsetListOrFN) == type(''):
        inFN = subsetListOrFN
        infile = open(UD + inFN, 'r')
        s = infile.read()
        infile.close()
        subsetListOrFN = s.split('\n')
    for fastaObj in inFasta:
        for item in subsetListOrFN:
            if item <> '':
                if fastaObj[0].find(item)<>-1:
                    outFasta.append(fastaObj)
    if Verbose == 1:
        print 'subsetFASTA: pulled ' + str(len(outFasta)) + ' files From:' #
        + str(inFN)
    return outFasta

def CleanFASTA(inFasta, Verbose = 0):
    """CleanFASTA
    |(inFasta, Verbose = 0)
    |Output: List of FASTA objects.
    |Takes in a list of FASTA objects of standard format ([label, seq, len])
    and makes the Sequence uppercase and free of non-base characters.||"""
    newseq = ''
    reportInterval = 10000
    report = 0
    i = 0
    outFasta = []
    for fastaObj in inFasta:
        seq = fastaObj[1].upper()
        seqLen = len(seq)
        for base in seq:
            i +=1
            if Verbose and (i > report):
                print 'CleanFASTA cleaned ' + str(i) + ' of ' + str(seqLen)
                report = report + reportInterval
            base = base.upper()
            if base in 'GATCNBDHKMNRSVWY':
                newseq = newseq + base
        fastaObj[1] = newseq
        newseq = ''
        outFasta.append(fastaObj)
    return outFasta

def IntervalFASTA(FN, intervalList, UD):
    """IntervalFASTA
    |(FN, intervalList, UD)
    |Output: Series of FASTA Files
    |Function creates a series of smaller FASTA files from a Larger one, based
    on a list of break points provided by the user.|
    intervals: [[2,56]]|"""
    inFasta = inFASTA(FN, UD)
    FN = FN[0:FN.find('.')]
    for interval in intervalList:
        outFa = inFasta[interval[0]:interval[1]]
        outFN = FN+str(interval[0]) + '-' + str(interval[1])+'.txt'
        s = OutFASTA(outFa, outFN, UD, '\n')

def RInFASTANL(infilename , userDir, Verbose = 0):
    """InFASTANL
    |Input:(infilename , userDir = None, Verbose = 0)
    |Output: List of FASTA Objects.
    |Either imports a FASTA multi-sequence file or a list of FASTA
    strings and converts to standard working format: [title, sequence,
    seqlength] If UD is sent in as 'lst' then it will operate on the infilename
    as if it was a file input stream, bypassing file opening. If None is sent as userDir
    Then it wil assume that a full path has been speced by the infilename
    |Input: File, Dir|"""
    splitChar = '\r' ###This is the only difference bewteen INFASTANL !!
    outFasta = []
    titleList = [['',0,0]]
    counter = 0
    if 'lst' == userDir:
        s = infilename
    else:
        infile = open(userDir + infilename, 'r')
        s = infile.read()
        infile.close()
    Li = s.split(splitChar)
    sequence = ''
    #titleList is [title,startseq, endseq]
    for piece in Li:
    #Kludgy way to deal with biog fasta files assumes '>' is true delimitter
        counter = counter + 1
        if '>' in piece:
            titleList[-1][2] = counter - 1
            titleList.append([piece[1:],counter,0])
    titleList[-1][2] = counter
    for title in titleList[1:]:
        for i in range(title[1],title[2],1):
            sequence = sequence + Li[i]
        outFasta.append([title[0],sequence,len(sequence)])
        sequence = ''
    if Verbose == 1:
        if userDir == 'lst': infilename = infilename[0:20]
        print 'inFASTA: pulled ' + str(len(outFasta)) + ' objects From: ' + str(infilename)
    return outFasta





def unitetest():
    interactive = None
    ud = os.getcwd() + '/test/'
    if interactive:
        print 'tabpath'
        tabpath = path()
        ifn_t = 'test_t.txt'
        print 'fapath'
        fapath = path()
        ifn_f = fapath.fn
    else:
        ifn_f = 'test_fa.txt'
        ifn_t = 'test_T.txt'
        
        tabpath = path([ud,ifn_t])
        fapath = path([ud,ifn_f])
    inDB= InFASTANL(ifn_f, ud)
    ifn_T = LM.OutArray(inDB,ifn_t, ud)
    fobj = FASTAObj(inDB[0][0],inDB[0][1])
    fLstobj = FASTALstObj(inList = inDB)
    print 'fLstobj has %i sequences' %fLstobj.n
    fLstobj.Add_T(tabpath)
    print 'fLstobj has %i sequences' %fLstobj.n
    fLstobj.Add_FA(fapath)
    print 'fLstobj has %i sequences' %fLstobj.n
    fLstobj.AddLst(inDB)
    print 'fLstobj has %i sequences' %fLstobj.n
    fLstobj.AddFAobj( [fobj])
    print 'fLstobj has %i sequences' %fLstobj.n
    fLstobj.AddLst(inDB)
    print 'fLstobj has %i sequences' %fLstobj.n
    outPath = newpath([ud,'test1.fasta'])
    fLstobj.Export_FA(outPath)
    print 'fLstobj has %i sequences' %fLstobj.n
    outPath = newpath([ud,'test1_T.txt'])
    fLstobj.Export_T(outPath)
    print 'fLstobj has %i sequences' %fLstobj.n
    testinDB = fLstobj.Lst()
    print 'testinDB has %i sequences' %len(testinDB)
    fLstobj.Validate()
    mappath = fLstobj.CreateIndexList(outpath = None)
    print 'Map file written to: ' + mappath.p
    i = fLstobj.DeIndexList(MapPath= None)
    print 'deindexed inDB has %i sequences' %fLstobj.n

if __name__ == '__main__':
    DEBUG = 1
    ud = os.getcwd() + '/test/'
    IFN = path([ud, 'tCBO_BI-fa-vs-ref_chr5.txt'],fileType = 'tab')
    testindb = IFN.inDB()
    indb = []
    totalProcessElements = IFN.EstimateNumLines()
    t = US.ReportTimer(totalProcessElements, ReportInterval = 1000, ProcessName = 'generic')
    iterator = 0
    for line in IFN:
        iterator +=1
        ret = t(iterator)##On each cycle send iterator
        if line:
            indb.append(line)
    inDB = IFN.inDB()
    before = len(testindb)
    by_iteration= len(indb)
    after= len(inDB)
    if (before == by_iteration)and (by_iteration ==after):
        print 'PASS:',
    else:
        print 'FAIL',
    print 'Size of retrieved items before: %i\nby iteration: %i\nafter: %i' %(len(testindb), len(indb),len(inDB))
    







