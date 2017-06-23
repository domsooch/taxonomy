print 'Import GiToTax.py'
import os
import sys



txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\taxdmp\\'
inpath = os.path.join(txDir,'gi_taxid_nucl.dmp')


if os.path.exists(inpath):
    print ' you are in a place where this taxid thing will work'
    print 'You have access to %s' %inpath
    goforlaunch = 1
else:
    print 'Sorry Charlie no access to ' + inpath
    goforlaunch = True
    
    
##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/


class BinarySearch:
    ##this assumes that the file is an ordered list of intergers or floats
    #it was designed for taxdump gi_taxid_nucl.dmp from ncbi that has all gi's mapped to taxid's
    def __init__(self, inpath ,Verbose = 0):
        self.Verbose = Verbose
        self.EstMaxLineLength = 100
        self.inpath = inpath
        self.fp = None
        self.fileSize = os.path.getsize(inpath)
        self.d = {}
    def fileopen(self):
        if not(self.fp) or self.fp.closed:
            print 'opening file . . . %i'%self.fileSize
            self.fp = open(self.inpath, 'r')
    def __getitem__(self, gi):
        self.fileopen()
        try:
            gi = int(gi)
        except:
            print 'Cannot convert this %s'%gi
        if gi in self.d:
            return self.d[gi]
        taxdict = self.BinaryChoice(gi)
        if gi in taxdict:
            self.d[gi] = taxdict[gi]
            return self.d[gi]##taxid
        else:
            return 'notfound'
    def BinaryChoice(self, gi):
        if self.Verbose: print 'BinaryChoice: gi %i' %(gi)
        delta = self.fileSize/2
        center_fpos = 0
        while abs(delta) > 1000:
            
            center_fpos = center_fpos + delta
            center_gi = self.giAtPos(center_fpos)[0]
            if self.Verbose: print 'BinaryChoice: gi %i center_fpos %i delta %i' %(gi, center_fpos, delta)
            if gi > center_gi:
                if self.Verbose: print '/tgi %i is > cebter_gi %i.' %(gi, center_gi)
                delta = abs(delta)/2 + self.EstMaxLineLength
            else:
                if self.Verbose: print '/tgi %i is < cebter_gi %i.' %(gi, center_gi)
                delta = -1*abs(delta)/2 - self.EstMaxLineLength
            
        left = center_fpos
        right = center_fpos + delta*2
        gidict = self.readInterval(left, right, gi)
        return gidict
    def readBuff(self, start, end):
        if self.Verbose: print 'readBuff(self, start_%i , end_%i)' %(start, end)
        self.fileopen()
        start = max(0, start- self.EstMaxLineLength)
        self.fp.seek(start)
        es = end-start
        readsize = max(1, end-start)
        #if self.Verbose: print 'After seek: Filepos: %i' %self.fp.tell()
        readBuff = self.fp.read(readsize)
        #if self.Verbose: print 'readBuff: ' + readBuff
        if not(readBuff): return '0\t0\n0\t0\n0\t0\n'
        if not(readBuff[-1] =='\n'):
            readBuff = readBuff + self.fp.readline()
            #if self.Verbose: print 'readBuff: read readline: ' + readBuff
        if self.Verbose: print 'After BuffRead: Filepos: %i bytesRead: %i' %(self.fp.tell(), readsize)
        buffLst = readBuff.split('\n')
        if buffLst[-1] =='': buffLst = buffLst[:-1]##because when you split a line with a trailing \n, you get an empty entr at the end
        return buffLst
    def readInterval(self, inleft, inright, gi):
        #print '\treadInterval: %i to %i' %(inleft, inright)
        left = min(inleft, inright)
        right = max(inleft, inright)
        left = max(left, 0)
        right = min(right, self.fileSize)
        #print 'Final readInterval: %i to %i' %(left, right)
        self.testdict = {}
        buffLst = self.readBuff(left, right)
        for line in buffLst:
            lineLst = line.split('\t')
            if len(lineLst) ==2:
                tgi, taxid = lineLst
                if tgi =='': continue
                if taxid =='':continue
                tgi = int(tgi)
                if tgi > gi:
                    if self.Verbose: print '\t\tStopped at %i because it is greater than %i Is gi in tempdicty?: %s' %(tgi, gi, str(gi in self.testdict))
                    break
                self.testdict[tgi] = int(taxid)
            else:
                if self.Verbose: print 'BadLine: ' + str(lineLst)
        return self.testdict    

    def readWholeLineAt(self, fpos):
        start = max(0, fpos - self.EstMaxLineLength)
        end = fpos
        line_Lst = self.readBuff(start, end)
        line = line_Lst[-1]
        fpos = self.fp.tell() - len(line)-2#one for the missing\n one for zerobased
        return line, fpos
    def LineToGi(self, line):
        if not(line): return 'notfound'
        line = line.replace('\n','')
        lineLst = line.split('\t')
        if len(line) <2: return 'notfound'
        gi = int(lineLst[0])
        tax  = int(lineLst[1])
        return gi
    def LineToTaxID(self, line):
        if line == '': return 'notfound'
        line = line.replace('\n', '')
        lineLst = line.split('\t')
        gi = int(lineLst[0])
        tax  = int(lineLst[1])
        return tax
    def giAtPos(self, fpos):
        line, fpos = self.readWholeLineAt(fpos)
        return self.LineToGi(line), fpos
    

        

if __name__ == '__main__':
    txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\taxdmp\\'
    inpath = os.path.join(txDir,'gi_taxid_nucl.dmp')
    yo = BinarySearch(inpath)
    yo.Verbose =0
    outdb = []
    totalLines = 10000#IFN.EstimateNumLines()
    #t = US.ReportTimer(totalLines, ReportInterval = 1000, ProcessName = 'taxindexer')
    #yo.Verbose = 1;
    yo[2516]
    
    for i in range(10000):
        #t(i)
        taxid = yo[i]
        print '%i \t has  a taxid \t%s' %(i, str(taxid))
        if taxid <> 'notfound':
            outdb.append([i,taxid, yo.fp.tell()])
    print 633258064, 'test >', yo[633258064]
else:
    #Instantiate Object
    if goforlaunch:
        bs = BinarySearch(inpath)
        print 'Access functionalitty by calling TX.bs[Number of gi]'



        