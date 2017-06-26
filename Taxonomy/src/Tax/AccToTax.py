print 'Import GiToTax.py'
import os
import sys

##This is waste!! The file from ncbi does not come sorted, so without sorting no binary search will work.

txDir = "J:/taxonomy/"
acctax_path = os.path.join(txDir,'nucl_gb.accession2taxid')


if os.path.exists(acctax_path):
    print ' you are in a place where this taxid thing will work'
    print 'You have access to %s' %acctax_path
    goforlaunch = 1
else:
    print 'Sorry Charlie no access to ' + acctax_path
    goforlaunch = False
    
    
##Make sure all the files especially acc_taxid_nucl are in taxdmp
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
    def EstimateLineLength(self, fp):
        fp.seek(0)
        fp.readline()
        b = fp.readline()
        self.EstMaxLineLength = len(b)+2
    def fileopen(self):
        if not(self.fp) or self.fp.closed:
            print 'opening file . . . %i'%self.fileSize
            self.fp = open(self.inpath, 'r')
            self.EstimateLineLength(self.fp)
    def __getitem__(self, acc):
        self.fileopen()
        acc = acc.split('.')[0]
        #acc = acc.replace('NZ_','')
        if acc in self.d:
            return self.d[acc]
        taxdict = self.BinaryChoice(acc)
        if acc in taxdict:
            self.d[acc] = taxdict[acc]
            return self.d[acc]##taxid
        else:
            return 'notfound'
    def BinaryChoice(self, acc):
        if self.Verbose: print 'BinaryChoice: acc: %s' %(acc)
        delta = self.fileSize/2
        center_fpos = 0
        while abs(delta) > 1000:
            center_fpos = center_fpos + delta
            center_acc = self.accAtPos(center_fpos)[0]
            if self.Verbose: print 'BinaryChoice: acc %s center_acc: %s center_fpos %i delta %i' %(acc, center_acc, center_fpos, delta)
            if acc > center_acc:
                if self.Verbose: print '\tacc %s is > center_acc %s.' %(acc, center_acc)
                delta = abs(delta)/2 + self.EstMaxLineLength
            else:
                if self.Verbose: print '\tacc %s is < center_acc %s.' %(acc, center_acc)
                delta = -1*abs(delta)/2 - self.EstMaxLineLength
            
        left = center_fpos
        right = center_fpos + delta*2
        accdict = self.readInterval(left, right, acc)
        return accdict
    def readBuff(self, start, end):
        #if self.Verbose: print 'readBuff(self, start_%i , end_%i)' %(start, end)
        self.fileopen()
        start = max(0, start- self.EstMaxLineLength*2)
        self.fp.seek(start)
        self.fp.readline()
        readsize = max(1, end-start)
        #if self.Verbose: print 'After seek: Filepos: %i' %self.fp.tell()
        readBuff = self.fp.read(readsize)
        #if self.Verbose: print 'readBuff: ' + readBuff
        if not(readBuff): return ''
        if not(readBuff[-1] =='\n'):
            readBuff = readBuff + self.fp.readline()
            #if self.Verbose: print 'readBuff: read readline: ' + readBuff
        if self.Verbose: print 'BuffRead: Filepos: %i bytesRead: %i' %(self.fp.tell(), readsize)#, readBuff
        buffLst = readBuff.split('\n')
        if buffLst[-1] =='': buffLst = buffLst[:-1]##because when you split a line with a trailing \n, you get an empty entr at the end
        return buffLst
    def readInterval(self, inleft, inright, acc):
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
            ##accession<TAB>accession.version<TAB>taxid<TAB>gi
            if len(lineLst) ==4:
                tacc, taccv, taxid, gi= lineLst
                if tacc =='': continue
                if taxid =='':continue
                
                if tacc > acc:
                    if self.Verbose: print '\t\tStopped at %s because it is greater than %s Is acc in tempdicty?: %s' %(tacc, acc, str(acc in self.testdict))
                    break
                try:
                    taxid = int(taxid)
                except:
                    print 'badline: ', lineLst
                    continue
                self.testdict[tacc] = int(taxid)
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
    def LineToacc(self, line):
        # accession<TAB>accession.version<TAB>taxid<TAB>gi
        if not(line): return 'notfound'
        line = line.replace('\n','')
        lineLst = line.split('\t')
        if len(line) <2: return 'notfound'
        acc = lineLst[0]
        return acc
    def LineToTaxID(self, line):
        if line == '': return 'notfound'
        line = line.replace('\n', '')
        lineLst = line.split('\t')
        gi = lineLst[0]
        tax  = int(lineLst[2])
        return tax
    def accAtPos(self, fpos):
        line, fpos = self.readWholeLineAt(fpos)
        return self.LineToacc(line), fpos
    def readfile(self):
        self.fp.seek(0)
        ret = ''
        while ret !='x':
            buff =self.fp.read(10000) + self.fp.readline()
            l= buff.split('\n')
            for r in l:
                print r
            ret = raw_input('x to quit:')
    

        

if __name__ == '__main__':
    testLst = ['NC_001911.1','CP009253.1', 
            'CP011299.1', 
            'CP013259.1', 
            'NC_002528.1', 
            'NC_004061.1', 
            'NC_004545.1', 
            'NC_008513.1', 
            'NC_011833.1', 
            'NC_011834.1', 
            'NC_015662.1', 
            'NC_017252.1', 
            'NC_017253.1', 
            'NC_017254.1', 
            'NC_017255.1', 
            'NC_017256.1', 
            'NC_017259.1', 
            'NZ_CP002697.1', 
            'NZ_CP002699.1', 
            'NZ_CP002701.1', 
            'NZ_CP002703.1', 
            'NZ_CP009253.1', 
            'NZ_CP011299.1', 
            'NZ_CP013259.1', 
            'NZ_LN890285.1', 
            'NZ_LT635893.1', 
            'NZ_LT667500.1', 
            'NZ_LT667503.1', 
            'NC_001910.1', 
            'NC_001911.1']
    
    acc2tax =  BinarySearch(acctax_path)
    acc2tax.Verbose =1
    outdb = []
    
    for acc in testLst:
        taxid = acc2tax[acc]
        print '%s \t has  a taxid \t%s' %(acc, str(taxid))
        if taxid <> 'notfound':
            outdb.append([acc,taxid, acc2tax.fp.tell()])
    print outdb
    acc2tax.readfile()
    
else:
    #Instantiate Object
    if goforlaunch:
        bs = BinarySearch(acctax_path)
        print 'Access functionalitty by calling TX.bs[Number of gi]'



        
