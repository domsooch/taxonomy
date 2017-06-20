'''
Created on Sep 6, 2012

@author: dominic
'''
import random
import GiToTax as T
import Util.FASTA as FA
import Nodes

def FASTAToGiLst(inpath):
    #GiToTaxDict = FASTA_giToTaxObjectDict(inpath, txDir)
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    #This is modfified to use the taxon object that can generate a label at the species level
    inFP = open(inpath, 'r')
    giLst = []
    c = 0
    next_inbuff = ''
    while 1:
        inbuff =  inFP.read(10000000) + inFP.readline()
        if not inbuff: break
        inLst = inbuff.split('>')
        for sObj in inLst:
            c+=1
            if not(sObj):continue
            label_pos = sObj.find('\n')
            label = sObj[:label_pos]
            try:
                gi = label.split('|')[1]
            except:
                gi = 'ukgi'
            giLst.append(gi)
    inFP.close()
    return giLst



def giToTax(tIFN, IFN, workingUD):
    #giToTax(tIFN, IFN, workingUD): Takes a file full of gi id's and creates a gi-tax_id file
    print "giToTax: %s using gi_tax: %s"%(IFN.p, tIFN.p)
    ifn = IFN.fn
    txDB = T.BinarySearch(tIFN)
    inFP = IFN.FP('r')
    oDB = []
    c = 0
    NotFound = []
    while 1:
        inbuff = inFP.read(10000000) + inFP.readline()
        print '%i records Read at : %s' %(len(oDB), str(inFP.tell()))
        if not(inbuff):
            inFP.close()
            break
        inLst = inbuff.split('\n')
        for gi in inLst:
            if not(gi):continue
            c+=1
            #print "gi:" + gi
            taxID = txDB[gi]
            if taxID == 'notfound':
                NotFound.append(gi)
            if random.random() > 0.99: print '%i %s \t\t %s' %(c, gi, taxID)
            oDB.append([gi, taxID])
    
    OFN = FA.path([workingUD,'gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(oDB)
    OFN = FA.path([workingUD,'NotFound-gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(NotFound)

def fagiToTax(tIFN, IFN, workingUD):
    ##This version just makes a new file with gi to taxID From Fasta File
    ifn = IFN.fn
    txDB = T.BinarySearch(tIFN)

    inFP = IFN.FP('r')
    oDB = []
    c = 0
    while 1:
        inbuff = inFP.read(10000000)
        print '%i records Read at : %s' %(len(oDB), str(inFP.tell()))
        if not(inbuff):
            inFP.close()
            break
        if not('>' in inbuff): continue
        inLst = inbuff.split('>')
        for sObj in inLst:
            c+=1
            if not(sObj):continue
            label = sObj[:sObj.find('\n')]
            if not('|' in label): continue
            gi = label.split('|')[1]
            
            taxID = txDB[gi]
            if random.random() > 0.99: print '%i: %s \t\t %s' %(c, label[:20], taxID)
            oDB.append([gi, taxID])
    
    OFN = FA.path([workingUD,'gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(oDB)    

    
def FASTA_giToTax(tIFN, IFN, workingUD):
    ##This makes new file with taxID as part of label
    ifn = IFN.fn
    txDB = T.BinarySearch(tIFN)

    inFP = IFN.FP('r')
    
    OFN = FA.path([workingUD,'txid_'+ifn],fileType = 'tab')
    outFP = OFN.FP('w')

    oDB = []
    c = 0
    while 1:
        inbuff = inFP.read(10000000)
        print '%i records Read at : %s' %(len(oDB), str(inFP.tell()))
        if not(inbuff):
            inFP.close()
            break
        if not('>' in inbuff):
            outFP.write(inbuff)
            continue
        inLst = inbuff.split('>')
        for sObj in inLst:
            c+=1
            if not(sObj):continue
            label = sObj[:sObj.find('\n')]
            if not('|' in label):
                outFP.write(sObj)
                continue
            gi = label.split('|')[1]
            
            taxID = txDB[gi]
            if random.random() > 0.99: print '%i: %s \t\t %s' %(c, label[:20], taxID)
            wbuff = '>%s|%s' %(taxID, sObj)
            outFP.write(wbuff)
            oDB.append([gi, taxID])
    
    OFN = FA.path([workingUD,'gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(oDB)
    
    inFP.close()
    outFP.close()
    
def FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD, SpecialInstructionDict = {}):
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    ifn = IFN.fn

    inFP = IFN.FP('r')
    
    OFN = FA.path([workingUD,'txid_'+ifn],fileType = 'tab')
    outFP = OFN.FP('w')
    
    LabelLst = []
    oDB = []
    c = 0
    while 1:
        inbuff = inFP.read(10000000)
        print '%i records Read at : %s' %(len(oDB), str(inFP.tell()))
        if not(inbuff):
            inFP.close()
            break
        if not('>' in inbuff):
            outFP.write(inbuff)
            continue
        inLst = inbuff.split('>')
        for sObj in inLst:
            c+=1
            if not(sObj):continue
            label = sObj[:sObj.find('\n')]
            if not('|' in label):
                outFP.write(sObj)
                continue
            gi = label.split('|')[1]
            if not(gi):
                taxID = 'unknown'
            else:
                taxID = GiToTaxDict[gi]
            if SpecialInstructionDict:
                si = SpecialInstructionDict[taxID]
                if si =='dump': continue
                else: taxID = si
            if random.random() > 0.99: print '%i: %s \t\t %s' %(c, label[:20], taxID)
            wbuff = '>%s|%s' %(taxID, sObj)
            newLabel = wbuff[:wbuff.find('\n')]
            LabelLst.append([c, gi, taxID, newLabel])
            outFP.write(wbuff)
            oDB.append([gi, taxID])
    
    OFN = FA.path([workingUD,'gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(oDB)
    OFN = FA.path([workingUD,'txLabelLst_'+ifn],fileType = 'tab')
    OFN.outDB(LabelLst)
    
    inFP.close()
    outFP.close()
