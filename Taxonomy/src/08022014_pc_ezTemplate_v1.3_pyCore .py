############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################

import os
import sys
import random
##Location
WORK = {'subRoutineDir':'c:/work/Python/New_Core/', #
        'WorkDir' : os.getcwd()+'/'}
CDrive = {'subRoutineDir':'C:\\Users\\dominic\\EclipseWorkspace\\new_pyCore\\src\\',
        'WorkDir' :os.getcwd()+'/',
        }
location = CDrive

subRoutineDir = location['subRoutineDir']

sys.path.append(subRoutineDir[:-1])

#Import Modules
import pyCore
from pyCore.pyLib.__init__ import ImportModulesExecString
execList = ImportModulesExecString()
for l in execList:
    print l
    exec l
    

import Tax.GiToTax as T
import Tax.LabelFasta as LF

location = LOC.LocationObj(locationdict =location)
#location.propagateLocation(modLst)
UD = location.UD
outUD = location.outUD
logUD = location.logUD
appDir = location.appDir



#Logging
errFP = open(logUD + 'log.txt','w')
screenWriter = sys.stdout
logObj = US.logger(errFP,screenWriter,Verbose =1)
sys.stdout = logObj #this diverts the print statement's stream




#############################     New ezTemplate-1.0  091509  #################################
#############################################################################################
#############################     end_EZTEMPLATE HEADER   ###################################
#############################################################################################


verbose = 1

###RUN CONTROL

RunControl = {
              'giToTax':1,
              'faGiToTax':1,
              'FASTA_giToTax':1,
              'FASTA_giToTaxDict':1,
              
              'Names':0,
              'Nodes':1,
              'MakeFinalDB':1,
              
              }


##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/
workingUD= 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy\\'
databaseUD = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy\\'


if RunControl['giToTax']:
    import Tax.LabelFasta as LF
    tIFN = FA.path([databaseUD + 'gi_taxid_nucl\\','gi_taxid_nucl.dmp'],fileType = 'tab')
    txDB = T.BinarySearch(tIFN)

    ifn = '1762_sequence.gi.txt'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    
    LF.giToTax(tIFN, IFN, workingUD)
    

if RunControl['faGiToTax']:
    ##This version just makes a new file with gi to taxID From Fasta File
    
    tIFN = FA.path([databaseUD + 'gi_taxid_nucl\\','gi_taxid_nucl.dmp'],fileType = 'tab')
    txDB = T.BinarySearch(tIFN)

    ifn = 'txID1762_sequence.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    
    LF.fagiToTax(tIFN, IFN, workingUD)
        

    
if RunControl['FASTA_giToTax']:
    ##This makes new file with taxID as part of label
    
    tIFN = FA.path([databaseUD + 'gi_taxid_nucl\\','gi_taxid_nucl.dmp'],fileType = 'tab')
    
    ifn = 'txId1762_sequence.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    
    LF.FASTA_giToTax(tIFN, IFN, workingUD)

if RunControl['FASTA_giToTaxDict']:
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    giTax_IFN = FA.path([databaseUD ,'gi-txid_1762_sequence.gi.txt'],fileType = 'tab')
    GiToTaxDict = giTax_IFN.inDict(0,1)
    ifn = 'txId1762_sequence.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    LF.FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD)


if RunControl['Names'] | RunControl['Nodes']:
    #BuildNames Names Data Structure
    namesIFN = FA.path([databaseUD + 'taxdmp\\','names.dmp'],fileType = 'tab')
    NamesDB = {} #taxID -> scientific name
    for line in namesIFN:
        LLst = line.split('\t')
        if len(LLst) < 6: continue
        if LLst[6] == 'scientific name':
            NamesDB[LLst[0]] = LLst[2]
            
            
if RunControl['Names']:
    #This part adds names to previously made list of this form: [num, gi, tx, newLabel]
    ifn = 'txLabelLst_txId1762_sequence.fasta' #num, gi, tx, newLabel
    IFN = FA.path([workingUD, ifn],fileType = 'tab')
    inDB = IFN.inDB()
    
    for j in range(len(inDB)):
        tx = inDB[j][2]
        if tx in NamesDB:
            name = NamesDB[tx]
        else:
            name = 'notfound'
        inDB[j].append(name)
    OFN = FA.path([workingUD,'txName-'+ifn],fileType = 'tab')
    OFN.outDB(inDB)



if RunControl['Nodes']:
    import Tax.Taxon as Txn
    
    nodesIFN = FA.path([databaseUD+ 'taxdmp\\','nodes.dmp'],fileType = 'tab')
    NodesDB = {}
    for line in nodesIFN:
        LLst = line.split('\t')##131573    |    294823    |    genus    |        |
        if len(LLst) < 5: continue
        NodesDB[LLst[0]] = [LLst[2], LLst[4]]
        
    ifn = 'taxID-nr.txt' #this file has all the taxID's I care about at this time
    IFN = FA.path([workingUD,ifn],fileType = 'tab')
    indb = IFN.inDB()
    TaxonDict = {}
    for txID in indb:
        txID = txID[0]
        if txID in TaxonDict: continue
        t = Txn.Taxon(txID)
        TaxonDict[txID] = t
        while 1:
            if txID == '1': break
            if not(txID) in NodesDB: break
            parentTx, rank = NodesDB[txID]
            name = NamesDB[txID]
            pLst = (txID, name, rank)
            t.process(pLst)
            txID = parentTx
    outLst = []
    for k in TaxonDict.keys():
        outLst.append(TaxonDict[k].display())
    OFN = FA.path([workingUD,'nrTaxNeighborhood.txt'],fileType = 'tab') #This file describes the Taxonomic Neighborhood
    OFN.outDB(outLst)
   
   
    

if RunControl['MakeFinalDB']:
    #This Module will take a Tax Neighborhood file that you curated by hand 
    # You took each taxID and decided wether you care about it or not
    #It will relabel your seqFile using this info
    ifn = '' #[tax    ToDo, other cols . . . . ]
    
    SI_IFN = FA.path([workingUD ,'1762_nrTaxNeighborhood.txt'],fileType = 'tab')
    SIDict = SI_IFN.inDict(0,1)
    
    giTax_IFN = FA.path([workingUD ,'gi-txid_1762_sequence.gi.txt'],fileType = 'tab')
    GiToTaxDict = giTax_IFN.inDict(0,1)
    
    ifn = 'txId1762_sequence.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    LF.FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD, SpecialInstructionDict = SIDict)
    
    
    
              
#######################
#######  CODA   #######
#######################
    
raw_input(' It appears this module ran ok:')        



