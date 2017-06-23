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
              'giToTax':0,
              'faGiToTax':0,#Gi to taxid file from fasta
              'FASTA_giToTax':0, #useless
              'FASTA_giToTaxDict':0, #useless
              
              'Names':0, #useless
              'Nodes':0,#this makes neighborhood file that you edit then you run . . . .
              'MakeFinalDB':1, #this guy
              
              }


##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/

TaxDB_UD = 'M:\\bacteria_refseq\\taxonomy\\'
workingUD = "M:\\bacteria_refseq\\"
outUD = "M:\\bacteria_refseq\\"+"\\Output\\"
gi_taxid_IFN = FA.path([TaxDB_UD, 'gi_taxid_nucl.dmp'], fileType = 'tab')

#databaseUD = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy\\'

if RunControl['giToTax']:
    import Tax.LabelFasta as LF
    

    ifn = 'test-giList.txt'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    #giToTax(tIFN, IFN, workingUD): Takes a file full of gi id's and creates a gi-tax_id file
    LF.giToTax(gi_taxid_IFN, IFN, outUD)
    

if RunControl['faGiToTax']:
    ##This version just makes a new file with gi to taxID From Fasta File
    
    txDB = T.BinarySearch(gi_taxid_IFN)

    ifn = 'test_fasta.fasta'
    
    ifn = 'all_safeIntervals.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    #fagiToTax: This version just makes a new file with gi to taxID From Fasta File
    LF.fagiToTax(gi_taxid_IFN, IFN, outUD)
        

    
if RunControl['FASTA_giToTax']:
    ##This makes new file with taxID as part of label

    ifn = 'alldezseqs.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    
    LF.FASTA_giToTax(gi_taxid_IFN, IFN, outUD)

if False and RunControl['FASTA_giToTaxDict']:
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    giTax_IFN = FA.path([databaseUD ,'gi-txid_1762_sequence.gi.txt'],fileType = 'tab')
    GiToTaxDict = giTax_IFN.inDict(0,1)
    ifn = 'txId1762_sequence.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    LF.FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD)


if RunControl['Names'] | RunControl['Nodes']:
    print 'Names or Nodes'
    #BuildNames Names Data Structure
    namesIFN = FA.path([TaxDB_UD + 'taxdump\\','names.dmp'],fileType = 'tab')
    NamesDB = {} #taxID -> scientific name
    for line in namesIFN:
        LLst = line.split('\t')
        if len(LLst) < 6: continue
        if LLst[6] == 'scientific name':
            NamesDB[LLst[0]] = LLst[2]
            
            
if False and RunControl['Names']:
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
    print 'Nodes:'
    import Tax.Taxon as Txn
    
    nodesIFN = FA.path([TaxDB_UD + 'taxdump\\','nodes.dmp'],fileType = 'tab')
    NodesDB = {}
    for line in nodesIFN:
        LLst = line.split('\t')##131573    |    294823    |    genus    |        |
        if len(LLst) < 5: continue
        NodesDB[LLst[0]] = [LLst[2], LLst[4]]
        
    ifn = 'gi-txid_all_safeIntervals.fasta' #this file has all the taxID's I care about at this time
    root = 'safeintervals'
    IFN = FA.path([outUD,ifn],fileType = 'tab')
    indb = IFN.inDB()
    TaxonDict = {}
    for txID in indb:
        txID = txID[1]
        if txID in TaxonDict: continue
        t = Txn.Taxon(txID)
        TaxonDict[txID] = t
        while 1:
            if txID == '1': break
            if not(txID) in NodesDB: 
                print 'NodesDB is Missing txID %s'%txID
                break
            parentTx, rank = NodesDB[txID]
            name = NamesDB[txID]
            pLst = (txID, name, rank)
            t.process(pLst)
            txID = parentTx
    outLst = []
    for k in TaxonDict.keys():
        outLst.append(TaxonDict[k].display())
    OFN = FA.path([outUD,'nrTaxNeighborhood_%s.txt'%root],fileType = 'tab') #This file describes the Taxonomic Neighborhood
    OFN.outDB(outLst)
   
   
    

if RunControl['MakeFinalDB']:
    #This Module will take a Tax Neighborhood file that you curated by hand 
    # You took each taxID and decided wether you care about it or not
    #It will relabel your seqFile using this info
    ifn = '' #[tax    ToDo, other cols . . . . ]
    
    SI_IFN = FA.path([outUD ,'edited_nrTaxNeighborhood_safeintervals.txt'],fileType = 'tab') #This is where you add the edited neighborhood file
    SIDict = SI_IFN.inDict(0,1)
    
    giTax_IFN = FA.path([outUD ,'gi-txid_all_safeIntervals.fasta'],fileType = 'tab')#add the mappings of gi to taxid
    
    GiToTaxDict = giTax_IFN.inDict(0,1)
    
    ifn = 'all_safeIntervals.fasta'  #'alldezseqs.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    LF.FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD, SpecialInstructionDict = SIDict)
    
    
    
              
#######################
#######  CODA   #######
#######################
    
raw_input(' It appears this module ran ok:')        



