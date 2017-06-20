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
Apple = {'subRoutineDir':'/Users/dominicsuciu/Dropbox/genarraytion/new_pyCore/src/',
         'WorkDir' :os.getcwd()+'/',}

location = CDrive
TaxDB_Dir = "C:/Users/dominic/Documents/Worki7/Taxonomy/" #this is where current tax db is located
Names_path = os.path.join(TaxDB_Dir, 'taxdmp/', 'names.dmp')
Nodes_path = os.path.join(TaxDB_Dir, 'taxdmp/', 'nodes.dmp')
gi_tax_path = os.path.join(TaxDB_Dir , "gi_taxid_nucl/", "gi_taxid_nucl.dmp")

subRoutineDir = location['subRoutineDir']
UD = location['WorkDir']

sys.path.append(subRoutineDir)

#Import Modules
import pyCore
from pyCore.pyLib.__init__ import ImportModulesExecString
execList = ImportModulesExecString()
for l in execList:
    print l
    exec l
    

import Tax.GiToTax as T
import Tax.LabelFasta as LF

#location = LOC.LocationObj(locationdict =location)
#location.propagateLocation(modLst)
UD = UD
outUD = UD
logUD = UD
appDir = UD



#Logging
if None:
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
              'genarraytion':0,
              'giToTax':0,
              'faGiToTax':0,
              'FASTA_giToTax':0,
              'FASTA_giToTaxDict':0,
              
              'Names':0,
              'Nodes':1,
              'MakeFinalDB':0,
              
              }


##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/
workingUD= 'M:/bacteria_refseq/genarraytion_taxid/'
databaseUD = "C:/Users/dominic/Documents/Worki7/Taxonomy/"
taxdmp_databaseUD = "C:/Users/dominic/Documents/Worki7/Taxonomy/taxdmp/"


#Take all gi's I care about: 
genarraytion_txid_path = FA.path("M:/bacteria_refseq/genarraytion_taxid/genarraytion_gi_taxid_nucl.dmp")
#Build names db

# using real gi_tax: build up full taxonomy dict
gi_tax_path = FA.path(os.path.join(databaseUD, "gi_taxid_nucl/gi_taxid_nucl.dmp"))

# Parse BO and determine real hits

# Generate intervals along with hit profile into tax neighborhood



if RunControl['genarraytion']:
    ##Fix unordered gi to taxid
    #tIFN = FA.path([databaseUD,'gi_taxid_nucl.dmp'],fileType = 'tab')
    #inLst = tIFN.inDB()
    #inLst.sort(lambda x,y:cmp(int(x[0]), int(y[0])))
    #tIFN.outDB(inLst)
            
    txDB = T.BinarySearch(genarraytion_txid_path)
    BO_IFN = FA.path(['/Users/dominicsuciu/src/EclipseWorkspace/SCpyTools/src/test_dez/', 'BO_8_W13.txt'],fileType = 'tab')
    for line in BO_IFN:
        LLst = line.split('\t')
        if len(LLst) < 10: continue
        query_gi = LLst[0].split('|')[1]
        subj_gi = LLst[1].split('|')[1]
        q_tx = txDB[query_gi]
        s_tx = txDB[subj_gi]
        print 'TaxID: ', query_gi, q_tx, '-->', subj_gi, s_tx
        

    


if RunControl['Names'] | RunControl['Nodes']:
    #BuildNames Names Data Structure
    namesIFN = FA.path([taxdmp_databaseUD,'names.dmp'],fileType = 'tab')
    NamesDB = {} #taxID -> scientific name
    for line in namesIFN:
        LLst = line.split('\t')
        if len(LLst) < 6: continue
        if LLst[6] == 'scientific name':
            NamesDB[LLst[0]] = LLst[2]



if RunControl['Nodes']:
    import Tax.Taxon as Txn
    print str(Txn)
    nodesIFN = FA.path([taxdmp_databaseUD,'nodes.dmp'],fileType = 'tab')
    NodesDB = {}
    for line in nodesIFN:
        LLst = line.split('\t')##131573    |    294823    |    genus    |        |
        if len(LLst) < 5: continue
        NodesDB[LLst[0]] = [LLst[2], LLst[4]]
        
    
    IFN = genarraytion_txid_path #this file has all the taxID's I care about at this time
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
                print txID, "Not in NodesDB"
                break
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
   
   

    
              
#######################
#######  CODA   #######
#######################
    
raw_input(' It appears this module ran ok:')        



