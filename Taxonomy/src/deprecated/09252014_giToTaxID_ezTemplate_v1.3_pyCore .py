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
              'ExtractFastasByTaxid':1,
              
              
              'Nodes':0,#this makes neighborhood file that you edit then you run . . . .
              'MakeFinalDB':0, #this guy
              
              }


##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/




workingUD = "M:/bacteria_refseq/09252014_Run/"
outUD =os.path.join(workingUD, "Output/")
if not(os.path.exists(outUD)):
    os.mkdir(outUD)


#databaseUD = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy\\'

txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\'


if RunControl['giToTax']:
    
    tax_inpath = os.path.join(txDir,'gi_taxid_nucl.dmp')
    import Tax.GiToTax as T
    for tx in range(100):
        print T.bs[str(tx)]

    ifn = 'test-giList.txt'
    IFN = FA.path(['testdata/',ifn],fileType = 'fasta')
    #giToTax(tIFN, IFN, workingUD): Takes a file full of gi id's and creates a gi-tax_id file
    LF.giToTax(tax_inpath, IFN, outUD)
    


def CloseToTaxID(tax_a, taxLst, dist = 0):
    node_db.process_txidLst([tax_a])
    for b in taxLst:
        if node_db.TaxonDict[str(tax_a)].Distance(node_db.TaxonDict[str(b)]) == dist:
            return b
    return False
                                             
if RunControl['ExtractFastasByTaxid']:
    ##This version just makes a new file with gi to taxID From Fasta File
    tax_inpath = os.path.join(txDir,'gi_taxid_nucl.dmp')
    txDB = T.BinarySearch(tax_inpath)
    
    
    TakeTaxIDs = [139]
    TakeTaxIDs = [str(x) for x in  TakeTaxIDs] 
    
    import Tax.Nodes as Nodes
    node_db = Nodes.NodesDB(txDir)
    node_db.process_txidLst(TakeTaxIDs)
                            
                            
    tdict = {'notfound':[0,0]};saveSeqs = []
    for t in TakeTaxIDs:
        tdict[str(t)] = [0, 0]
        
    ifn = 'all_tickseqs.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    faDB = IFN.inDB()
    scrap = []
    for fa_obj in faDB:
        label, seq, length = fa_obj
        if not(label):
            print 'huh? ', fa_obj
            continue
        gi = label.split('|')[1]
        taxid = txDB[gi]
        node_db.process_txidLst([taxid])
        if taxid == 'notfound':
            scrap.append(fa_obj)
            continue
        
        mapped_taxid = CloseToTaxID(taxid, TakeTaxIDs, dist = 0)
        if mapped_taxid:
            taxid = str(taxid)
            tdict[mapped_taxid][0] +=1
            tdict[mapped_taxid][1] +=length
            nlabel = taxid + '_' + mapped_taxid +'-' + label
            saveSeqs.append([nlabel, seq, length])
    for t in tdict.keys():
        print t, tdict[t]
    OFN = FA.path([workingUD,'DEZ_bergirf_seqs.fasta'],fileType = 'fasta')
    OFN.outDB(saveSeqs)
    
    OFN = FA.path([workingUD,'scrap_seqs.fasta'],fileType = 'fasta')
    OFN.outDB(scrap)

    

            
if  RunControl['Nodes']:
    import Tax.Nodes as Nodes
    node_db = Nodes.NodesDB(txDir)
    
        
    ifn = 'gi-txid_all_safeIntervals.fasta' #this file has all the taxID's I care about at this time
    root = 'safeintervals'
    IFN = FA.path([outUD,ifn],fileType = 'tab')
    txid_Lst = IFN.inDB()
    node_db.process_txidLst([t[0] for t in txid_Lst])
    
    opath = os.path.join(outUD,'nrTaxNeighborhood_%s.txt'%root)
    node_db.export_TaxonDict(self, outpath=opath)
    
   
   
    

if RunControl['MakeFinalDB']:
    #This Module will take a Tax Neighborhood file that you curated by hand 
    # You took each taxID and decided wether you care about it or not
    #It will relabel your seqFile using this info
    ifn = '' #[tax    ToDo, other cols . . . . ]
    
    SI_IFN = FA.path([outUD ,'edited_nrTaxNeighborhood_%s.txt'%root],fileType = 'tab') #This is where you add the edited neighborhood file
    SIDict = SI_IFN.inDict(0,1)
    
    giTax_IFN = FA.path([outUD ,'gi-txid_all_safeIntervals.fasta'],fileType = 'tab')#add the mappings of gi to taxid
    
    GiToTaxDict = giTax_IFN.inDict(0,1)
    
    ifn = 'all_safeIntervals.fasta'  #'alldezseqs.fasta'
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    LF.FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD, SpecialInstructionDict = SIDict)
    
    
    
        



