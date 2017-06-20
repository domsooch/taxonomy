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
              'ExtractFastasByTaxid':0,
              'Nodes':1,#this makes neighborhood file that you edit then you run . . . .
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
    
    #0 New form it updates TaxDB by itself
    import Tax.GiToTax as T
    txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\'
    gitax_path = os.path.join(txDir,'gi_taxid_nucl.dmp')
    gi_tax = T. BinarySearch(gitax_path)
    #giToTax(tIFN, IFN, workingUD): Takes a file full of gi id's and creates a gi-tax_id file
    LF.giToTax(tax_inpath, IFN, outUD)
    


def CloseToTaxID(tax_a, taxLst,  DistanceToAccept = 0):
    node_db.process_txidLst([tax_a])
    for b in taxLst:
        d = node_db.TaxonDict[str(tax_a)].Distance(node_db.TaxonDict[str(b)])
        #print '%s to %s dist = %s'%(tax_a, b, d)
        if  DistanceToAccept == d:
            return b
    return False
                                             
if RunControl['ExtractFastasByTaxid']:
    import Tax.Nodes as Nodes
    ##This version just makes a new file with gi to taxID From Fasta File
    tax_inpath = os.path.join(txDir,'gi_taxid_nucl.dmp')#0. Add location of taxdmp
    txDB = T.BinarySearch(tax_inpath)
    
    ifn = 'DEZ_Borellia.fasta'#1. Fasta Input File
    DistanceToAccept = 2#2. set distance to accept 0=species, 1 = genus, 2 = family, 3 = order
    ofn = 'DEZ_%s'%ifn
    tx_ofn = 'txidLst_%s'%ifn
    
    TakeTaxIDs = [138]#2. Add taxid's you want from input file
    TakeTaxIDs = [str(x) for x in  TakeTaxIDs] 
    
    node_db = Nodes.NodesDB(txDir)
    node_db.process_txidLst(TakeTaxIDs)
                            
                            
    tdict = {'notfound':[0,0]};saveSeqs = []
    for t in TakeTaxIDs:
        tdict[str(t)] = [0, 0]
        
    
    IFN = FA.path([workingUD,ifn],fileType = 'fasta')
    faDB = IFN.inDB()
    scrap = []
    tax_idLst = []
    for fa_obj in faDB:
        label, seq, length = fa_obj
        if len(seq) < 10:
            print 'NoSeq for %s'%label
        if not(label):
            print 'huh? ', fa_obj
            continue
        gi = label.split('|')[1]
        taxid = txDB[gi]
        node_db.process_txidLst([taxid])
        if taxid == 'notfound':
            scrap.append(fa_obj)
            continue
        if not(taxid in tax_idLst):
            tax_idLst.append(taxid)
        mapped_taxid = CloseToTaxID(taxid, TakeTaxIDs,  DistanceToAccept =  DistanceToAccept)
        if mapped_taxid:
            taxid = str(taxid)
            tdict[mapped_taxid][0] +=1
            tdict[mapped_taxid][1] +=length
            nlabel = taxid + '_' + mapped_taxid +'-' + label
            saveSeqs.append([nlabel, seq, length])
    for t in tdict.keys():
        print t, tdict[t]
    OFN = FA.path([workingUD,ofn],fileType = 'fasta')
    OFN.outDB(saveSeqs)
    
    OFN = FA.path([workingUD,'scrap_%s'%ofn],fileType = 'fasta')
    OFN.outDB(scrap)
    
    ofp = open(os.path.join(workingUD, tx_ofn), 'w')
    ofp.write('\n'.join([str(t) for t in tax_idLst]))
    ofp.close()

    

            
if  RunControl['Nodes']:
    import Tax.Nodes as Nodes
    node_db = Nodes.NodesDB(txDir)
    
    tx_ofn = 'txidLst_DEZ_Borellia.fasta'#1. this file has all the taxID's I care about at this time
    root = 'borellia'
    ifn = tx_ofn 
    IFN = FA.path([workingUD,ifn],fileType = 'tab')
    txid_Lst = IFN.inDB()
    #print txid_Lst
    node_db.process_txidLst([t[0] for t in txid_Lst])
    
    opath = os.path.join(workingUD,'nrTaxNeighborhood_%s.txt'%root)
    node_db.export_TaxonDict(outpath=opath)
    
   
   
    

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
    
    
    
        



