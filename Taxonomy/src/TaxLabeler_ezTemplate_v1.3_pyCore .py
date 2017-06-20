#!/usr/bin/python
############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################

import os, fnmatch
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
              'prep_taxneigh':1,
              'faGiToTax':0,#Gi to taxid file from fasta
              }


##Make sure all the files especially gi_taxid_nucl are in taxdmp
##Get them from here:  ftp://ftp.ncbi.nih.gov/pub/taxonomy/

txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\taxdmp\\'


#databaseUD = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy\\'

if RunControl['prep_taxneigh']:
    #0 New form it updates TaxDB by itself
    import Tax.GiToTax as T
    #CUSTOMIZE
    
    gitax_path = os.path.join(txDir,'gi_taxid_nucl.dmp')
    gi_tax = T. BinarySearch(gitax_path)
    
        
    #1 Nodes allow distances to be computed between taxa
    import Tax.Nodes as Nodes
    TaxDB = Nodes.GiToTaxon(gi_tax, Nodes.NodesDB(txDir))
    #Nodes
    #taxneigh_path = 'M:/bacteria_refseq/11052014_Run/nrTaxNeighborhood_rickettsBart.txt'
    #node_db.import_TaxonDict(taxneigh_path)
    
    #dont need the whole bacterial world
    #taxneigh_path = 'M:/bacteria_refseq/genarraytion_taxid/nrTaxNeighborhood.txt'
    #TaxDB.import_TaxonDict(taxneigh_path)
    
    species_Dict = {}
    tDB = {} 
    #CUSTOMIZE
    taxneigh_path = None#'M:/bacteria_refseq/11052014_Run/edited_nrTaxNeighborhood_leptospira.txt'
    if taxneigh_path:
        print 'tax_neigh import from: ', taxneigh_path
        for line in fileinput.input(taxneigh_path):
            tobj = TaxObj(line)
            tDB[str(tobj.taxid)] = tobj
            species = tobj.species
            if not(species in species_Dict):
                species_Dict[species] = 0
            species_Dict[species] +=1
        species_index_Lst = species_Dict.keys()
        species_index_Lst.sort()
        for t in tDB.keys():
            tDB[t].RecordSpeciesIndex(species_index_Lst)
        TaxDB.TaxonDict.update(tDB)
   

def FASTA_giToTaxDict(GiToTaxDict, IFN, workingUD, SpecialInstructionDict = {}, Tax_Level='species'):
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    #This is modfified to use the taxon object that can generate a label at the species level
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
                tax_label = 'unknown'
                taxID = -1
            else:
                taxon = GiToTaxDict[int(gi)]
                if Tax_Level == 'no rank':
                    tax_label = taxon.no_rankLabel()
                elif Tax_Level == 'subspecies':
                    tax_label = taxon.subspeciesLabel()
                else:
                    tax_label = taxon.speciesLabel()
                taxID = taxon.species
            if SpecialInstructionDict:
                si = SpecialInstructionDict[taxID]
                if si =='dump': 
                    continue
                else: 
                    tax_label = si
                    taxID = si
                
            if random.random() > 0.99: print '%i: %s \t\t %s' %(c, label[:20], tax_label)
            wbuff = '>%s|%s' %(tax_label, sObj)
            newLabel = wbuff[:wbuff.find('\n')]
           
            LabelLst.append([c, gi, taxID, newLabel])
            outFP.write(wbuff)
            oDB.append([gi, taxID])
    
    OFN = FA.path([workingUD,'gi-txid_'+ifn],fileType = 'tab')
    OFN.outDB(oDB)
    print 'yo'+ workingUD+'txLabelLst_'+ifn
    OFN = FA.path([workingUD,'txLabelLst_'+ifn],fileType = 'tab')
    OFN.outDB(LabelLst)

if __name__ == '__main__':
    
    if len(sys.argv) == 1:
        cmd = """x.python M:   """
        fnLst = [
                 #"M:/bacteria_refseq/112914_Norovirus/norovirus_bgdb.fasta"
                #"M:/bacteria_refseq/11052014_Run/InCommonintervals_all_Leptospira-0.30.fasta"
                "M:/GGenomics/PaulSchaud_150917/39824_sequence.fasta"
        ]
    else:
        fnLst = sys.argv[1:]
        
    ofn = open('M:/GGenomics/PaulSchaud_150917/BGDB_allvir.fa', 'w')
    for root, dir, files in os.walk("M:/GGenomics/PaulSchaud_160110/all.fna/"):
        print 'root: ', root
        print ""
        for items in fnmatch.filter(files, "*"):
                print "..." + items
                if '.fna' in items:
                    b = open(os.path.join(root, items)).read()
                    ofn.write(b)
                    
        print ""
    ofn.close()
    
    fnLst = ['M:/GGenomics/PaulSchaud_150917/BGDB_allvir.fa']
    for inpath in fnLst:
        IFN= FA.path(inpath)
        workingUD = os.path.dirname(inpath)+'/'
        FASTA_giToTaxDict(TaxDB, IFN, workingUD, SpecialInstructionDict = {}, Tax_Level='species')
    

        

 
    
    
    
              
#######################
#######  CODA   #######
#######################
    
raw_input(' It appears this module ran ok:')        



