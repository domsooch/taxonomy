############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################


##need to add application directory references
"""
Looks at all seqs and looks for missing gi's


"""

import os
import fileinput, random
import sys
##Location
WORK = {'subRoutineDir':'c:/work/Python/New_Core/', #
        'WorkDir' : os.getcwd()+'/'}
CDrive = {'subRoutineDir':'C:\\Users\\dominic\\EclipseWorkspace\\new_pyCore\\src\\',
        'WorkDir' :os.getcwd()+'/',
        }
location = CDrive

subRoutineDir = location['subRoutineDir']

sys.path.append(subRoutineDir[:-1])
sys.path.append('/Users/dominicsuciu/Dropbox/genarraytion/new_pyCore/src/')
#Import Modules
import pyCore
from pyCore.pyLib.__init__ import ImportModulesExecString
execList = ImportModulesExecString()
for l in execList:
    print l
    exec l
from Util import *

#location = LOC.LocationObj(locationdict =location)
#location.propagateLocation(modLst)
inUD = os.getcwd()
UD = inUD
outUD = inUD
logUD = inUD
appDir = inUD



#Logging
errFP = open(logUD + 'log.txt','w')
screenWriter = sys.stdout
logObj = US.logger(errFP,screenWriter,Verbose =1)
sys.stdout = logObj #this diverts the print statement's stream

##ERROR Logging
#fsock = open(logUD + '_ezTemplate_error.log', 'w')
#sys.stderr = fsock

#############################     New ezTemplate-1.0  091509  #################################
#############################################################################################
#############################     end_EZTEMPLATE HEADER   ###################################
#############################################################################################


verbose = 1

###RUN CONTROL

RunControl = {'prep_taxneigh':1,
              'ProcessBlastRuns':1,
              }


inUD = 'M:\\bacteria_refseq\\11052014_Run'
outUD = inUD






def FASTA_giToTaxDict(GiToTaxDict, inpath, UsedTaxDict, SaveSeq=False):
    ##This makes new file with taxID as part of label But it Uses a pre-determined giToTaxID dictionary
    #This is modfified to use the taxon object that can generate a label at the species level
    inFP = open(inpath, 'r')
    c = 0
    next_inbuff = ''
    while 1:
        inbuff = next_inbuff + inFP.read(10000000)
        next_inbuff = ''
        while True:
            line = inFP.readline()
            if not(line): break
            if '>' in line:
                next_inbuff = line
                break
            inbuff = inbuff + line
        print '%i records Read at : %s' %(len(UsedTaxDict), str(inFP.tell()))
        if not(inbuff):
            inFP.close()
            break
        inLst = inbuff.split('>')
        for sObj in inLst:
            c+=1
            if not(sObj):continue
            label_pos = sObj.find('\n')
            seq = sObj[label_pos:].replace('\n','')
            label = sObj[:label_pos]
            seq_len = len(seq)
            #FIX Label
            if not('|' in label):
                continue
            gi_str = label.split('|')[1]
            if '-' in gi_str and '_' in gi_str:
                gi = gi_str.split('-')[-1]
                label = label.replace(gi_str, gi)
            #/FIXLabel
            gi = label.split('|')[1]
            if not(gi):
                tax_label = 'unknown'
                taxID = -1
                raw_input('Not TaxID for this label:', label)
            else:
                print label
                taxon = GiToTaxDict[gi]
                if taxon.species_name == 'uk':
                    NotFoundLst.append(label)
                tax_label = taxon.speciesLabel()
                taxID = taxon.species
            if 'genus' in label[:10]:
                newlabel = label
            elif 'unknown' in label[:20]:
                newlabel = label
            else:
                newlabel = "%s|%s"%(tax_label, label[3:])
            faObj = [newlabel, seq, seq_len]
            if not(tax_label in UsedTaxDict):
                UsedTaxDict[tax_label] =[0, 0, []]
            UsedTaxDict[tax_label][0] +=1
            UsedTaxDict[tax_label][1] +=seq_len
            if SaveSeq:
                UsedTaxDict[tax_label][2].append(faObj)
            if random.random() > 0.999: print '%i: %s \t\t %s' %(c, label[:20], tax_label)
           
def DownSample(FALst, MaxSeq_space):
    outFALst = []
    numSeqs = len(FALst)
    seqSpace = 0
    for s in FALst:
        seqSpace += s[2]
    if seqSpace < MaxSeq_space:
        print 'DownSample has nothing to do here: %s'%FALst[0][0]
        for f in range(len(FALst)):
            label = FALst[f][0]
            FALst[f][0] = label.split(' ')[0].replace(',','')
        return FALst
    print 'DownSample %i seqs for %i seq_space to MaxSeq_space: %i for: %s'%(numSeqs, seqSpace, MaxSeq_space, FALst[0][0])
    random.shuffle(FALst)
    seqSpaceTaken = 0
    for faobj in FALst:
        label, seq, sl = faobj
        label = label.split(' ')[0]
        label = label.replace(',','')
        seqSpaceTaken += sl
        outFALst.append([label, seq, sl])
        if seqSpaceTaken >MaxSeq_space:
            break
    return outFALst

if RunControl['prep_taxneigh']:
    #0 New form it updates TaxDB by itself
    import Tax.GiToTax as T
    #CUSTOMIZE
    txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\'
    gitax_path = os.path.join(txDir,'gi_taxid_nucl.dmp')
    gi_tax = T.BinarySearch(gitax_path)
    
        
    #1 Nodes allow distances to be computed between taxa
    import Tax.Nodes as Nodes
    GiToTaxon = Nodes.GiToTaxon(gi_tax, Nodes.NodesDB(txDir))
    #Nodes
    #taxneigh_path = 'M:/bacteria_refseq/11052014_Run/nrTaxNeighborhood_rickettsBart.txt'
    #node_db.import_TaxonDict(taxneigh_path)
    
    #dont need the whole bacterial world
    #taxneigh_path = 'M:/bacteria_refseq/genarraytion_taxid/nrTaxNeighborhood.txt'
    #TaxDB.import_TaxonDict(taxneigh_path)



a = ['M:\\bacteria_refseq\\11052014_Run\\DEZ_leptospira_genus.fasta',
         'M:\\bacteria_refseq\\11052014_Run\\DEZ_rickettsBart.fasta',
         'M:\\bacteria_refseq\\09252014_Run\\DEZ_Borellia.fasta',
         'M:\\bacteria_refseq\\09252014_Run\\all_tickseqs.fasta',]

fnLst = [
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\DEZcombined_.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\babesiacanis_5867.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\rickettsiaHelvetica_35789.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\babesia_microti_5868.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\borellia_afzelli_29518.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\borellia_miyamotoi_47466.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\anaplasmaphagocty_948.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\InCommonintervals_anaplasigenusInCommon-0.40.fasta',
         
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\erlich_1423420.fasta',
         'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\deartickvirus_58535.fasta',
         
         ]
indir = 'M:\\bacteria_refseq\\0925and1105_ProbeDesigns\\'

if True:
    #Survey
    UsedTaxDict = {};NotFoundLst= []
    for inpath in fnLst:
        print inpath
        FASTA_giToTaxDict(GiToTaxon, inpath, UsedTaxDict, SaveSeq=True)
    oLst = []
    faLst = []
    for it in UsedTaxDict.items():
        spec_label = it[0]
        num = it[1][0]
        seq_space =  it[1][1]
        spec_FaLst = it[1][2]
        faLst.extend(DownSample(spec_FaLst, 300000))
        #genusTxID[773]_Bartonella_alsatica[52764]|395409949|gb|AIME00000000.1|AIME01000000 Bartonella alsatica IBS 382
        if '_' in spec_label:
            p = spec_label.find('_')
            genus_num = spec_label[:p].replace('genusTxID[', '').replace(']','')
            species = spec_label[p:]
            species_num = species.split('[')[1].replace(']', '')
            l = [spec_label, genus_num, species, species_num, str(num), str(seq_space)]
        oLst.append('\t'.join(l))
    #oLst = ['\t'.join([str(x) for x in list(it)]) for it in UsedTaxDict.items()]
    ofp = open('used_taxIDs.txt', 'w')
    ofp.write('\n'.join(oLst))
    ofp.write('\n'.join(NotFoundLst))
    ofp.close()
    
    fLst = oLst = [','.join(['index',    'seq_accession',    'seq_custom',    'targetname',    'probes',    'replicates',    'notes'])]
    i= 1
    seqSpace = 0
    for f in range(len(faLst)):
        faObj = faLst[f]
        l = [i, '', faObj[1], faObj[0], min(200, faObj[2]/100), 1, '' ]
        seqSpace += faObj[2]
        i+=1
        fLst.append(','.join([str(c) for c in l]))
    IT_path = os.path.join(indir, '1114_IncludedTargets.csv')
    ofp = open(IT_path, 'w')
    ofp.write('\n'.join(fLst))
    ofp.close()
    print 'TotalSeqSpace: %i'%seqSpace
    
    
    
 
            
        
    



    








                
#######################
#######  CODA   #######
#######################
    
#raw_input(' It appears this module ran ok:')        
sys.stdout = logObj.uncouple()

#######  EOM   #######
