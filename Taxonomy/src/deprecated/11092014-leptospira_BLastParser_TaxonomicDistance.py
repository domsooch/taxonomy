############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################


##need to add application directory references
"""
This program takes a Blast In file as well as a Blast Out file and computes resultant sequence intervals in two distinct modes:
InCommon: looks for regions that are conserved. This is appropriate if you are trying to create genus-level probes

Upgrades:
1. Uses a dynamic Taxonomic Nodes database that expands based on usage
2. You have the option of setting taxon mappings by porviding your own tax neighborhood file
3. Uses compact BitOperations.BitArray() to keep track of in-common hits


"""

import os
import fileinput
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

#RUN Settings
DISTANCE_THRESH = 1 #genus  if s_tobj.Distance(q_tobj) <= DISTANCE_THRESH:
MIN_SEQ_SIZE = 100 #
HIT_THRESHOLD = 0.3 #if nh < hit_threshold:
RUN_LIMIT = None
MODE = 'InCommon'#MODE='Specific'
root = 'all_Leptospira-%0.2f'%HIT_THRESHOLD
inUD = 'M:\\bacteria_refseq\\11052014_Run'
outUD = inUD

fnLst = [['%i-BO_Lepto-vs-Lepto_W15.txt'%f,'%i-leptospira_genus.fasta'%f] for f in range(1,7,1)]




class Hit:
    def __init__(self, hitLst):
        #query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
        self.b = {'query':0,'subject':1,'%id':2,'align':3,'mm':4,'gap':5,'qstart':6,'qend':7,'sstart':8,'send':9,'e':10,'score':11}
        self.hit = hitLst
        self.query = self.hit[0]
        self.subject = self.hit[1]
        self.percid = float(self.hit[2])
        self.align = int(self.hit[3])
        self.mm = int(self.hit[4])
        self.gap = int(self.hit[5])
        self.qstart = int(self.hit[6])
        self.qend = int(self.hit[7])
        self.sstart = int(self.hit[8])
        self.send = int(self.hit[9])
        self.e  = float(self.hit[10])
        self.score = float(self.hit[11])
        if (self.qstart - self.qend)*(self.sstart-self.send) > 0:
            self.sense = 1
        else: self.sense = 0


class TaxObj:
    def __init__(self, line):
        line = line.replace('\n', '')
        self.Lst = line.split('\t')
        self.taxid = self.Lst.pop(0)
        self.species = None
        self.genus = None
        self.family = None
        self.order = None
        self.species_index = -1 #This is used for the In-Common analysis to keep track of part_species_Hit
        for x in self.Lst:
            #print x
            try:
                txid, name, rank = x.split('__')
            except:
                print 'Cant do anything with this: ', x
            if rank == 'species':
                self.species = txid
                continue
            if rank == 'genus':
                self.genus = txid
                continue
            if rank == 'family':
                self.family = txid
                continue
            if rank == 'order':
                self.order = txid
                continue
    def Distance(self, taxobj):
        if self.species == taxobj.species:
            return 0
        if self.genus == taxobj.genus:
            return 1
        if self.family == taxobj.family:
            return 2
        if self.order == taxobj.order:
            return 3
        return 10
    def RecordSpeciesIndex(self, species_index_Lst):
        self.species_index = species_index_Lst.index(self.species)

def TestSelfHit(subj, query, tax_dist_threshold):
    #Determines if a seq that is hit is the same or nearly the same taxon
    #Requires a dict from gi to tax_id
    #Requires a dictionary of class TaxObj keyed to the taxon
    if subj == query:
        return True
    s_gi = subj.split('|')[1]
    q_gi = query.split('|')[1]
    try:
        s_txid = gi_tax[s_gi].replace('\n','')
    except:
        print 'No taxid: s_gi', s_gi
        return False
    try:
        q_txid = gi_tax[q_gi].replace('\n','')
    except:
        print 'No taxid: q_gi', q_gi
        return False
    if s_txid == 'notfound':
        print 'notfound taxid for: subj', s_gi
        return False
    if q_txid == 'notfound':
        print 'notfound taxid for: query', q_gi
        return False
    if s_txid == q_txid:
        return True
    s_tobj = TaxDB[s_txid]
    q_tobj = TaxDB[q_txid]
    if s_tobj.Distance(q_tobj) <= tax_dist_threshold:
        return True
    return False
    

def processHits_Specific(faObj, hitLst, ofp, hit_threshold = 1, tax_dist_threshold = 1, verbose = False):
    takehittableintoconsideration_ZeroThresh = 0.0
    #Process hits (hitLst) for a given query sequence called faObj
    query, seq, slen = faObj
    print 'processHits_Specific: %i hits for query %s  %s'%(len(hitLst), query, hitLst[0].query)
    if query != hitLst[0].query:
        raw_input('ERR')    
    seqIndexLst = [0 for x in range(slen+1)]
    
    for h in hitLst:
        if h.subject == h.query:
            for i in range(h.qstart, h.qend, 1):
                seqIndexLst[i] += takehittableintoconsideration_ZeroThresh
            continue
        #Tests nearness of hit to query taxon Ignore hit if it is close to query taxon
        if TestSelfHit(h.subject, h.query, tax_dist_threshold):
            if verbose:print 'Self Hit: ', h.subject , h.query
            continue
        #Test if start and ends are compatible
        if h.qstart > slen or h.qend > slen:
            if verbose:print 'BAD HIT: Writin %i to %i on seq of len %i  %s %s'%(h.qstart, h.qend, slen, query, str(h.hit))
            continue
        #print 'Writing hit onto hitLst: %i to %i'%(h.qstart, h.qend)
        #Register this hit al;ong numberline representing sequence
        for i in range(h.qstart, h.qend, 1):
            seqIndexLst[i] +=1
    
    oLst = []
    min_seq_size = MIN_SEQ_SIZE
    seq_saved = 0
    start = None
    sLst = []
    num_seqs = 0
    #print seqIndexLst[:1000]
    #Extract clean seq intervals
    for s in range(len(seqIndexLst)):
        nh = seqIndexLst[s]
        if nh < hit_threshold and nh >= takehittableintoconsideration_ZeroThresh: #Need regions that hit self
            if start != None: 
                if s%50000 == 0: 
                    if verbose:print nh ,
                continue
            else:
                if verbose:print '\nstart interval at %i'%s,
                start = s
        else: #This indicates you have hit a non-unique region you need to now start collecting the interval
            if start !=None:
                seq_len = s - start -1
                if  seq_len > min_seq_size:
                    num_seqs +=1
                    seq_saved += seq_len
                    save_label = ">%s|HT%i-TT%i-ZT%.1f|%i-%i|"%(query, hit_threshold, tax_dist_threshold, takehittableintoconsideration_ZeroThresh,  start, s)
                    save_seq = seq[start:s]+'\n'
                    if save_seq:
                        sLst.append(save_label)
                        sLst.append(save_seq)
                    if verbose:print 'saved at %i  len=%i'%(s, seq_len)
                    start = None
                else:
                    if verbose:print 'Stop'
                    start = None
            else:
                #print nh
                continue
    #Mop-Up the last interval
    if start != None:
        seq_len = s - start -1
        if  seq_len > min_seq_size:
            num_seqs +=1
            seq_saved += seq_len
            save_label = ">%s|HT%i-TT%i-ZT%.1f|%i-%i|"%(query, hit_threshold, tax_dist_threshold, takehittableintoconsideration_ZeroThresh,  start, s)
            save_seq = seq[start:s]+'\n'
            if save_seq:
                sLst.append(save_label)
                sLst.append(save_seq)
            if verbose:print 'saved %i'%seq_len
    ofp.write('\n'.join(sLst))
    print '\n', [query, num_seqs, seq_saved], '\n\n\n'
    return [query, num_seqs, seq_saved]

def ProcessBLASTRun_Specific(f, faObj, procLst, ofp, LL=None):
    if LL:
        return ['f', 'taxid', 'query', 'hits_processed', 'num_seqs', 'seq_saved', 'orig_seqLen', 'hit_threshold', 'tax_dist_threshold', 'MIN_SEQ_SIZE', 'species', '==> taxonomy']
    #DISTANCE_THRESH = 1 #genus  if s_tobj.Distance(q_tobj) <= DISTANCE_THRESH:
    #HIT_THRESHOLD = 1 #if nh < hit_threshold:
    tax_dist_threshold = DISTANCE_THRESH
    hit_threshold = HIT_THRESHOLD
    
    while True:
        #processHits(faObj, hitLst, ofp, hit_threshold = 1, tax_dist_threshold = 1, verbose = False):
        query, num_seqs, seq_saved = processHits_Specific(faObj, procLst, ofp, hit_threshold, tax_dist_threshold)
        gi = query.split('|')[1]
        species = 'notfound'
        try:
            taxid = gi_tax[gi].replace('\n','')
            if taxid in TaxDB:
                tobj = TaxDB[taxid]
                species = tobj.species
                specLst = tobj.Lst
            else:
                specLst = ['nospec']
        except:
            taxid = 'notfound'
            specLst = ['nospec']
        results = [f, taxid, query, len(procLst), num_seqs, seq_saved, faObj[2], hit_threshold, tax_dist_threshold, MIN_SEQ_SIZE, '==>', species] + specLst
        if num_seqs > 0:
            return results
        else:
            if tax_dist_threshold >1:
                tax_dist_threshold -=1
                hit_threshold +=1
            else:
                tax_dist_threshold +=1
            if hit_threshold > 4:
                print 'Failed for: ', query 
                return results
            print 'Trying again query: %s with hitThreshold: %i TaxDistThreshold: %i' %(query, hit_threshold, tax_dist_threshold)


def processHits_InCommon(faObj, hitLst, ofp, hit_min_threshold = 0.7, NumSpeciesToTrac=1, verbose = False):
    #Process hits (hitLst) for a given query sequence called faObj
    #This version looks for regions that are in common within a genus Exactly the opposite of the previous analysis
    query, seq, slen = faObj
    print 'processHits_InCommon: %i hits for query %s  %s'%(len(hitLst), query, hitLst[0].query)
    if query != hitLst[0].query:
        raw_input('ERR')    
    seqIndexLst = [BitOperations.BitArray() for x in range(slen+1)]
    
    query_gi = query.split('|')[1]
    query_txid = gi_tax[query_gi]
    #print query, query_txid, TaxDB
    query_species_index = TaxDB[str(query_txid)].species_index

    for h in hitLst:
        if h.subject == h.query:
            species_index = query_species_index
            for i in range(h.qstart, h.qend, 1):
                seqIndexLst[i].setBit(species_index)
            continue
        #Find species_index of this hit
        subj_gi = h.subject.split('|')[1]
        subj_txid = gi_tax[subj_gi]
        species_index = TaxDB[str(subj_txid)].species_index
        #Test if start and ends are compatible
        if h.qstart > slen or h.qend > slen:
            if verbose:print 'BAD HIT: Writin %i to %i on seq of len %i  %s %s'%(h.qstart, h.qend, slen, query, str(h.hit))
            continue
        #print 'Writing hit onto hitLst: %i to %i'%(h.qstart, h.qend)
        #Register this hit al;ong numberline representing sequence
        for i in range(h.qstart, h.qend, 1):
            seqIndexLst[i].setBit(species_index)
    
    for x in range(slen+1):
        score = float(len(seqIndexLst[x]))/NumSpeciesToTrack
        seqIndexLst[x] = score
        
    oLst = []
    min_seq_size = MIN_SEQ_SIZE
    seq_saved = 0
    start = None
    sLst = []
    num_seqs = 0
    #print seqIndexLst[:1000]
    #Extract clean seq intervals
    for s in range(len(seqIndexLst)):
        nh = seqIndexLst[s]
        if nh > hit_min_threshold: #Need regions that hit self
            if start != None: 
                if s%50000 == 0: 
                    if verbose:print nh ,
                continue
            else:
                if verbose:print '\nstart interval at %i'%s,
                start = s
        else: #This indicates you have hit a non-unique region you need to now start collecting the interval
            if start !=None:
                seq_len = s - start -1
                if  seq_len > min_seq_size:
                    num_seqs +=1
                    seq_saved += seq_len
                    ss= sum(seqIndexLst[start:s]); avgHitScore = float(ss)/seq_len
                    save_label = ">%s|HT%0.2f-HS%0.2f|%i-%i|"%(query, hit_min_threshold, avgHitScore,  start, s)
                    save_seq = seq[start:s]+'\n'
                    if save_seq:
                        sLst.append(save_label)
                        sLst.append(save_seq)
                    if verbose:print 'saved at %i  len=%i'%(s, seq_len)
                    start = None
                else:
                    if verbose:print 'Stop'
                    start = None
            else:
                #print nh
                continue
    #Mop-Up the last interval
    if start != None:
        seq_len = s - start -1
        if  seq_len > min_seq_size:
            num_seqs +=1
            seq_saved += seq_len
            ss= sum(seqIndexLst[start:s]);avgHitScore = float(ss)/seq_len
            save_label = ">%s|HT%0.2f-HS%0.2f|%i-%i|"%(query, hit_threshold, avgHitScore,  start, s)
            save_seq = seq[start:s]+'\n'
            if save_seq:
                sLst.append(save_label)
                sLst.append(save_seq)
            if verbose:print 'saved %i'%seq_len
    ofp.write('\n'.join(sLst))
    print '\n', [query, num_seqs, seq_saved], '\n\n\n'
    return [query, num_seqs, seq_saved]

def ProcessBLASTRun_InCommon(f, faObj, procLst, ofp, NumSpeciesToTrack = 0, LL=None):
    #f, faObj, procLst, ofp, NumSpeciesToTrack, LL=None
    if LL:
        return ['f', 'taxid', 'query', 'hits_processed', 'num_seqs', 'seq_saved', 'orig_seqLen', 'hit_threshold', 'tax_dist_threshold', 'MIN_SEQ_SIZE', 'species', '==> taxonomy']
    #DISTANCE_THRESH = 1 #genus  if s_tobj.Distance(q_tobj) <= DISTANCE_THRESH:
    #HIT_THRESHOLD = 1 #if nh < hit_threshold:
    tax_dist_threshold =-1
    hit_min_threshold = HIT_THRESHOLD
    
    #processHits(faObj, hitLst, ofp, hit_threshold = 1, tax_dist_threshold = 1, verbose = False):
    query, num_seqs, seq_saved = processHits_InCommon(faObj, procLst, ofp, hit_min_threshold, NumSpeciesToTrack)
    gi = query.split('|')[1]
    
    try:
        taxid = str(gi_tax[gi])
        tobj = TaxDB[taxid]
        species = tobj.species
        specLst = tobj.Lst
    except:
        taxid = 'notfound'
        species = 'notfound'
        specLst = ['nospec']
    results = [f, taxid, query, len(procLst), num_seqs, seq_saved, faObj[2], hit_min_threshold, tax_dist_threshold, MIN_SEQ_SIZE, '==>', species] + specLst
    return results
        

        

if RunControl['prep_taxneigh']:
    #0 New form it updates TaxDB by itself
    import Tax.GiToTax as T
    #CUSTOMIZE
    txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\'
    gitax_path = os.path.join(txDir,'gi_taxid_nucl.dmp')
    gi_tax = T. BinarySearch(gitax_path)
    
        
    #1 Nodes allow distances to be computed between taxa
    import Tax.Nodes as Nodes
    TaxDB = Nodes.NodesDB(txDir)
    #Nodes
    #taxneigh_path = 'M:/bacteria_refseq/11052014_Run/nrTaxNeighborhood_rickettsBart.txt'
    #node_db.import_TaxonDict(taxneigh_path)
    
    #dont need the whole bacterial world
    #taxneigh_path = 'M:/bacteria_refseq/genarraytion_taxid/nrTaxNeighborhood.txt'
    #TaxDB.import_TaxonDict(taxneigh_path)
    
    species_Dict = {}
    tDB = {} 
    #CUSTOMIZE
    taxneigh_path = 'M:/bacteria_refseq/11052014_Run/edited_nrTaxNeighborhood_leptospira.txt'
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
    NumSpeciesToTrack = len(species_index_Lst)
    print 'NumSpeciesToTrack: %i'%NumSpeciesToTrack
    for s in range(len(species_index_Lst)):
        print species_index_Lst[s], s
        
            
        
    

if RunControl['ProcessBlastRuns']:
    seq_out_path = os.path.join(outUD, '%sintervals_%s.fasta'%(MODE, root))
    ofp = open(seq_out_path, 'w')
    #hitLst = []
    out_results_path = os.path.join(outUD, 'results-%s.txt'%root)
    outresults_fp = open(out_results_path, 'w')
    LL = ProcessBLASTRun_Specific('dummy', [], [], 'dummy', LL=1)#Labels
    outresults_fp.write('\t'.join(LL)+'\n')
              
    
    
    for blast_fn_set in fnLst:
        i = 0
        bo_fn , bi_fn = blast_fn_set
        
        bo_inpath = os.path.join(inUD, bo_fn)
        bi_inpath = os.path.join(inUD, bi_fn)
        infp = open(bi_inpath, 'r')
        b = infp.read()
        
        infp.close()
        fbuffLst = b.split('>')
        
        #Build FADict
        FADict = {};results = []
        for fobj in fbuffLst:
            fobj = fobj.split('\n')
            label = fobj.pop(0).split(' ')[0]
            seq = ''.join(fobj)
            FADict[label] = [label, seq, len(seq)]
            #print label, len(seq), seq[:100]
        infp = 1
        query_num = 0
        n = 0;q=0
        current_query = None
        #OutputFile
        print bo_inpath
        BO_INP = fileinput.input(bo_inpath)
        for line in BO_INP:
            #print line
            n+=1
            if n%1000 == 0:
                print 'fn: %s Processed %i hits %i queries saved %i result seqs'%(bo_inpath, n, query_num, len(results))
            try:
                h = Hit(line.split('\t'))
            except:
                #print 'Nothing to do with this: %s'%line
                continue
            if current_query == None:
                current_query = h.query
                try:
                    faObj = FADict[h.query]
                except:
                    current_query = None
                    continue
                procLst = []
            if current_query == h.query:
                procLst.append(h)
            else:
                if MODE == 'Specific':
                    results = ProcessBLASTRun_Specific(f, faObj, procLst, ofp, LL=None)
                elif MODE == 'InCommon':
                    NumSpeciesToTrack = len(species_index_Lst)
                    results = ProcessBLASTRun_InCommon(f, faObj, procLst, ofp, NumSpeciesToTrack, LL=None)
                if results: 
                    results = [str(r) for r in results]
                    outresults_fp.write('\t'.join(results)+'\n')
                #print 'New accumulation: ', h.query
                current_query = h.query
                faObj = FADict[h.query]
                procLst = [h]
                query_num +=1
            if RUN_LIMIT and query_num>RUN_LIMIT: break
        BO_INP.close()
        if procLst:
            if MODE == 'Specific':
                results = ProcessBLASTRun_Specific(f, faObj, procLst, ofp, LL=None)
            elif MODE == 'InCommon':
                NumSpeciesToTrack = len(species_index_Lst)
                results = ProcessBLASTRun_InCommon(f, faObj, procLst, ofp, NumSpeciesToTrack, LL=None)
            if results: 
                results = [str(r) for r in results]
                outresults_fp.write('\t'.join(results)+'\n')
        
    ofp.close()

    outresults_fp.close()
                         
        
    








                
#######################
#######  CODA   #######
#######################
    
#raw_input(' It appears this module ran ok:')        
sys.stdout = logObj.uncouple()

#######  EOM   #######
