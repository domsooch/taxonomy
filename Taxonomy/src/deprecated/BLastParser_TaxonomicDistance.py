############################################################################################
############################# EZTEMPLATE HEADER ############################################
############################################################################################


##need to add application directory references
# This program Blasts Sequences against a bgfile and eliminates regions that hit other regions from taxonomically distant organisms
#Components:
#1. class TaxObj:
#2. def TestSelfHit(subj, query, tax_dist_threshold)
#
#Run: ProcessBLASTRun(f, faObj, procLst, ofp, LL=None) 

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
              'mod1':1,
              }

DISTANCE_THRESH = 1 #genus  if s_tobj.Distance(q_tobj) <= DISTANCE_THRESH:
MIN_SEQ_SIZE = 100 #
HIT_THRESHOLD = 1 #if nh < hit_threshold:


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
    

def processHits(faObj, hitLst, ofp, hit_threshold = 1, tax_dist_threshold = 1, verbose = False):
    #Process hits (hitLst) for a given query sequence called faObj
    query, seq, slen = faObj
    print 'processHits: %i hits for query %s  %s'%(len(hitLst), query, hitLst[0].query)
    if query != hitLst[0].query:
        raw_input('ERR')    
    seqIndexLst = [0 for x in range(slen+1)]
    takehittableintoconsideration_ZeroThresh = 0.0
    for h in hitLst:
        if h.subject == h.query:
            takehittableintoconsideration_ZeroThresh = 0.5
            for i in range(h.qstart, h.qend, 1):
                seqIndexLst[i] +=0.5
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



def ProcessBLASTRun(f, faObj, procLst, ofp, LL=None):
    if LL:
        return ['f', 'taxid', 'query', 'hits_processed', 'num_seqs', 'seq_saved', 'orig_seqLen', 'hit_threshold', 'tax_dist_threshold', 'MIN_SEQ_SIZE', 'species', '==> taxonomy']
    #DISTANCE_THRESH = 1 #genus  if s_tobj.Distance(q_tobj) <= DISTANCE_THRESH:
    #HIT_THRESHOLD = 1 #if nh < hit_threshold:
    tax_dist_threshold = DISTANCE_THRESH
    hit_threshold = HIT_THRESHOLD
    
    while True:
        #processHits(faObj, hitLst, ofp, hit_threshold = 1, tax_dist_threshold = 1, verbose = False):
        query, num_seqs, seq_saved = processHits(faObj, procLst, ofp, hit_threshold, tax_dist_threshold)
        gi = query.split('|')[1]
        species = 'notfound'
        if gi in gi_tax:
            taxid = gi_tax[gi].replace('\n','')
            if taxid in TaxDB:
                tobj = TaxDB[taxid]
                species = tobj.species
                specLst = tobj.Lst
            else:
                specLst = ['nospec']
        else:
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


class TaxObj:
    def __init__(self, line):
        line = line.replace('\n', '')
        self.Lst = line.split('\t')
        self.taxid = self.Lst.pop(0)
        self.species = None
        self.genus = None
        self.family = None
        self.order = None
        for x in self.Lst:
            #print x
            txid, name, rank = x.split('__')
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
        

if RunControl['prep_taxneigh']:
    
    gi_tax = {}
    gitax_path = "M:/bacteria_refseq/genarraytion_taxid/genarraytion_gi_taxid_nucl.dmp"
    print 'gi_tax import from: ', gitax_path
    for line in fileinput.input(gitax_path):
        gi, txid = line.split('\t')
        gi_tax[gi] = txid.replace('\n', '')
        

    TaxDB = {} 
    taxneigh_path = 'genarraytion_taxid/nrTaxNeighborhood.txt'
    print 'tax_neigh import from: ', taxneigh_path
    for line in fileinput.input(taxneigh_path):
        tobj = TaxObj(line)
        TaxDB[tobj.taxid] = tobj
        
            
        
    

if RunControl['mod1']:

    #inUD = 'M:\\bacteria_refseq\\'
    hitLst = []
    oLst = [ProcessBLASTRun('dummy', [], [], 'dummy', LL=1)]#Labels
    for f in range(11):
        print f
        bo_fn ='BO_%i_W13.txt'%f
        bi_fn = '%i-BI.fasta'%f
        
        bo_inpath = os.path.join(inUD, bo_fn)
        bi_inpath = os.path.join(inUD, bi_fn)
        infp = open(bi_inpath, 'r')
        b = infp.read()
        
        infp.close()
        fbuffLst = b.split('>')
        
        #Build FADict
        FADict = {}
        for fobj in fbuffLst:
            fobj = fobj.split('\n')
            label = fobj.pop(0).split(' ')[0]
            seq = ''.join(fobj)
            FADict[label] = [label, seq, len(seq)]
            #print label, len(seq), seq[:100]
       
        infp = 1
        query_num = 0
        n = 0
        current_query = None
        
        #OutputFile
        ofp = open('safeintervals_dist-%i_%s'%(DISTANCE_THRESH,bi_fn), 'w')
        for line in fileinput.input(bo_inpath):
            #print line
            n+=1
            #if n%1000 == 0: break
            try:
                h = Hit(line.split('\t'))
                #hitLst.append(h)
            except:
                #print 'Nothing to do with this: %s'%line
                continue
            if current_query == None:
                current_query = h.query
                faObj = FADict[h.query]
                procLst = []
            if current_query == h.query:
                procLst.append(h)
            else:
                results = ProcessBLASTRun(f, faObj, procLst, ofp, LL=None)
                oLst.append(results)
                current_query = h.query
                faObj = FADict[h.query]
                procLst = [h]
                query_num +=1
            #if query_num>2: break
            
        if procLst:
            results = ProcessBLASTRun(f, faObj, procLst, ofp, LL=None)
            oLst.append(results)
        
        ofp.close()
        ofp = open('results-dist%i.txt'%DISTANCE_THRESH, 'w')
        buff = '\n'.join(['\t'.join([str(x) for x in obj]) for obj in oLst])
        ofp.write(buff)
        ofp.close()
                         
        
    








                
#######################
#######  CODA   #######
#######################
    
#raw_input(' It appears this module ran ok:')        
sys.stdout = logObj.uncouple()

#######  EOM   #######
