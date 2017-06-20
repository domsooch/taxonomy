
import os, sys, fileinput

import Taxon as Txn
    

class GiToTaxon:
    def __init__(self, giToTaxID, Nodesdb):
        self.giToTaxID = giToTaxID
        self.Nodes = Nodesdb
    def __getitem__(self, gi):
        taxID = self.giToTaxID[gi]
        return self.Nodes[taxID]
        

            
    
    
class NodesDB:
    def __init__(self, taxonomy_dbDir):
        self.taxonomy_dbDir = taxonomy_dbDir
        self.names_path = os.path.join(self.taxonomy_dbDir, 'names.dmp')
        self.nodes_path = os.path.join(self.taxonomy_dbDir, 'nodes.dmp')
        self.NodesDB = {}
        self.NamesDB = {}
        self.TaxonDict = {}
        self.process_Nodes()
        self.process_Names()
    def __getitem__(self, txID):
        t = txID
        if not(t in self.TaxonDict):
            self.TaxonDict[t] = self.process_txid(t)
        return self.TaxonDict[t]
    def process_Nodes(self):
        i = 0
        for line in fileinput.input(self.nodes_path):
            i+=1
            if i%100000 == 0: print 'process_Nodes: %i'%i
            LLst = line.split('\t')##131573    |    294823    |    genus    |        |
            #if len(LLst) < 5: continue
            self.NodesDB[LLst[0].strip()] = [LLst[2], LLst[4]]
    def process_Names(self):
        i = 0
        for line in fileinput.input(self.names_path):
            i+=1
            LLst = line.split('\t')
            if len(LLst) < 6: continue
            if LLst[6] == 'scientific name':
                self.NamesDB[LLst[0]] = LLst[2]
            if i%100000 == 0:
                if len(LLst) >2:
                    n = LLst[2]
                else:
                    n = 'noname'
                print 'process_Names: %i  %s'%(i, n)
    def process_txidLst(self, txidLst):
        for txID in txidLst:
            txID = str(txID)
            if txID in self.TaxonDict: continue
            t = Txn.Taxon(txID)
            self.TaxonDict[txID] = t
            while 1:
                if txID == '1': break
                if not(txID in self.NodesDB): 
                    print 'Nodes is Missing txID %s %s'%(txID, str(type(txID)))
                    break
                parentTx, rank = self.NodesDB[txID]
                name = self.NamesDB[txID]
                pLst = (txID, name, rank)
                t.process(pLst)
                txID = parentTx
    def process_txid(self, txID):
        #print  'process_txid(self, txID: %s)'%str(txID)
        txID = str(txID)
        if txID in self.TaxonDict:
            return self.TaxonDict[txID]
        t = Txn.Taxon(txID)
        self.TaxonDict[txID] = t
        while 1:
            if txID == '1': break
            if not(txID in self.NodesDB): 
                print 'Nodes is Missing txID %s %s'%(txID, str(type(txID)))
                break
            parentTx, rank = self.NodesDB[txID]
            name = self.NamesDB[txID]
            pLst = (txID, name, rank)
            t.process(pLst)
            txID = parentTx
        return t
    def import_TaxonDict(self, inpath):
        for line in fileinput.input(inpath):
            tobj = Txn.TaxObj(line)
            self.TaxonDict[tobj.taxid] = tobj
    def export_TaxonDict(self, outpath='nrtxNeigh.txt'):
        outLst = []
        for k in self.TaxonDict.keys():
            r = self.TaxonDict[k].display()
            outLst.append('\t'.join([str(x) for x in r]))
        ofp = open(outpath, 'w')
        ofp.write('\n'.join(outLst))
        ofp.close()
        return outLst
        
if __name__ == '__main__':
    #txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\'
    txDir = 'C:\\Users\\dominic\\Documents\\Worki7\\Taxonomy_pydev\\taxonomy_db\\taxdmp\\'
    node_db = NodesDB(txDir)
    
    txLst = [
                498736,
                1491,
                1496,
                1502,
                ]
    
    print 'Process set of TaxID\'s'
    txLst = [str(t) for t in txLst]
    node_db.process_txidLst(txLst)
    for t in txLst:
        print node_db.TaxonDict[t].display()
    a = txLst[0]
    b = txLst[1]
    print 'Tax Distance between %s and %s is %i'%(a, b, node_db.TaxonDict[a].Distance(node_db.TaxonDict[b]) )
    node_db.export_TaxonDict('../test_txdict.txt')
    