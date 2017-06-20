'''
Created on Sep 6, 2012

@author: dominic
'''

class TaxObj:
    #This one is simply r-enstated from a line output
    #ToDO melt this back into the other one
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

class Taxon:
    def __init__(self, txID):
        #print "Taxon: %s"%str(txID)
        self.tx = txID
        self.txLst = []
        self.nameLst = []
        self.rankLst = []
        self.n = 0
        self.genusLoc = None
        self.no_rank = None
        self.subspecies = None
        self.species = None
        self.genus = None
        self.family = None
        self.order = None
        self.species_name = 'uk'
        self.subspecies_name = 'uk'
        self.no_rank_name = 'uk'
        self.genus_name = 'uk'
    def process(self, inLst):
        txid, name, rank = inLst[:3]
        self.txLst.append(txid)
        self.nameLst.append(name)
        self.rankLst.append(rank)
        self.n +=1
        if self.no_rank == None and rank == 'no rank':
            self.no_rank = txid
            self.no_rank_name = name.replace('  ', ' ').replace(' ', '_')
        if rank == 'subspecies':
            self.subspecies = txid
            self.subspecies_name = name.replace('  ', ' ').replace(' ', '_')
        if rank == 'species':
            self.species = txid
            self.species_name = name.replace('  ', ' ').replace(' ', '_')
        elif rank == 'genus':
            self.genus = txid
            self.genus_name = name
        elif rank == 'family':
            self.family = txid
        elif rank == 'order':
            self.order = txid
    def speciesLabel(self):
        l = "genusTxID[%s]_%s[%s]"%(str(self.genus), self.species_name, str(self.species))
        return l
    def subspeciesLabel(self):
        if self.subspecies:
            l = "genusTxID[%s]_%s[%s]_%s[%s]"%(str(self.genus), self.species_name, str(self.species), self.subspecies_name, str(self.subspecies))
            return l
        else:
            return self.speciesLabel()
    def no_rankLabel(self):
        if self.no_rank:
            l = "genusTxID[%s]_%s[%s]_%s[%s]_%s[%s]"%(str(self.genus), self.species_name, str(self.species), self.subspecies_name, str(self.subspecies), self.no_rank_name, str(self.no_rank))
            return l
        else:
            return self.subspeciesLabel()
    def process_line(self, line):
        #This allows you to import this from an tax neighborhood archive
        line = line.replace('\n', '')
        Lst = line.split('\t')
        self.tx = self.Lst.pop(0)
        self.process(Lst)
    def gLoc(self):
        if 'genus' in self.rankLst:
            self.genusLoc = self.rankLst.index('genus')
        else:
            self.genusLoc = -1
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
    def display(self):
        oLst = []
        for i in range(self.n):
            oLst.append('%s__%s__%s' %(self.txLst[i], self.nameLst[i], self.rankLst[i]))
        self.gLoc()
#        if self.genusLoc != -1:
#            tabAdd = max(0, 10-self.genusLoc)
#            for i in range(tabAdd):
#                oLst = ['']  + oLst    
        oLst = [self.tx]  + oLst    
        return oLst
    


if __name__ == '__main__':
    txID = 832
    t = Taxon(txID)
    

    