

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