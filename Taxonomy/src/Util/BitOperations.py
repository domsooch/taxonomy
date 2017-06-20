
# testBit() returns a nonzero result, 2**offset, if the bit at 'offset' is one.

def testBit(int_type, offset):
    mask = 1 << offset
    return(int_type & mask)

# setBit() returns an integer with the bit at 'offset' set to 1.

def setBit(int_type, offset):
    mask = 1 << offset
    return(int_type | mask)

# clearBit() returns an integer with the bit at 'offset' cleared.
def clearBit(int_type, offset):
    mask = ~(1 << offset)
    return(int_type & mask)

# toggleBit() returns an integer with the bit at 'offset' inverted, 0 -> 1 and 1 -> 0.

def toggleBit(int_type, offset):
    mask = 1 << offset
    return(int_type ^ mask)

class BitArray:
    def __init__(self):
        self.v = 0
        self.max = 0
    def setBit(self, offset):
        mask = 1 << offset
        if offset > self.max:
            self.max = offset
        self.v = self.v | mask
    def clearBit(self, offset):
        mask = ~(1 << offset)
        self.v = self.v & mask
    def testBit(self, offset):
        mask = 1 << offset
        r = self.v & mask
        if r !=0:
            return 1
        else:
            return 0
    def __len__(self):
        s = 0
        for i in range(self.max+1):
            s+= self.testBit(i)
        return s
    def list(self):
        oLst = []
        for i in range(self.max+1):
            if self.testBit(i):
                oLst.append(i)
        return oLst
    
        
        
        
        
if __name__ == "__main__":
    vals = [12, 10, 11, 3, 2]
    a = 0
    for v in vals:
        a = setBit(a, v)
    
    for i in range(max(vals)+1):
        print i, testBit(a, i)
        
    b = BitArray()
    for v in vals:
        b.setBit(v)
    print 'count: ', len(b)
    print 'list: ', b.list()