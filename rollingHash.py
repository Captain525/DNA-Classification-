
class RollingHashFamily:
    def __init__(self, m, k, sequence, pList):
        
        self.hashFamily = []
        for j in len(pList):
            self.hashFamily.append(iter(RollingHash(m, pList[j], sequence, k)))
    def __iter__(self):
        return self
    def __next__(self):
        hashValues = []
        for iterable in self.hashFamily:
            hashValues.append(next(iterable))
        return hashValues
class RollingHash:
    """
    Split it up into kmers and hash each kmer immediately
    """
    def __init__(self,m,p, sequence, k):
        self.n = len(sequence)
        self.sequence = sequence.lower()
        self.index = 0
        self.stepSize = 1
        self.Hprev = 0
        self.Hcurrent = 0
        self.m = m
        self.p = p
        self.k = k
    def __iter__(self):
        return self



    def __next__(self):
        print("index is: ", self.index)
        if(self.index==0):
            self.Hcurrent = self.calcHash( self.sequence[0:self.k])
            self.Hprev = self.Hcurrent
            self.index = self.index + 1
            return self.Hcurrent
        if(self.index+self.k>self.n):
            return StopIteration
        else:
            oldValue = (ord(self.sequence[self.index - 1]) - ord('a'))*pow(self.p, self.k-1)
            newValue = ord(self.sequence[self.index + self.k-1]) - ord('a')
            self.Hcurrent = ((self.Hprev - oldValue)*self.p + newValue)%self.m
            self.Hprev = self.Hcurrent
            self.index = self.index + 1
            return self.Hcurrent
            
    def calcHash(self,substring):
        pValues = [pow(self.p, k) for k in range(0, self.k)]
        sumValues = [(ord(substring[i])-ord('a'))*pValues[i] for i in range(self.k)]
        sumVals = sum(sumValues)%self.m
        return sumVals
def tryOut():
    sequence = "ACAGGACAGTCTGCC"
    print(ord("A".lower()))
    print(ord("C".lower()))
    print(ord("T".lower()))
    print(ord("G".lower()))
    k = 3
    hash = iter(RollingHash(10, 31, sequence, k))
    
    for i in range(0, 20):
        print(next(hash))
    hash2 = iter(RollingHash(10, 35, sequence, k))
    for i in range(0, 20):
        print(next(hash2))
#tryOut()    