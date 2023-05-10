from helperMethods import *
import numpy as np
import mmh3
from Bio import SeqIO
from Bio.Seq import Seq
import time
from rollingHash import RollingHash, RollingHashFamily
def checkClassesForExamplesAlternate(bearType, classSequences):
    """
    Instead of using minhash, this instead just iterates through the sequences we have of a given bear type and checks if 
    there are duplicates in the other classes. 
    """
    dataFiles = loadInputDataLocations()
    arrayIn = np.zeros(shape = (len(classSequences), 4))
    arrayIn[:, bearType] = np.ones(shape = (len(classSequences)))
    print(arrayIn)
    for i in range(len(dataFiles)):
        if(i == bearType):
            continue
        fasta_sequences = SeqIO.parse(open(dataFiles[i]), 'fasta')
        for j in range(len(classSequences)):
            found = False
            if(j%100 == 0):
                    print("on sequence: ", j,"/", len(classSequences),"\n")
            classSeq = classSequences[j]
            for sequenceRecord in fasta_sequences:
                sequence = sequenceRecord.seq
                if(sequence.find(classSeq) !=-1):
                    print("found it")
                    found = True
                    arrayIn[j, i] = 1
                    break
    return arrayIn
def generatorSequence(sequence, k):
    """
    Yields the kmers when we want them, so we don't need to store them all in memory. 
    """
   
    stepsize = 1
    
    for i in range(0, (len(sequence) - k)):
        if(i%10000000 ==0):
            print("I value: ", i," / ", len(sequence) - k)
        yield sequence[stepsize*i:stepsize*i+k]

def minhashWithRollingHash(fasta_sequences, sequenceSize, n, k):
    timeStart = time.time()
    pList = range(20, 20+k)
    bin = np.zeros(shape =(n, ), dtype=bool)
    count = 0
    for sequenceF in fasta_sequences:
        print("count: ", count)
        sequence = str(sequenceF.seq)
        hashFamily = iter(RollingHash(n, pList[0], sequence, sequenceSize))
        count = count+1
        try:
            while(True):
                bin[next(hashFamily)]=1
        except StopIteration:
            break
        
    timeEnd = time.time()
    print("minhash with roll took: ", timeEnd-timeStart)
    return bin, pList    

           
def minhash(fasta_sequences, sequenceSize,n, k):
    """
    Goal of this is to keep a storage of all possible kmers of a certain size, ie the size of sequences we want to draw, in a given set. 
    Need a way to be memory efficient, so treating it like a bloom filter. 

    n - number of bits. 
    k - number of hash functions, independent. 

    Ignoring reverse complements for now. 
    """
    print("starting minhash")
    startTime = time.time()
    hashList = range(0, k)
    bin = np.zeros(shape =(n, ), dtype=bool)
    count = 0
    for sequence in fasta_sequences:
        startTimeMini = time.time()
        print("On sequence num: ", count)
        sequenceString = str(sequence.seq)
        print("length: ", len(sequenceString))
        gen = generatorSequence(sequenceString, sequenceSize)
        while(True):
            try:
                kmer = next(gen)

                #generates values in huge range. 
                hashValues = hashKmer(kmer, hashList, n)
                bin[hashValues] = 1
            except StopIteration:
                break
        count = count+1
        endTimeMini = time.time()
        print("Time for one sequence: ", endTimeMini- startTimeMini)
    print("done minhash")
    endTime = time.time()
    print("Total time: ", endTime - startTime)
    return bin, hashList

def checkHashes(sequences, binList, hashListList, index, checkDuplicates):
    """
    Check if eqch element of the sequence is in other classes or not. 
    """
    Y = np.zeros(shape = (len(sequences), 4))
    Y[:, index] = np.ones(shape = (len(sequences)))
    if(not checkDuplicates):
        return Y
    for i in range(len(binList)):
        if(i == index):
            continue
        bin = binList[i]
        hashList = hashListList[i]
        Y[:, i] = checkHashesBin(sequences, bin, hashList)
    return Y
def hashesBinGenerator(sequence, bin, hashList):
    k = len(hashList)
    n = bin.shape[0]
    hashValues = hashKmer(sequence, hashList, n)
    for i in range(0, k):
        yield bin[hashValues[i]]
def checkHashesBin(sequences, bin, hashList):
    """
    Checks if each sequence element is in the given bin, ie in the set. 
    """
    Y = np.zeros(shape = (len(sequences), ))
    for j in range(len(sequences)):
        sequence = sequences[j]
        try:
            gen = hashesBinGenerator(sequence, bin, hashList)
            while(next(gen)):
                continue
        except StopIteration:
            Y[j] = 1
    return Y


def hashKmer(kmer, hashList, n):
    hashValues = [mmh3.hash64(kmer, hashList[i])[0]%n for i in hashList]
    return hashValues
