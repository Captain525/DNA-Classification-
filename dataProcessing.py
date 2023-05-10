from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random
import os.path
import time
import mmh3
from sklearn.feature_extraction.text import CountVectorizer
from chooseMethods import *
from rollingHash import RollingHash, RollingHashFamily
from helperMethods import *
from encoding import *
from duplicateMethods import *


def loadDataGeneral(k, n, numData, sequenceSize, checkDuplicates, loadPreviousData, saveData):
    """
    General method to ALWAYS call when loading data, has booleans which can set exactly what will happen. 
    """
    dataFiles = loadInputDataLocations()
    #file locations where data is stored. 
    listDataStored = ["trainX.npy", "trainY.npy", "validationX.npy", "validationY.npy", "testX.npy", "testY.npy"]
    
    #if the boolean and the files exist. 
    if(loadPreviousData and np.all(np.array([os.path.isfile(value) for value in listDataStored]))):
        print("Loading prior data")
        return loadData()
    listSequences = []
    cutoff = []
    listSplits = []
    binList = []
    hashListList = []
    for i in range(0, len(dataFiles)):
        #iterator over sequence objects. Each sequence object is a given size. 
        fasta_sequences = SeqIO.parse(open(dataFiles[i], 'r'), 'fasta')

        #returns training, validation, test sequences. Also returns bins and hash functions if we have the checkDuplicates variable activated.  /
        trainSeq, valSeq, testSeq, bin, hashList = genericPickMethod(fasta_sequences, numData, sequenceSize, checkDuplicates)
        binList.append(bin)
        hashListList.append(hashList)
        indexTrain = len(trainSeq)
        indexVal = len(valSeq) + indexTrain
        sequences = trainSeq+ valSeq + testSeq
        
        listSequences.append(sequences)
        previous = 0
        if(i>0):
            previous = cutoff[i-1]
        cutoff.append(len(sequences) + previous)
        listSplits.append([previous, previous + indexTrain, previous + indexVal, cutoff[i]])
    yData = genY(listSequences, hashListList, binList, checkDuplicates)
    #makes kmers. 
    xData,__, vectorizer = makeCorpusBagOfWords(listSequences, k, n)
    return splitXY(xData, yData, listSplits, saveData)

def genericPickMethod(fasta_sequences, numberToSample, sequenceLength, checkDuplicates):
    """
    A generic pick method which organizes all types of pick methods which will be called. 
    If checkDuplicates is true, it runs minhash to see if everything works. 
    """
    n = 10000
    k = 1
    bin = None
    hashList = None
    simplified = True
    if(simplified):
        fasta_sequences = [fasta_sequences.__next__()]
    if(checkDuplicates):
        
        bin, hashList = minhash(fasta_sequences, sequenceLength,n, k)
    #pick simplified works way worse than pick alternate. 
    if(simplified):
        trainSeq, valSeq, testSeq = pickSimplified(fasta_sequences, numberToSample, sequenceLength)
    else:
         trainSeq, valSeq, testSeq = pickAlternate(fasta_sequences, numberToSample, sequenceLength)
    return trainSeq, valSeq, testSeq, bin, hashList