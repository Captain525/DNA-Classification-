from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random
import os.path
import time
import mmh3
from sklearn.feature_extraction.text import CountVectorizer
from chooseMethods import pickAlternate
from rollingHash import RollingHash, RollingHashFamily


def loadDataGeneral(k, n, numData, sequenceSize, checkDuplicates, loadPreviousData, saveData):
    dataFiles = loadInputDataLocations()
    listDataStored = ["trainX.npy", "trainY.npy", "validationX.npy", "validationY.npy", "testX.npy", "testY.npy"]
    
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
        #could maybe use minhash to see which sequences already in data? 

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
    yData = genY(listSequences, hashListList, binList)
    xData,__, vectorizer = makeCorpusBagOfWords(listSequences, k, n)
    return splitXY(xData, yData, listSplits)
def genericPickMethod(fasta_sequences, numberToSample, sequenceLength, replaceLater):
    n = 10000
    k = 100
    bin, hashList = minhash(fasta_sequences, sequenceLength,n, k)
    return pickAlternate(fasta_sequences, numberToSample, sequenceLength, replaceLater), bin, hashList

def genY(listSequences, hashListList, binList):
    listY = []
    for index in range(len(listSequences)):
        sequences = listSequences[index]
        listY.append(checkHashes(sequences, binList, hashListList, index))
    Y = np.vstack(listY)
    return Y
def splitXY(xData, yData, listSplits):
    n = len(listSplits)
    trainXList = []
    trainYList = []
    valXList = []
    valYList = []
    testXList = []
    testYList = []
    for i in range(0, n):
        splits = listSplits[i]
        trainData = xData[splits[0]:splits[1], :]
        trainXList.append(trainData)
        trainDataY = yData[splits[0]:splits[1], :]
        
        trainYList.append(trainDataY)
        valData = xData[splits[1]:splits[2], :]
        valXList.append(valData)
        valDataY = yData[splits[1]:splits[2],:]
        valYList.append(valDataY)
        testData = xData[splits[2]:splits[3], :]
        testXList.append(testData)
        testDataY = yData[splits[2]:splits[3],:]
        testYList.append(testDataY)
    xTrain = np.vstack(trainXList)
    yTrain = np.vstack(trainYList)
    xVal = np.vstack(valXList)
    yVal = np.vstack(valYList)
    xTest = np.vstack(testXList)
    yTest = np.vstack(testYList)
    if(saveData):
        saveData(xTrain,xVal, xTest, yTrain, yVal, yTest)
    return xTrain,yTrain, xVal, yVal, xTest, yTest
def customKmerEncode(k):
    def kmerEncode(sequence):
        return " ".join([sequence[x:x+k].lower() for x in range(len(sequence) - k + 1)])
    return kmerEncode
def kmerEncoding(sequence, k):
    #break up a single sequence into subsequences of size size. 
    return [sequence[x:x+k].lower() for x in range(len(sequence) - k + 1)]
def oneHot(sequence):
    return
def ordinalEncoding(sequence):
    #make the A, C, T, G have related numerical values. 

    return

def encodeSequenceData():
    #ordinal encoding - encode each nitrogen base as ordinal value. 
    #One hot encoding - 4 categories. 
    #Kmer counting - divide it up to get strings of uniform length. 
    return
def loadInputDataLocations():
    with open("dataFiles.txt") as f:
        text = f.readlines()
    links = [link.strip() for link in text]
    print(links)
    return links

def makeCorpusBagOfWords(listFileSequences, k, n):
    vectorizer = CountVectorizer(input="content", ngram_range= (n,n), dtype= np.int8)
    listVocab = []
    listkmers =[]
    for sequenceList in listFileSequences:
        kmerFunction = customKmerEncode(k)
        print("Sequence length: ", len(sequenceList))
        kmerSequence = list(map(kmerFunction, sequenceList))
        assert(len(sequenceList)==len(kmerSequence))
        listVocab+=kmerSequence
        listkmers.append(kmerSequence)
    vectorizer.fit(listVocab)
    listArrays = []
    listClassVectors = []
    print("length of kmer list: ", len(listkmers))
    for i in range(0, len(listkmers)):
        sequenceList = listkmers[i]
        #should be numSequences by vocabSize
        print("size of seq list: ", len(sequenceList))
        X = vectorizer.transform(sequenceList).toarray()
        print("Post vectorized sequence list length: ", X.shape)
        classVector = i*np.ones(shape=(X.shape[0], ))
        listArrays.append(X)
        listClassVectors.append(classVector)
        print("X size is: ", X.shape)
    xArray = np.vstack(listArrays)
    rowSum = np.sum(xArray, axis=1)
    #normalization
    xArray = xArray/rowSum[:, np.newaxis]
    yArray = np.concatenate(listClassVectors)
    return xArray, yArray, vectorizer
def checkClassesForExamples(bearType, classSequences):
    """
    Given a sequence of data from a class, checks if those sequences occur in any other
    classes, to check if need multiple labels. 

    Classequences is a list of sequences, where each sequence is a string. 

    """
    dataFiles = loadInputDataLocations()
    arrayIn = np.zeros(shape = (len(classSequences), 4))
    arrayIn[:, bearType] = np.ones(shape = (len(classSequences)))
    for i in range(len(dataFiles)):
        if(i == bearType):
            continue
        fasta_sequences = SeqIO.parse(open(dataFiles[i]), 'fasta')
        for sequenceRecord in fasta_sequences:
            sequence = sequenceRecord.seq
            print("sequence length: ", len(sequence))
            found = False
            for j in range(len(classSequences)):
                if(found):
                    break
                if(j%100 == 0):
                    print("on sequence: ", j,"/", len(classSequences),"\n")
                if(sequence.find(classSequences[j])!=-1):
                    arrayIn[j, i] = 1
                    found = True
    return arrayIn
def checkClassesForExamplesAlternate(bearType, classSequences):
    """
    MUCH FASTER THAN standard. Standard straight up doesn't finish. This one 
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
    stepsize = 4
    for i in range(0, (len(sequence) - k)//stepsize):
        print("I value: ", i)
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
        print("On sequence num: ", count)
        sequenceString = str(sequence.seq)
        print("length: ", len(sequenceString))
        gen = generatorSequence(sequenceString, sequenceSize)
        while(True):
            try:
                kmer = next(gen)
               # print("Kmer: ", kmer, "\n")

                #generates values in huge range. 
                hashValues = hashKmer(kmer, hashList, n)
                #print("hash values: ", hashValues)
                bin[hashValues] = 1
            except StopIteration:
                break
        count = count+1
    print("done minhash")
    endTime = time.time()
    print("Total time: ", endTime - startTime)
    return bin, hashList

def checkHashes(sequences, binList, hashListList, index):
    Y = np.zeros(shape = (len(sequences), 4))
    Y[:, index] = np.ones(shape = (len(sequences)))
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
        return bin[hashValues[i]]
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

def saveData(trainX, validationX, testX, trainY, validationY, testY):
    print("Shape of training data: ", trainX.shape)
    print("shape of Y data: ", trainY.shape)
    np.save("trainX", trainX)
    np.save("validationX", validationX)
    np.save("testX", testX)
    np.save("trainY", trainY)
    np.save("validationY", validationY)
    np.save("testY", testY)
def loadData():
    trainX = np.load("trainX.npy")
    validationX = np.load("validationX.npy")
    testX = np.load("testX.npy")
    trainY = np.load("trainY.npy")
    validationY = np.load("validationY.npy") 
    testY = np.load("testY.npy")
    assert(trainY.shape)
    return trainX, trainY, validationX, validationY, testX, testY
def callLoad():
    trainX, valX, testX, trainY, valY, testY = loadData()
    print("train x shape: ", trainX.shape)
    print("trainy shape: ", trainY.shape)
    print(trainX[0:5])
    print(trainY[0:5])
    avgNumberClasses = np.mean(np.sum(trainY, axis=1))
    print(avgNumberClasses)
#callLoad()