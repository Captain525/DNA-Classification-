from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random

from sklearn.feature_extraction.text import CountVectorizer

def loadDataSimpler(k,n, numData, sequenceSize):
    dataFiles = loadInputDataLocations()
    listFileSequences = []
    listSplits = []
    cutoff=[]
    for i in range(0, len(dataFiles)):
        print("before sequence")
        fasta_sequences = SeqIO.parse(open(dataFiles[i]), 'fasta')
        trainSequences, valSequences, testSequences = pickSimpler(fasta_sequences, numData, sequenceSize, True)
        indexTrain = len(trainSequences)
        indexVal = len(valSequences) + indexTrain
        sequences = trainSequences + valSequences + testSequences
        #FORGOT THE + PREVIOUS HERE PROBLEMS. 
        
        previous = 0
        if(i>0):
            previous = cutoff[i-1]
        cutoff.append(len(sequences) + previous)
        listFileSequences.append(sequences)
        listSplits.append([previous, previous + indexTrain, previous + indexVal, cutoff[i]])
        
    
    xData, yData, vectorizer = makeCorpusBagOfWords(listFileSequences, k, n)
    trainXList = []
    trainYList = []
    valXList = []
    valYList = []
    testXList = []
    testYList = []
    for i in range(0, len(dataFiles)):
        splits = listSplits[i]
        trainData = xData[splits[0]:splits[1], :]
        trainXList.append(trainData)
        trainDataY = yData[splits[0]:splits[1]]
        trainYList.append(trainDataY)
        valData = xData[splits[1]:splits[2], :]
        valXList.append(valData)
        valDataY = yData[splits[1]:splits[2]]
        valYList.append(valDataY)
        testData = xData[splits[2]:splits[3], :]
        testXList.append(testData)
        testDataY = yData[splits[2]:splits[3]]
        testYList.append(testDataY)
    xTrain = np.vstack(trainXList)
    
    yTrain = np.concatenate(trainYList)
    xVal = np.vstack(valXList)
    yVal = np.concatenate(valYList)
    xTest = np.vstack(testXList)
    yTest = np.concatenate(testYList)
    firstValueData = xData[0,:]
    print("length of xdata: ", firstValueData.shape)
    print("first value of xdata: ", firstValueData)
    sumBoolVals = np.sum(firstValueData.astype(bool))
    print("sum of bool vals: ",sumBoolVals)
    firstLineVectorizer = vectorizer.inverse_transform(xData[0].reshape(1,-1))[0]
    print("size of line vectorizer: ", firstLineVectorizer.shape)
    print("vectorizer first line: ", firstLineVectorizer)
    return xTrain, yTrain, xVal, yVal, xTest, yTest, vectorizer

        
def pickSimpler(fasta_sequences, numberToSample, randomSequenceLength, replaceLater):
    """
    Idea: Have a certain amount of training, val, test data we want. 
    Iterate through. Assign each sequence as a train val or test sequence. 
    sample from that sequence. 
    Once we get past the value we have, we can continue and replace previous data we already sampled to make sure everything
    is truly random and not biased towards the first half of the dataset. 
    """
    #1 over avg desired num of subsequences from a given sequence. 
    pGeom = .05
    probReplace = .5
    percentTrain = .7
    percentVal = .2
    percentTest = .1
    trainNum = percentTrain*numberToSample
    valNum = percentVal * numberToSample
    testNum = percentTest*numberToSample
    trainSeqs = []
    valSeqs = []
    testSeqs = []
    seqTuple = (trainSeqs, valSeqs, testSeqs)
    numTuple = (trainNum, valNum, testNum)
    for sequenceF in fasta_sequences:
        #max number of possible unique subsequences one could take. 
        maxNumberOfSubseqs = len(sequenceF)- randomSequenceLength
        if(maxNumberOfSubseqs<=0):
            continue
        #print("max subsequence length", maxNumberOfSubseqs)
        geomValue = np.random.geometric(pGeom)
        
        value = random.random()
        #0 if train, 1 if val, 2 if test. Also works for val size of 0. 
        choice = int(value>=percentTrain) + int(value>=percentTrain+percentVal)
        #if have enough and don't want to replace anything, skip
        if(numTuple[choice]<= len(seqTuple[choice]) and not replaceLater):
            continue
        #want to take a certain number of subsequences from this sequence. 
        sequence = str(sequenceF.seq)
        #can't take more than max number, no reason to take more than you want. 
        if(maxNumberOfSubseqs>=2147483647):
            maxNumberOfSubseqs = 2147483647-1
        numberSubsequences = min(geomValue, maxNumberOfSubseqs, numTuple[choice])
        
        #pick numSubsequence random subsequences now, of fixed size. 
        #what to do if size bigger than an int? 
        indexList = np.random.randint(0, maxNumberOfSubseqs, numberSubsequences)
        #don't know faster way. 
        for i in range(0, indexList.shape[0]):
            subseq = sequence[indexList[i]:indexList[i] + randomSequenceLength]
            #in this case, replace a random element. Also there is the option to not replace it. 
            if(numTuple[choice]<= len(seqTuple[choice])):
                if(random.random()<probReplace):

                    randIndex = np.random.randint(0, len(seqTuple[choice]))
                    seqTuple[choice][randIndex] = subseq
            else:
                seqTuple[choice].append(subseq)
    return seqTuple[0], seqTuple[1], seqTuple[2]       

    
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
    yArray = np.concatenate(listClassVectors)
    """
    indices = np.arange(0, xArray.shape[0])
    np.random.shuffle(indices)
    xShuffled = xArray[indices, :]
    yShuffled = yArray[indices]
    return xShuffled, yShuffled, vectorizer
    """
    return xArray, yArray, vectorizer
