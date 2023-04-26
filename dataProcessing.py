from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random
import time

from sklearn.feature_extraction.text import CountVectorizer

def loadDataSimpler(k,n, numData, sequenceSize):
    dataFiles = loadInputDataLocations()
    listFileSequences = []
    listSplits = []
    cutoff=[]
    listY = []
    for i in range(0, len(dataFiles)):
        print("before sequence")
        fasta_sequences = SeqIO.parse(open(dataFiles[i]), 'fasta')
        trainSequences, valSequences, testSequences = pickAlternate(fasta_sequences, numData, sequenceSize, True)
        indexTrain = len(trainSequences)
        indexVal = len(valSequences) + indexTrain
        sequences = trainSequences + valSequences + testSequences
        #FORGOT THE + PREVIOUS HERE PROBLEMS. 
      
        """
        trainY = checkClassesForExamplesAlternate(i, trainSequences)
        valY = checkClassesForExamplesAlternate(i, valSequences)
        testY = checkClassesForExamplesAlternate(i, testSequences)
        """
        #ignoring this problem for now. 
        #"""
        trainY = np.zeros(shape = (len(trainSequences), 4))
        trainY[:, i] = np.ones(shape = (len(trainSequences)))
        valY = np.zeros(shape = (len(valSequences), 4))
        valY[:, i] = np.ones(shape = (len(valSequences)))
        testY = np.zeros(shape = (len(testSequences), 4))
        testY[:, i] = np.ones(shape = (len(testSequences)))
        #"""
        listY.append([trainY, valY, testY])
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
        #trainDataY = yData[splits[0]:splits[1]]
        trainDataY = listY[i][0]
        trainYList.append(trainDataY)
        valData = xData[splits[1]:splits[2], :]
        valXList.append(valData)
        #valDataY = yData[splits[1]:splits[2]]
        valDataY = listY[i][1]
        valYList.append(valDataY)
        testData = xData[splits[2]:splits[3], :]
        testXList.append(testData)
        #testDataY = yData[splits[2]:splits[3]]
        testDataY = listY[i][2]
        testYList.append(testDataY)
    xTrain = np.vstack(trainXList)
    
    #yTrain = np.concatenate(trainYList)
    yTrain = np.vstack(trainYList)
    xVal = np.vstack(valXList)
    #yVal = np.concatenate(valYList)
    yVal = np.vstack(valYList)
    xTest = np.vstack(testXList)
    #yTest = np.concatenate(testYList)
    yTest = np.vstack(testYList)
    saveData(xTrain,xVal, xTest, yTrain, yVal, yTest)
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
    pGeom = .001
    probReplace = .1
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
        indexList = np.random.randint(0, maxNumberOfSubseqs, size = numberSubsequences)
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
def pickAlternate(fasta_sequences, numberToSample, randomSequenceLength, replaceLater):
    trainPercent = .7
    valPercent = .2
    numSequenceSample = int(numberToSample/randomSequenceLength)
    listSeqs = []
    for sequenceF in fasta_sequences:
        sequence = str(sequenceF.seq)
        randomIndex = np.random.randint(0, len(sequence) - randomSequenceLength)
        subseq = sequence[randomIndex: randomIndex + randomSequenceLength]
        listSeqs.append(subseq)
    n = len(listSeqs)
    random.shuffle(listSeqs)
    trainSeq = listSeqs[0:int(n*trainPercent)]
    valSeq = listSeqs[int(n*trainPercent):int(n*(trainPercent + valPercent))]
    testSeq = listSeqs[int(n*(trainPercent + valPercent)):]
    return trainSeq, valSeq, testSeq
    
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
    """
    indices = np.arange(0, xArray.shape[0])
    np.random.shuffle(indices)
    xShuffled = xArray[indices, :]
    yShuffled = yArray[indices]
    return xShuffled, yShuffled, vectorizer
    """
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
    return trainX, validationX, testX, trainY, validationY, testY
def callLoad():
    trainX, valX, testX, trainY, valY, testY = loadData()
    print("train x shape: ", trainX.shape)
    print("trainy shape: ", trainY.shape)
    print(trainX[0:5])
    print(trainY[0:5])
    avgNumberClasses = np.mean(np.sum(trainY, axis=1))
    print(avgNumberClasses)
#callLoad()