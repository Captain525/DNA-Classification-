from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import random

from sklearn.feature_extraction.text import CountVectorizer
def loadData():
    dataFiles = loadInputDataLocations()
    listFileSequences = []
    for i in range(0, len(dataFiles)):
        print("before sequence")
        listSequencesDataFiles = loadFasta(dataFiles[i])
        listFileSequences.append(listSequencesDataFiles)
    #keep vectorizer in case want to interpret the data. 
    #lose the sequential nature of given subsequences doing it this way. 
    xData, yData, vectorizer = makeCorpusBagOfWords(listFileSequences, 4, 4)
    print(xData.shape)
    print(yData.shape)
    return xData, yData, vectorizer
def loadFasta(input_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    
    #Should be a list of strings. 
    listSequences = pickFromEachRandomly(fasta_sequences, 10000, 1000)
    return listSequences


def pickFromEachRandomly(fasta_sequences, numberToSample, randomSequenceLength):
    """
    Iterate through the list of fasta sequences, go through each and pick a random 
    subset of that sequence of a certain size. 
    Can't go through all of them? 
    randomSequenceLength - how long we want each sequence to be. 
    
    """
    p = .2
    pGeom = .05
    count = 0
    numSamples = 0
    listSequences = []
    for fasta in fasta_sequences:
        if(numSamples>=numberToSample):
            break
        #with probability 1-p skip this one. 
        dontChooseThis = random.random()>p
        if(dontChooseThis):
            continue
        name, sequence = fasta.id, str(fasta.seq)
        #sequenceArray = np.array(list(sequence), dtype=str)
        maxNumberOfSubseqs = len(sequence)- randomSequenceLength
        #print("Max number of subseq: ", maxNumberOfSubseqs)
        #number of trials UNTIL a success, so want pGeom to be small. 
        #print("random geometric value: ", np.random.geometric(pGeom))
        numberOfSubsequences = min(np.random.geometric(pGeom), maxNumberOfSubseqs)
        #print("number of subsequences: ", numberOfSubsequences)
        #generate a random integer from 0 to maxNumberOfSubseqs exclusive. This is valid starting indices. 
        arrayIndices = np.random.randint(0, maxNumberOfSubseqs, numberOfSubsequences)
        sequenceList = []
        for i in range(0, arrayIndices.shape[0]):
            sequenceList.append(sequence[arrayIndices[i]:arrayIndices[i] + randomSequenceLength])
        listSequences+=sequenceList

        #can't think of a faster way to do this for now
        count+=1
        numSamples+=numberOfSubsequences
    avgSubstringLength = numSamples/count
    print("average substring length is: ", avgSubstringLength)
    return listSequences
def pickSimpler(fasta_sequences, numberToSample, randomSequenceLength):
    size = len(fasta_sequences)
    listSequences = random.shuffle(fasta_sequences)
    percentTrain = .8
    pGeom = .05
    trainCutoff = int(size*percentTrain)
    testCutoff = size - trainCutoff
    #divides the total number of desired subsequences by the expected value. 
    numberSequences = min(int(numberToSample*pGeom), trainCutoff)
    trainSubset = listSequences[0:numberSequences]
    testSize = testCutoff
    testSubset = listSequences[-testSize:]
    trainSequences = helperSampleFromSequence(trainSubset, randomSequenceLength, numberToSample, pGeom)

    testSequences = helperSampleFromSequence(testSubset, randomSequenceLength,numberToSample, pGeom)
    return trainSequences, testSequences
    
def helperSampleFromSequence(subset, randomSequenceLength, numberToSample, pGeom):
    numSamples = 0
    count = 0
    listSequences = []
    for sequence in subset:
        if(numSamples>=numberToSample):
            break
        #with probability 1-p skip this one. 
       
        sequence = str(sequence.seq)
        
        maxNumberOfSubseqs = len(sequence)- randomSequenceLength
        numberOfSubsequences = min(np.random.geometric(pGeom), maxNumberOfSubseqs)
        #generate a random integer from 0 to maxNumberOfSubseqs exclusive. This is valid starting indices. 
        arrayIndices = np.random.randint(0, maxNumberOfSubseqs, numberOfSubsequences)
        sequenceList = []
        for i in range(0, arrayIndices.shape[0]):
            sequenceList.append(sequence[arrayIndices[i]:arrayIndices[i] + randomSequenceLength])
        listSequences+=sequenceList

        #can't think of a faster way to do this for now
        count+=1
        numSamples+=numberOfSubsequences
    avgSubstringLength = numSamples/count
    print("average substring length is: ", avgSubstringLength)
    return listSequences
def loadFastaSimple(input_file):
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    transformedSequenceData = transformEach(fasta_sequences)

def transformEach(sequences):
    """
    Alternative where it does the transformations to the data immediately. 
    """
    listJoined = []
    for sequence in sequences:
        #join a sequence of words into one long string with spaces to make it so that you can extract. 
        kmerString = " ".join(kmerEncoding(sequence,4))
        listJoined.append(kmerString)
    vectorizer = CountVectorizer()
    X = vectorizer.fit( raw_documents = listJoined)
    print("shape of X: ", X.shape)
    return X
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
    print("done fit")
    listArrays = []
    listClassVectors = []
    
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
    indices = np.arange(0, xArray.shape[0])
    np.random.shuffle(indices)
    xShuffled = xArray[indices, :]
    yShuffled = yArray[indices]
    return xShuffled, yShuffled, vectorizer

