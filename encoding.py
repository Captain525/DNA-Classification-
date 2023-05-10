
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer

def makeCorpusBagOfWords(listFileSequences, k, n):
    """
    Encodes the sequences into KMERS. 
    """
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