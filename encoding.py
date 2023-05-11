
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import LabelBinarizer
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
def makeOneHot(lb):
    def oneHot(sequence):
        return lb.transform(list(sequence))
    return oneHot
def ordinalEncoding(sequence):
    #make the A, C, T, G have related numerical values. 
    #Let A = 1 C = 2 T = 3 G = 4
    return np.array(list(map(ordinalEncode, list(sequence.lower()))))
    
def ordinalEncode(c):
    if(c=="a"):
        return 1
    elif(c=="c"):
        return 2
    elif(c=="t"):
        return 3
    elif(c=="g"):
        return 4
    else:
        return 0

def encodeSequenceData(listFileSequences, choice, k):
    #ordinal encoding - encode each nitrogen base as ordinal value. 
    #One hot encoding - 4 categories. 
    #Kmer counting - divide it up to get strings of uniform length. 
    if(choice == 0):
        #one hot encoding. 
        lb = LabelBinarizer()
        lb.fit(['a', 'c', 't', 'g'])
        oneHotMethod = makeOneHot(lb)
        arrayList =[]
        for sequenceList in listFileSequences:
            array = np.row_stack(list(map(oneHotMethod, sequenceList)))
            arrayList.append(array)
        encodedSequences = np.concatenate(arrayList, axis=0)
        print("encoded sequence size: ", encodedSequences.shape)
        return encodedSequences
    elif(choice == 1):
        #ordinal encoding
        arrayList = []
        for sequenceList in listFileSequences:
            array = np.row_stack(list(map(ordinalEncoding, sequenceList)))
            arrayList.append(array)
        encodedSequences = np.en
    return