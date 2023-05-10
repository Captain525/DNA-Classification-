import random
import numpy as np
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
def pickAlternate(fasta_sequences, numberToSample, randomSequenceLength):
    """
    Just picks one subsequence form each sequence. 
    """
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

def pickSimplified(fasta_sequences, numberToSample, randomSequenceLength):
    """
    A pick method designed to lessen the amount of data we need to deal with. One way to do this is only sample from the first
    contig for each sequence. 
    """
    sequence = str(fasta_sequences[0].seq)
    trainPercent = .7
    valPercent = .2
    listSeqs = []
    for i in range(numberToSample):
        randomIndex = np.random.randint(0, len(sequence) - randomSequenceLength)
        subseq = sequence[randomIndex: randomIndex + randomSequenceLength]
        listSeqs.append(subseq)
    n = len(listSeqs)
    random.shuffle(listSeqs)
    trainSeq = listSeqs[0:int(n*trainPercent)]
    valSeq = listSeqs[int(n*trainPercent):int(n*(trainPercent + valPercent))]
    testSeq = listSeqs[int(n*(trainPercent + valPercent)):]
    return trainSeq, valSeq, testSeq