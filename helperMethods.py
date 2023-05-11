import numpy as np
from duplicateMethods import *
def genY(listSequences, hashListList, binList, checkDuplicates, doMultiple):
    """
    Generates Y data. Allows for versatility with the minhash algorithm, as well as doing it a slightly different way. 
    """
    listY = []
    #goes across the whole length of total sequences. 
    for index in range(len(listSequences)):
        sequences = listSequences[index]
        if(checkDuplicates):
            listY.append(checkHashes(sequences, binList, hashListList, index, checkDuplicates))
        elif(doMultiple):
            listY.append(checkClassesForExamplesAlternate(index, sequences))
        else:
            #just one label. 
            zeros = np.zeros(shape=(len(sequences), len(listSequences)))
            zeros[:, index] = np.ones(shape = (len(sequences), ))
            listY.append(zeros)
    Y = np.vstack(listY)
    return Y
def splitXY(xData, yData, listSplits, saveDataBool):
    """
    Splits X data and Y data into training, validation, and test according to the splits generated earlier. 
    Assumes the size of xDAta is numExamplesxdims, yData is numExamples x 4. 

    """
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
    xTrain, yTrain = shuffleData(xTrain, yTrain)
    xVal = np.vstack(valXList)
    yVal = np.vstack(valYList)
    xVal, yVal = shuffleData(xVal, yVal)
    xTest = np.vstack(testXList)
    yTest = np.vstack(testYList)
    xTest, yTest = shuffleData(xTest, yTest)
   
    if(saveDataBool):
        saveData(xTrain,xVal, xTest, yTrain, yVal, yTest)
    return xTrain,yTrain, xVal, yVal, xTest, yTest
def shuffleData(x,y):
    shuffleIndices = np.arange(0, x.shape[0])
    np.random.shuffle(shuffleIndices)
    newX = x[shuffleIndices, :]
    newY = y[shuffleIndices,:]
    return newX, newY
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
def loadInputDataLocations():
    with open("dataFiles.txt") as f:
        text = f.readlines()
    links = [link.strip() for link in text]
    print(links)
    return links
def callLoad():
    trainX, valX, testX, trainY, valY, testY = loadData()
    print("train x shape: ", trainX.shape)
    print("trainy shape: ", trainY.shape)
    print(trainX[0:5])
    print(trainY[0:5])
    avgNumberClasses = np.mean(np.sum(trainY, axis=1))
    print(avgNumberClasses)

