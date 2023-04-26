from dataProcessing import *
import tensorflow as tf
from tensorflow.python.client import device_lib
from sklearn.model_selection import train_test_split, KFold
from model import SimpleModel
print(device_lib.list_local_devices())
def runModel():
    #X, Y, vectorizer, dataIndices = loadDataSimpler(4, 4, 10000, 1000)
    
    #n_split = 2
    """
    for train_index, test_index in KFold(n_split).split(X):
        x_train, x_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]

        model = SimpleModel(potential_labels.shape[0])
    
        model.compile(loss="sparse_categorical_crossentropy", optimizer='adam', metrics=['accuracy'])
        model.fit(x_train, y_train, epochs=5)

        print("Model evaluation: ", model.evaluate(x_test, y_test))
    """
    #loadDataSimpler(6,1,10000, 10)
    #xTrain,xVal, xTest, yTrain, yVal, yTest= loadData()
    xTrain, yTrain, xVal, yVal, xTest, yTest, vectorizer = loadDataGeneral(6, 1, 10000, 50, True, False, True)
    print("train size: ", xTrain.shape, yTrain.shape)
    print("val size: ", xVal.shape, yVal.shape)
    print("test size: ", xTest.shape, yTest.shape)
    #maybe should shuffle data here. 
    #hopefully none of the classes aren't ABSENT from this one. 
    numLabels = 4
    model = SimpleModel(numLabels)
    model.compile(loss="binary_crossentropy", optimizer='adam',metrics=['accuracy'])
    print("Y train shape: ", yTrain.shape)
    model.fit(xTrain, yTrain, validation_data = (xTest, yTest), epochs=5)
    print("Model evaluation: ", model.evaluate(xVal, yVal))

runModel()
