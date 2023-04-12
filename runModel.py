from dataProcessing import *
import tensorflow as tf
from sklearn.model_selection import train_test_split, KFold
from model import SimpleModel

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
    xTrain, yTrain, xVal, yVal, xTest, yTest, vectorizer = loadDataSimpler(5,1,10000, 1000)
    print("train size: ", xTrain.shape, yTrain.shape)
    print("val size: ", xVal.shape, yVal.shape)
    print("test size: ", xTest.shape)
    print("Y test size: ", yTest.shape)
    #maybe should shuffle data here. 
    #hopefully none of the classes aren't ABSENT from this one. 
    potential_labels = np.unique(yTrain)
    model = SimpleModel(potential_labels.shape[0])
    model.compile(loss="sparse_categorical_crossentropy", optimizer='adam', metrics=['accuracy'])
    model.fit(xTrain, yTrain, validation_data = (xTest, yTest), epochs=5)
    print("Model evaluation: ", model.evaluate(xVal, yVal))

runModel()
