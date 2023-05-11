from dataProcessing import *
import tensorflow as tf
from tensorflow.python.client import device_lib
from sklearn.model_selection import train_test_split, KFold
from model import SimpleModel
print(device_lib.list_local_devices())
def runModel():
    xTrain, yTrain, xVal, yVal, xTest, yTest = loadDataGeneral(6, 1, 8000, 100, checkDuplicates=False, 
                                                               loadPreviousData = False, saveData = True, 
                                                               doMultiple = False, simplified = True)
    print("train size: ", xTrain.shape, yTrain.shape)
    print("val size: ", xVal.shape, yVal.shape)
    print("test size: ", xTest.shape, yTest.shape)
    
    numLabels = 4
    model = SimpleModel(numLabels)
    model.compile(loss="binary_crossentropy", optimizer='adam',metrics=['accuracy'])
    print("Y train shape: ", yTrain.shape)
    model.fit(xTrain, yTrain, validation_data = (xVal, yVal), epochs=5)
    print("Model evaluation: ", model.evaluate(xTest, yTest))

runModel()
