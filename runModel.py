from dataProcessing import *
import tensorflow as tf
from sklearn.model_selection import train_test_split, KFold
from model import SimpleModel

def runModel():
    X, Y, vectorizer = loadData()
    potential_labels = np.unique(Y)
    n_split = 2
    for train_index, test_index in KFold(n_split).split(X):
        x_train, x_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]

        model = SimpleModel(potential_labels.shape[0])
    
        model.compile(loss="sparse_categorical_crossentropy", optimizer='adam', metrics=['accuracy'])
        model.fit(x_train, y_train, epochs=5)

        print("Model evaluation: ", model.evaluate(x_test, y_test))


runModel()
