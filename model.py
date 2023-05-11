import tensorflow as tf

from tensorflow.python.client import device_lib

class SimpleModel(tf.keras.Model):
    def __init__(self, numClasses):
        super().__init__()
        hiddenSize1 = 1000
        hiddenSize2= 100
        self.dense1 = tf.keras.layers.Dense(hiddenSize1, activation = tf.nn.relu)
        self.dense2 = tf.keras.layers.Dense(hiddenSize2, activation = tf.nn.relu)
        self.dense3 = tf.keras.layers.Dense(numClasses, activation = "sigmoid")
        #self.dense2 = tf.keras.layers.Dense(numClasses, "softmax")
    def call(self, inputs, training=False):
        print("input shape: ", inputs.shape)
        x = self.dense1(inputs)
        x = self.dense2(x)
        x = self.dense3(x)
        print("x shape: ", x. shape)
        return x