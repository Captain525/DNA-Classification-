import tensorflow as tf
class SimpleModel(tf.keras.Model):
    def __init__(self, numClasses):
        super().__init__()
        hiddenSize1 = 100
        self.dense1 = tf.keras.layers.Dense(hiddenSize1, activation = tf.nn.relu)
        self.dense2 = tf.keras.layers.Dense(numClasses, activation = tf.nn.softmax)
    def call(self, inputs, training=False):
        x = self.dense1(inputs)
        x = self.dense2(x)
        return x