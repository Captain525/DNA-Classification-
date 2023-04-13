import tensorflow as tf
import torch
from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())
print(tf.config.list_physical_devices())
#print (torch.version.cuda)
#print([torch.cuda.device(i) for i in range(torch.cuda.device_count())])
class SimpleModel(tf.keras.Model):
    def __init__(self, numClasses):
        super().__init__()
        hiddenSize1 = 10000
        self.dense1 = tf.keras.layers.Dense(hiddenSize1, activation = tf.nn.relu)
        self.dense2 = tf.keras.layers.Dense(numClasses, activation = tf.nn.softmax)
    def call(self, inputs, training=False):
        x = self.dense1(inputs)
        x = self.dense2(x)
        return x