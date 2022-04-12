from tensorflow import keras
(x_train, y_train), (x_test, y_test) = keras.datasets.mnist.load_data()
x_train = x_train[:num_train,:,:].astype("float32") / 255
y_train = y_train[:num_train].astype("float32") / 255
x_test = x_test[:num_test,:,:]
y_test = y_test[:num_test]
