import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from vae.vae_class import VAE

#print(data.shape)
#print(vae.encoder.summary())
data = np.expand_dims(data, -1).astype("float32") #/ 255
z_mean, z_log_var, z = vae.encoder.predict(data)