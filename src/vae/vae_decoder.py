import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from vae.vae_class import VAE

z = np.atleast_2d(z)
print('z shape = ' + str(z.shape))
#print(vae.decoder.summary())
x = vae.decoder.predict(z)

