import numpy as np
from tensorflow import keras
from vae.vae_class import VAE

print('Loading model...')
encoder = keras.models.load_model(enc_save_folder)
decoder = keras.models.load_model(dec_save_folder)

vae = VAE(encoder, decoder)
print('Compiling model...')
vae.compile(optimizer=keras.optimizers.Adam())
print('Model loaded.')