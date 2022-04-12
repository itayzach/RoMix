import numpy as np
from tensorflow import keras
from vae.vae_class import VAE

encoder = keras.models.load_model(enc_save_folder)
decoder = keras.models.load_model(dec_save_folder)

vae = VAE(encoder, decoder)
vae.compile(optimizer=keras.optimizers.Adam())