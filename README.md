# RoMix
This is the official repo to reproduce the results of our paper "Graph Signals Interpolation and Extrapolation Using RoMix - Reproducing Kernel Hilbert Space over Manifold of Gaussian Mixture"

## Clone
This repo contains [GSPBox](https://epfl-lts2.github.io/gspbox-html/) as a submodule. Be sure to clone with `--recursive` or ``git submodule update`` after you clone.

## Setup
All simulations were conducted with MATLAB 2021b on Windows 10. Other OS's or MATLAB versions might not be supported.

For MNIST simulation:
* Training VAE was done with Python 3.9.0 and the packages in requirements.txt. You may install them by (Anaconda 2.1.1):
  ```
  conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
  pip install -r requirements.txt
  ```
* Verify:
  ```
  # CPU:
  python -c "import tensorflow as tf; print(tf.reduce_sum(tf.random.normal([1000, 1000])))"
  # GPU:
  python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
  ```
* You must link your MATLAB to your python installation.

  
## How to run simulations
Your working directory should be `src\`. 

On MATLAB startup, run `startup.m` for all required addpaths.

* For a single preset, run `Main.m` with the desired preset struct returned from each `Get<>Preset` function, e.g.:
  ```
  Main(GetBulgariBeaconsPreset);
  ```
* Reproduce paper simulations:
   ```
   RunAllPaperPresets;
   ```
* Run all presets:
   ```
   RunAllPresets;
   ```
