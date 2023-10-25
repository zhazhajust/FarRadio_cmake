#!/usr/bin/env python3
import os
import yaml
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyfaradio.faradio import SpheDetector
from numpy import pi

def load_screen(filename):
    with open(filename, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    return data

def divide_chunk(rank, size, shape):
    chunksize = shape // size
    return rank*chunksize, min(shape, (rank + 1)*chunksize)

def main(filename: str):
    screen = load_screen("./screen.yaml")
    os.environ['OMP_NUM_THREADS'] = str(screen["omp_num_threads"])
    det = SpheDetector(screen["dmin"], screen["dmax"], screen["nf"])
    det.set_approx(screen["if_approx"])
    axis_time = np.linspace(screen["dmin"][0], screen["dmax"][0], screen["nf"][0])
    axis_phi = np.linspace(screen["dmin"][2], screen["dmax"][2], screen["nf"][2])

    print("Script name:", sys.argv[0])
    # Print the command-line arguments
    print("Arguments:", sys.argv[1:])

    # Read Data
    with h5py.File(filename, "r") as f:
        time = f["time"][...]

    with h5py.File(filename, "r") as f:
        position_cur = f["position_cur"][...]
        position_prev = f["position_prev"][...]
        beta_cur = f["beta_cur"][...]
        beta_prev = f["beta_prev"][...]
        time = f["time"][...]
        charge = f["charge"][...]
        dt = f["dt"][...]

    # MPI arange data
    for i in range(time.shape[0]):
        det.cmp_emf(position_cur[i], position_prev[i], beta_cur[i], beta_prev[i], time[i], charge, dt)
    det.reduce()

    # Get Data
    field3d = det.get_emf()
    screen = det.get_screen_potisions()
    data = np.array(field3d.to_memoryview(), copy = False)

    lambda_L = 0.8
    wavelength = 2 * pi
    c = 1
    um = wavelength/lambda_L
    fs = 0.3 * um/c
    vmax = 1e-9 #np.max(np.abs(data[:, 0, :]))
    vmin = -vmax
    plt.figure(figsize=[4, 3])
    plt.pcolormesh(axis_phi*180/3.14, axis_time/fs, 
                   data[:, 0, :], 
                   vmax = vmax, vmin = vmin,
                   cmap = "seismic")
    plt.colorbar()
    plt.xlabel("phi [degree]")
    plt.ylabel("time [fs]")
    plt.tight_layout()
    plt.savefig("test.jpg")

    plt.figure(figsize=[4, 3])
    plt.plot(data[:, 0, 50])
    plt.savefig("test2.jpg")

if __name__ == "__main__":
    import sys
    filename = "./data.h5"
    args = sys.argv[1:]
    if len(args) >= 1:
        filename = args[0]
    main(filename)