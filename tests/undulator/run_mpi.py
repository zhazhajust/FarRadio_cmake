#!/usr/bin/env python3
import os
import yaml
import h5py
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from pyfaradio.faradio import SpheDetector, FaradioMPI

def load_screen(filename):
    with open(filename, "r") as f:
        data = yaml.load(f, Loader=yaml.FullLoader)
    return data

def divide_chunk(rank, size, shape):
    chunksize = shape // size
    return rank*chunksize, min(shape, (rank + 1)*chunksize)

def main(filename: str, savefile: str):
    screen = load_screen("/home/yujq/users/caijie/PostProcessing/FarRadio_cmake/tests/screen.yaml")
    os.environ['OMP_NUM_THREADS'] = str(screen["omp_num_threads"])
    det = SpheDetector(screen["dmin"], screen["dmax"], screen["nf"])
    det.set_approx(screen["if_approx"])
    axis_time = np.linspace(screen["dmin"][0], screen["dmax"][0], screen["nf"][0])
    axis_phi = np.linspace(screen["dmin"][2], screen["dmax"][2], screen["nf"][2])
    axis_theta = np.linspace(screen["dmin"][1], screen["dmax"][1], screen["nf"][1])

    # MPI init
    faradio_mpi: FaradioMPI = det.get_mpi()
    if faradio_mpi.rank() == 0:
        # Print the script name
        print("Script name:", sys.argv[0])
        # Print the command-line arguments
        print("Arguments:", sys.argv[1:])
        # Save filename
        print("savefile:", sys.argv[2:])

    # Read Data
    with h5py.File(filename, "r") as f:
        time = f["time"][...]
    shape = time.shape[0]
    rank = faradio_mpi.rank()
    size = faradio_mpi.size()
    idxl, idxr = divide_chunk(rank, size, shape)

    with h5py.File(filename, "r") as f:
        position_cur = f["position_cur"][idxl:idxr]
        position_prev = f["position_prev"][idxl:idxr]
        beta_cur = f["beta_cur"][idxl:idxr]
        beta_prev = f["beta_prev"][idxl:idxr]
        time = f["time"][idxl:idxr]
        # charge = f["charge"][...]
        charge = f["charge"][idxl:idxr]
        dt = f["dt"][...]

    # MPI arange data
    for i in range(time.shape[0]):
        # det.cmp_emf(position_cur[i], position_prev[i], beta_cur[i], beta_prev[i], time[i], charge, dt)
        det.cmp_emf(position_cur[i], position_prev[i], beta_cur[i], beta_prev[i], time[i], charge[i], dt)
    det.reduce()

    # Get Data
    field3d = det.get_emf()
    screen = det.get_screen_potisions()
    data = np.array(field3d.to_memoryview(), copy = False)

    #Constant
    lambda_L = 0.8
    wavelength = 2 * pi
    c = 1
    um = wavelength/lambda_L
    fs = 0.3 * um/c
    
    # Plot Data
    if faradio_mpi.rank() == 0:
        # save data
        with h5py.File(savefile, "w") as f:
            f.create_dataset("data", data = data)
            f.create_dataset("axis_time", data = axis_time)
            f.create_dataset("axis_phi", data = axis_phi)
            f.create_dataset("axis_theta", data = axis_theta)
            # f.create_dataset("screen", data = screen)


if __name__ == "__main__":
    import sys
    # filename = "/home/yujq/users/caijie/PostProcessing/FarRadio_cmake/tests/data.h5"
    filename = "/home/yujq/users/caijie/PostProcessing/smilei/new_data.h5"
    savefile = "/home/yujq/users/caijie/PostProcessing/smilei/save_data.h5"

    args = sys.argv[1:]
    if len(args) >= 1:
        filename = args[0]
    if len(args) >= 2:
        savefile = args[1]

    main(filename, savefile)