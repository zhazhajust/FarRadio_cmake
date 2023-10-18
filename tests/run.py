#!/usr/bin/env python3
import os
import yaml
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pyfaradio.faradio import SpheDetector, FaradioMPI

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
    faradio_mpi: FaradioMPI = det.get_mpi()
    det.set_approx(screen["if_approx"])
    with h5py.File(filename, "r") as f:
        time = f["time"][...]
    shape = time.shape[0]
    rank = faradio_mpi.rank()
    size = faradio_mpi.size()
    idxl, idxr = divide_chunk(rank, size, shape)

    with h5py.File(filename, "r") as f:
        #print(rank, idxl, idxr)
        position_cur = f["position_cur"][idxl:idxr]
        position_prev = f["position_prev"][idxl:idxr]
        beta_cur = f["beta_cur"][idxl:idxr]
        beta_prev = f["beta_prev"][idxl:idxr]
        time = f["time"][idxl:idxr]
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

    # Plot Data
    if faradio_mpi.rank() == 0:
        plt.pcolormesh(data[:, :, 0], cmap = "seismic")
        plt.colorbar()
        plt.savefig("test.jpg")

if __name__ == "__main__":
    import sys
    filename = "./data.h5"
    # Print the script name
    print("Script name:", sys.argv[0])
    # Print the command-line arguments
    print("Arguments:", sys.argv[1:])
    args = sys.argv[1:]
    if len(args) >= 1:
        filename = args[0]
    main(filename)