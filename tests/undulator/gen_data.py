#!/usr/bin/env python3
import h5py
import numpy as np

pi = np.pi
lambda_L = 0.8
wavelength = 2 * pi
c = 1
um = wavelength/lambda_L
fs = 0.3 * um/c

v = 0.995
T = 0.8 * um/(c - v)
#T0 = 0.8*um/c

def get_position(time, x0 = 0, t0 = 0):
    y = np.sin((time - t0) * 2 * pi/T) * 0.05
    #z = np.cos(time * 2 * pi/T) * 0.05
    x = (time - t0) * v + x0
    z = np.zeros_like(x)
    return np.array([x, y, z]).T.reshape(-1, 1, 3)

def generate_data():
    time = np.arange(0, 10000, 5)
    dt = time[1] - time[0]

    position = get_position(time)
    delta = 0.1
    # for i in range(1, 20):
    #     x0 = -i * T * delta
    #     t0 = 0 #x0
    #     position = np.concatenate((position, get_position(time, x0, t0)), axis=1)
    beta = np.diff(position, axis=0)/dt
    position = position[:-1]

    time = time[1:-1]
    beta_prev = beta[:-1, :]
    beta_cur = beta[1:, :]
    position_prev = position[:-1, :]
    position_cur = position[1:, :]
    charge = np.ones([time.shape[0], beta_cur.shape[1]], dtype=np.float64)

    with h5py.File("data.h5", "w") as f:
        f.create_dataset("time", data=time)
        f.create_dataset("position_cur", data=position_cur)
        f.create_dataset("position_prev", data=position_prev)
        f.create_dataset("beta_cur", data=beta_cur)
        f.create_dataset("beta_prev", data=beta_prev)
        f.create_dataset("charge", data=charge)
        f.create_dataset("dt", data=dt)

    return

if __name__ == "__main__":
    generate_data()