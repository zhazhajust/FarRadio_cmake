import h5py
import numpy as np
pi = np.pi
lambda_L = 0.8
wavelength = 2 * pi
c = 1
um = wavelength/lambda_L
fs = 0.3 * um/c

v = 0.995

T = 0.8*um/c
print(T)
def get_position(time):
    z = np.sin(time * 2 * pi/2)# * 0.01
    y = np.cos(time * 2 * pi/2)# * 0.01
    x = time * v
    return np.array([x, y, z]).T

def generate_data():
    time = np.arange(0, 6000, 0.2)
    dt = time[1] - time[0]
    charge = 1

    position = get_position(time)
    beta = np.diff(position, axis=0)/dt
    position = position[:-1]

    time = time[1:-1]
    beta_prev = beta[:-1, :]
    beta_cur = beta[1:, :]
    position_prev = position[:-1, :]
    position_cur = position[1:, :]

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