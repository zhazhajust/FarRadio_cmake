import matplotlib.pyplot as plt
import h5py
import numpy as np
from pypic.constants import Constants

savefile = "./save_data.h5"

with h5py.File(savefile, "r") as f:
    data = f["data"][...]
    axis_time = f["axis_time"][...]
    axis_phi = f["axis_phi"][...]
    axis_theta = f["axis_theta"][...]

#Constant
lambda_L = 0.8
wavelength = 2 * np.pi
c = 1
um = wavelength/lambda_L
fs = 0.3 * um/c

vmax = 1e-9 #np.max(np.abs(data[:, 0, :]))
vmin = -vmax
# plt.figure(figsize=[4, 3])
# plt.pcolormesh(axis_phi*180/3.14, axis_time/fs, data[:, 0, :], 
#             # vmax = vmax, vmin = vmin,
#             cmap = "seismic")
# plt.colorbar()
# plt.savefig("test1.jpg")

plt.figure(figsize=[4, 3])
plt.pcolormesh(axis_theta*180/3.14, axis_time[::20]/fs, data[::20, :, 0], 
            # vmax = vmax, vmin = vmin,
            cmap = "seismic")
plt.colorbar()
plt.savefig("./test_3.jpg")

plt.figure(figsize=[4, 3])
plt.plot(axis_time/fs, data[:, 10, 0])
plt.savefig("./test_4.jpg")

# plt.figure(figsize=[4, 3])
# plt.plot(axis_time/fs, data[:, 0, 10])
# plt.savefig("test2.jpg")
