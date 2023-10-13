import numpy as np
import matplotlib.pyplot as plt
from pyfaradio.faradio import Particle, SpheDetector, Field3D

dt = 0.1
time = 100
charge = 1
iter_id = 200

x = np.linspace(0, 1, 100)
y = np.cos(x)
z = np.zeros_like(x)
position = np.array([x, y, z]).T
beta_x = np.diff(x)/dt
beta_y = np.diff(y)/dt
beta_z = np.diff(z)/dt
beta = np.array([beta_x, beta_y, beta_z]).T
position = position[:-1, :]
position_prev = position[:-1, :]
position_cur = position[1:, :]
beta_prev = beta[:-1, :]
beta_cur = beta[1:, :]

from numpy import pi

det = SpheDetector([0, 0, 0], [1000, pi, pi], [10, 10, 10])

#det.cmp_emf(tracer)
det.cmp_emf(position_cur, position_prev, beta_cur, beta_prev, time, dt, charge, iter_id)

#field3d = det.get_emf()
screen = det.get_screen_potisions()
#screen = det.get_screen_x()
#data = np.array(field3d.to_memoryview(), copy = False)

#data = np.array(screen, copy = True)

#plt.pcolormesh(data[:, :, 0], cmap = "jet")
#plt.colorbar()
#plt.savefig("test.jpg")