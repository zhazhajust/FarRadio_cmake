from pyfaradio.faradio import Particle, Field3D
from pyfaradio.faradio import SpheDetector as _SpheDetector
from numpy import ndarray
from numpy import double

class SpheDetector:
    def __init__(self, dmin: list[double], dmax: list[double], nf: list[double]) -> None:
        self.det = _SpheDetector(dmin, dmax, nf)
        return
    
    def cmp_emf(self, position_cur: ndarray[double], position_prev: ndarray[double], 
                beta_cur: ndarray[double], beta_prev: ndarray[double], 
                time: double, charge: double, particle_nums: double, dt: double):
        self.det.cmp_emf(position_cur, position_prev, beta_cur, beta_prev, time, charge, particle_nums, dt)
        return
    