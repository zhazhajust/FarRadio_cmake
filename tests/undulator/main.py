#/usr/bin/env python3
from pyfaradio import FaradioSim

if __name__ == "__main__":
    filename = "./data.h5"
    savefile = "./save_data.h5"
    configfile = "./screen.yaml"
    simulation = FaradioSim(filename, savefile, configfile)
    simulation()
