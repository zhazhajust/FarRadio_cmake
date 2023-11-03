import h5py
from pypic.font import setFont
from pypic.constants import Constants
from pypic.plot import plot_lines
from pypic.colorbar import getCmap

constants = Constants(0.8)
um = constants.um
fs = constants.fs

def read(filename):
    with h5py.File(filename, "r") as f:
        position_cur = f["position_cur"][...]
        time = f["time"][...]
    return time, position_cur

def main():
    time, position_cur = read("data.h5")
    import matplotlib.pyplot as plt

    xmax, ymax = 0, 0
    xmin, ymin = 2<<8, 2<<8

    fig, axs = plt.subplots(figsize=[4, 3])
    for i in range(0, position_cur.shape[1], 4):
        line = plot_lines(axs, position_cur[:, i, 0]/um, position_cur[:, i, 1]/um,
                    time/fs, lw=0.5, clim = [0, time.max()/fs], cmap = "jet")
        xmax = max(xmax, position_cur[:, i, 0].max()/um)
        xmin = min(xmin, position_cur[:, i, 0].min()/um)
        ymax = max(ymax, position_cur[:, i, 1].max()/um)
        ymin = min(ymin, position_cur[:, i, 1].min()/um)
    plt.xlim([1.2*xmin, 1.2*xmax])
    plt.ylim([1.2*ymin, 1.2*ymax])
    plt.colorbar(line)
    plt.xlabel("x [um]")
    plt.ylabel("y [um]")

    plt.tight_layout()
    plt.savefig("position_cur.jpg", dpi=300)

if __name__ == "__main__":
    setFont()
    main()