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
    filename = "/home/yujq/users/caijie/PostProcessing/smilei/new_data.h5" #"data.h5"
    time, position_cur = read(filename)
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(figsize=[4, 3])
    for i in range(0, position_cur.shape[1], position_cur.shape[1]//20):
    # for i in [0, -5]:
        line = plot_lines(axs, position_cur[:, i, 0]/um, position_cur[:, i, 1]/um,
                    time/fs, lw=0.5, clim = [0, time.max()/fs], cmap = "jet")
    plt.xlim([-10, 5000])
    plt.ylim([-10, 10])
    plt.colorbar(line)
    plt.xlabel("x [um]")
    plt.ylabel("y [um]")

    plt.tight_layout()
    plt.savefig("position_cur.jpg", dpi=300)

if __name__ == "__main__":
    setFont()
    main()