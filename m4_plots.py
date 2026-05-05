import numpy as np
import matplotlib.pyplot as plt
import scipy as sp



def plot_style():

    # Set the style preferences:
    plt.style.use("seaborn-v0_8-darkgrid")
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 12,
        "figure.figsize": (8,6),
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "text.usetex": True,
        "xtick.labelsize": 12,      
        "ytick.labelsize": 12,      
        "xtick.major.size": 6,      
        "xtick.minor.size": 3,      
        "ytick.major.size": 6,      
        "ytick.minor.size": 3,      
        "xtick.major.width": 1.0,   
        "xtick.minor.width": 0.75,  
        "ytick.major.width": 1.0,   
        "ytick.minor.width": 0.75,  
        "xtick.direction": "out",   
        "ytick.direction": "out",   
        "xtick.color": "black",     
        "ytick.color": "black",
    })


def load_power_spectrum_data(filename):
    # Load the data from power_spectrum.txt
    data = np.loadtxt(filename)
    ell  = data[:, 0]
    Cell = data[:, 1]

    return ell, Cell

def plot_Cell():
    ell, Cell = load_power_spectrum_data("cells.txt")

    plt.figure()
    plt.plot(ell, Cell, label=r"$C_\ell^{TT}$", color="blue")
    # plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"Multipole moment $\ell$")
    plt.ylabel(r"Power spectrum $\frac{\ell(\ell+1)}{2\pi}C_\ell^{TT}$")
    plt.title("CMB Temperature Power Spectrum")
    plt.legend()
    plt.tight_layout()
    plt.show()




### calling the plots ###

if __name__ == "__main__":
    plot_style()
    plot_Cell()
