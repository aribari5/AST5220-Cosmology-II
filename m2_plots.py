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

def load_recombination_data():
    # Load the data from recombination.txt
    filename = "./recombination.txt"
    data = np.loadtxt(filename)

    return data

def plot_Xe():
    data        = load_recombination_data()
    x           = data[:, 0]
    Xe          = data[:, 1]

    plt.figure()
    plt.semilogy(x, Xe, label=r'$X_e(x)$', color='green')

    plt.xlabel(r'$x=\ln a$')
    plt.ylabel(r'Ionization Fraction $X_e$')
    plt.legend()
    plt.xlim(-12,0)
    plt.ylim(-1e4,1e1)
    plt.show()





def plot_optical_depth_taus():

    data        = load_recombination_data()
    x           = data[:, 0]
    tau         = data[:, 3]
    tau_deriv   = data[:, 4]
    tau_deriv2  = data[:, 5]


    plt.figure()
    plt.semilogy(x, tau, label=r'$\tau(x)$', color='blue')
    plt.semilogy(x, -tau_deriv, label=r"$-\tau'(x)$", color='green')
    plt.semilogy(x, tau_deriv2, label=r"$\tau''(x)$", color='red')

    plt.xlabel(r'$x=\ln a$')
    plt.ylabel(r'Optical Depth $\tau$')
    plt.legend()
    plt.ylim(-1e8, 1e8)
    plt.xlim(-12,0)
    plt.show()


### calling the plots ###

if __name__ == "__main__":
    plot_style()
    plot_optical_depth_taus()
    plot_Xe()