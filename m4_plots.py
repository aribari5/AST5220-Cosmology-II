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


def load_cell_data(filename):
    # Load the data from cells.txt

    data = np.loadtxt(filename)
    ell  = data[:, 0]
    Cell = data[:, 1]

    return ell, Cell

def load_power_spectrum_data(filename):
    # Load the data from powerspectrum.txt

    data = np.loadtxt(filename)
    k    = data[:, 0]
    Pk   = data[:, 1]

    return k, Pk


def load_Theta_ell_data(filename):
    # Load the data from Theta_ell_of_k.txt
    # k 0

    data            = np.loadtxt(filename)
    k               = data[0,1:]
    Theta_ell_of_k  = data[1:,:]

    return k, Theta_ell_of_k





def plot_Cell():
    ell, Cell = load_cell_data("cells.txt")

    plt.figure()
    plt.plot(ell, Cell, label=r"$C_\ell^{TT}$", color="blue")
    plt.yscale("log")
    plt.xscale("log")

    plt.xlabel(r"Multipole $\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}C_\ell^{TT}\left[\mu K^2\right]$")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_Pk():
    k, Pk = load_power_spectrum_data("powerspectrum.txt")

    plt.figure()
    plt.plot(k, Pk, label=r"$P(k)$", color="red")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$k [h/Mpc^{-1}$]")
    plt.ylabel(r"$P(0,k) [Mpc^3/h^3]$")
    plt.legend()
    plt.tight_layout()
    plt.show()

"""
def plot_Source_function():

    data = np.loadtxt("perturbations_k.01.txt")

    x                                       = data[:, 0]

    Source_func_full                        = data[:,11]
    Source_func_j_ell_5                     = data[:,12]
    Source_func_j_ell_50                    = data[:,13]
    Source_func_j_ell_500                   = data[:,14]

    Source_func_SW_contribution             = data[:,15]
    Source_func_ISW_contribution            = data[:,16]
    Source_func_Doppler_contribution        = data[:,17]
    Source_func_Polarization_contribution   = data[:,18]
    
    plt.figure()

    plt.plot(x, Source_func_full, label=r"$S(k,x)$",ls='solid',color="black")
    plt.plot(x, Source_func_SW_contribution, label=r"SW",ls='dashed',color="blue")
    plt.plot(x, Source_func_ISW_contribution, label=r"ISW",ls='dashed',color="red")
    plt.plot(x, Source_func_Doppler_contribution, label=r"Doppler",ls='dashed',color="green")
    plt.plot(x, Source_func_Polarization_contribution, label=r"Polarization",ls='dashed',color="orange")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$x\?$")
    plt.title(r"Source Function Contributions for $k=0.01/Mpc$")
    plt.legend()
    plt.tight_layout()
    plt.show()
"""

def plot_Theta_ells():

    ells = [15, 100, 500, 1000, 1500, 2000] 
    k, Theta_ell_of_k = load_Theta_ell_data("Theta_ell_of_k.txt")

    ell_array = np.array(len (Theta_ell_of_k[0,:]))

    colors = ["BF1A2F", "F0A202", "98CE00", "16E0BD", "454E9E", "F98404"]

    
    for ell in ells:
        plt.figure()

        normfactor = np.sqrt(ell*(ell+1))
        plt.plot(ell_array, normfactor*Theta_ell_of_k[ell,:], label=rf"$\ell={ell}$", color=colors[ells.index(ell)])
        # plt.xscale("log")
        # plt.yscale("log")
        plt.xlabel(r"$k [h/Mpc^{-1}$]")
        plt.ylabel(r"$\Theta_\ell(k)$")
    
        plt.legend()
    plt.tight_layout()
    plt.show()




### calling the plots ###

if __name__ == "__main__":
    plot_style()
    # plot_Cell()
    # plot_Pk()
    plot_Theta_ells()
    # plot_Source_function()