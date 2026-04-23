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

def load_perturbation_data(filename):
    # Load the data from perturbations.txt
    filename = "./perturbations.txt"
    data = np.loadtxt(filename)

    return data

def plot_delta_gamma_cdm_b():

    k_values = [0.1, 0.01, 0.001]
    for k in k_values:

        data            = load_perturbation_data(f'perturbations_k{k}.txt')
        x               = np.linspace(-18.0, 0.0, len(data[:,0]))
        Theta0          = data[:,5]
        delta_gamma     = 4*Theta0
        delta_cdm       = data[:,0]
        delta_b         = data[:,1]

        plt.figure()
        plt.semilogy(x, delta_gamma, label=r'$\delta_\gamma$', color='purple')
        plt.semilogy(x, delta_cdm, label=r'$\delta_{cdm}$', color='black')
        plt.semilogy(x, delta_b, label=r'$\delta_b$', color='red')
        plt.show()
### blablabla fix later



if __name__ == "__main__":
    plot_style()
    plot_delta_gamma_cdm_b()