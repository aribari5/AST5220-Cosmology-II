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
    data                 = np.loadtxt(filename)
    x                    = data[:, 0]
    Theta0               = data[:, 1]
    Theta1               = data[:, 2]
    Theta2               = data[:, 3]
    Phi                  = data[:, 4]
    Psi                  = data[:, 5]
    Pi                   = data[:, 6]

    delta_cdm            = data[:, 7]
    delta_b              = data[:, 8]

    v_cdm                = data[:, 9]
    v_b                  = data[:, 10]

    Source_T             = data[:, 11]
    Source_T_j_ell_5     = data[:, 12]
    Source_T_j_ell_50    = data[:, 13]
    Source_T_j_ell_500   = data[:, 14]
    

    return x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500

def plot_delta_gamma_cdm_b():

    k_values     = [0.1, 0.01, 0.001]
    k_val_colors = ['green', 'red', 'blue']

    linestyles   = {'dotted': r'$\delta_\gamma$',
                    '--'    : r'$\delta_\mathrm{CDM}$',
                    'solid' : r'$\delta_b$'}

    
    plt.figure()


    for k in k_values:
        
        file_path = f"./perturbations_k{k}.txt"
        x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500 = load_perturbation_data(file_path)
        delta_gamma = 4*Theta0 

        colour = k_val_colors[k_values.index(k)]
       
        plt.semilogy(x, delta_gamma, ls = 'dotted', color=colour)
        plt.semilogy(x, delta_cdm,   ls = '--',color=colour)
        plt.semilogy(x, delta_b,     ls = 'solid',color=colour, label = f'$k={k}/Mpc$')

    k_val_legend = plt.legend(loc='upper left')
    plt.gca().add_artist(k_val_legend)

    linestyle_legend_handles = [
        plt.Line2D([0], [0], color='black', linestyle=ls, label=label)
        for ls, label in linestyles.items()
    ]
    plt.legend(handles=linestyle_legend_handles, loc='upper left', bbox_to_anchor=(0, 0.85))

    plt.xlabel(r'$x=\ln a$')
    plt.ylabel(r'Perturbation Amplitude')
    plt.ylim(1e-1,1e5)

    plt.show()



if __name__ == "__main__":
    plot_style()
    plot_delta_gamma_cdm_b()