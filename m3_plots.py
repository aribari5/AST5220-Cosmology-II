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

def plot_v_gamma_cdm_b():

    k_values     = [0.1, 0.01, 0.001]
    k_val_colors = ['green', 'red', 'blue']

    linestyles   = {'dotted': r'$v_\gamma$',
                    '--'    : r'$v_\mathrm{CDM}$',
                    'solid' : r'$v_b$'}

    
    plt.figure()


    for k in k_values:
        
        file_path = f"./perturbations_k{k}.txt"
        x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500 = load_perturbation_data(file_path)
        v_gamma = -3.0*Theta1

        colour = k_val_colors[k_values.index(k)]
       
        plt.semilogy(x, v_gamma, ls = 'dotted', color=colour)
        plt.semilogy(x, v_cdm,   ls = '--',color=colour)
        plt.semilogy(x, v_b,     ls = 'solid',color=colour, label = f'$k={k}/Mpc$')

    k_val_legend = plt.legend(loc='upper left')
    plt.gca().add_artist(k_val_legend)

    linestyle_legend_handles = [
        plt.Line2D([0], [0], color='black', linestyle=ls, label=label)
        for ls, label in linestyles.items()
    ]
    plt.legend(handles=linestyle_legend_handles, loc='upper left', bbox_to_anchor=(0, 0.85))

    plt.xlabel(r'$x=\ln a$')
    plt.ylabel(r'Perturbation Amplitude')
    plt.ylim(1e-6,1e2)

    plt.show()


def plot_Theta0_Theta1():

    k_values     = [0.1, 0.01, 0.001]
    k_val_colors = ['green', 'red', 'blue']

    linestyles   = {'solid': r'$\Theta_0$',
                    'solid' : r'$\Theta_1$'}

    
    
    fig, (ax1, ax2) = plt.subplots(1, 2)

    for k in k_values:
        
        file_path = f"./perturbations_k{k}.txt"
        x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500 = load_perturbation_data(file_path)
        

        colour = k_val_colors[k_values.index(k)]
       
        ax1.plot(x, Theta0, ls = 'solid', color=colour, label = f'$k={k}/Mpc$')
        ax1.set_xlabel(r'$x=\ln a$')
        ax1.set_ylabel(r'$\Theta_0$')
        ax1.legend(loc='upper left')

        ax1.set_ylim(-0.5,1)

        ax2.plot(x, Theta1, ls = 'solid', color=colour, label = f'$k={k}/Mpc$')
        ax2.set_xlabel(r'$x=\ln a$')
        ax2.set_ylabel(r'$\Theta_1$')
        ax2.legend(loc='upper left')

        ax2.set_ylim(-0.4,0.5)





    k_val_legend = plt.legend(loc='upper left')
    plt.gca().add_artist(k_val_legend)

    # linestyle_legend_handles = [
    #     plt.Line2D([0], [0], color='black', linestyle=ls, label=label)
    #     for ls, label in linestyles.items()
    # ]
    # plt.legend(handles=linestyle_legend_handles, loc='upper left', bbox_to_anchor=(0, 0.85))

    plt.tight_layout()
    plt.show()

def plot_Phi():
    k_values     = [0.1, 0.01, 0.001]
    k_val_colors = ['green', 'red', 'blue']

    plt.figure()


    for k in k_values:
        
        file_path = f"./perturbations_k{k}.txt"
        x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500 = load_perturbation_data(file_path)
        

        colour = k_val_colors[k_values.index(k)]
       
        plt.plot(x, Phi, ls = 'solid',color=colour, label = f'$k={k}/Mpc$')

    k_val_legend = plt.legend(loc='lower left')
    plt.gca().add_artist(k_val_legend)

    linestyle_legend_handles = [
        plt.Line2D([0], [0], color='black', linestyle='solid', label=r"$\Phi$")
    ]
    plt.legend(handles=linestyle_legend_handles, loc='upper left', bbox_to_anchor=(0, 0.85))

    plt.xlabel(r'$x=\ln a$')
    plt.ylabel(r'Perturbation Amplitude')
    # plt.ylim(0.1,0.8)

    plt.show()


def plot_Phi_PhiplusPsi():

    k_values     = [0.1, 0.01, 0.001]
    k_val_colors = ['green', 'red', 'blue']

    linestyles   = {'solid': r'$\Theta_0$',
                    'solid' : r'$\Theta_1$'}

    
    
    fig, (ax1, ax2) = plt.subplots(2, 1)

    for k in k_values:
        
        file_path = f"./perturbations_k{k}.txt"
        x,Theta0,Theta1,Theta2,Phi,Psi,Pi,delta_cdm,delta_b,v_cdm,v_b,Source_T,Source_T_j_ell_5,Source_T_j_ell_50,Source_T_j_ell_500 = load_perturbation_data(file_path)
        
        Phi_plus_Psi = Phi + Psi

        colour = k_val_colors[k_values.index(k)]
       
        ax1.plot(x, Phi, ls = 'solid', color=colour, label = f'$k={k}/Mpc$')
        ax1.set_xlabel(r'$x=\ln a$')
        ax1.set_ylabel(r'$\Phi$')
        ax1.legend(loc='upper left')

        ax1.set_ylim(-0.5,1)

        ax2.plot(x, Phi_plus_Psi, ls = 'solid', color=colour, label = f'$k={k}/Mpc$')
        ax2.set_xlabel(r'$x=\ln a$')
        ax2.set_ylabel(r'$\Phi+\Psi$')
        ax2.legend(loc='upper left')

        ax2.set_ylim(-0.4,0.5)





    k_val_legend = plt.legend(loc='upper left')
    plt.gca().add_artist(k_val_legend)

    # linestyle_legend_handles = [
    #     plt.Line2D([0], [0], color='black', linestyle=ls, label=label)
    #     for ls, label in linestyles.items()
    # ]
    # plt.legend(handles=linestyle_legend_handles, loc='upper left', bbox_to_anchor=(0, 0.85))

    plt.tight_layout()
    plt.show()





if __name__ == "__main__":
    plot_style()
    plot_delta_gamma_cdm_b()
    plot_Theta0_Theta1()
    plot_v_gamma_cdm_b()
    plot_Phi()
    plot_Phi_PhiplusPsi()