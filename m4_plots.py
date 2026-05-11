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

    file_ells       = data[0, 1:].astype(int)
    k               = data[1:,0]
    Theta_ell_of_k  = data[1:,1:]

    return k,file_ells, Theta_ell_of_k

def load_low_ell_TT_data(filename):
    # Load the data from planck_cell_low.txt

    data                = np.loadtxt(filename,skiprows=1)

    ell                 = data[:, 0]
    Dell_planck         = data[:, 1]
    err_up_planck       = data[:, 2]
    err_down_planck     = data[:, 3]


    return ell, Dell_planck, err_up_planck, err_down_planck



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
    
    # see final output from RecombinationHistory.cpp
    k_eq = 0.0201474           
    plt.figure()

    plt.vlines(k_eq,ymin=1e2,ymax=1e5, colors="green", linestyles="dashed", label=r"$k_{eq}$")

    plt.plot(k, Pk, label=r"$P(k)$", color="red")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$k [h/Mpc^{-1}$]")
    plt.ylabel(r"$P(0,k) [Mpc^3/h^3]$")
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_Cell_contributions():

    # load the full data
    data                = np.loadtxt("cells.txt")

    # contributions
    data_SW             = np.loadtxt("cells_SW.txt")
    data_ISW            = np.loadtxt("cells_ISW.txt")
    data_Doppler        = np.loadtxt("cells_Doppler.txt")
    data_Polarization   = np.loadtxt("cells_Polarization.txt")

    sources = ["SW", "ISW", "Doppler", "Polarization"]

    ell  = data[:, 0]

    Cell_full                       = data[:,1]
    Cell_SW_contribution            = data_SW[:,1]
    Cell_ISW_contribution           = data_ISW[:,1]
    Cell_Doppler_contribution       = data_Doppler[:,1]
    Cell_Polarization_contribution  = data_Polarization[:,1]



    
    plt.figure()

    plt.plot(ell, Cell_full, label=r"$S(k,x)$",ls='solid',color="black")
    plt.plot(ell, Cell_SW_contribution, label=r"SW",ls='dashed',color="orange")
    plt.plot(ell, Cell_ISW_contribution, label=r"ISW",ls='dashed',color="green")
    plt.plot(ell, Cell_Doppler_contribution, label=r"Doppler",ls='dashed',color="blue")
    plt.plot(ell, Cell_Polarization_contribution, label=r"Polarization",ls='dashed',color="red")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(r"$\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}C_\ell^{TT}\left[\mu K^2\right]$")
    plt.title(r"$k=0.01/Mpc$")
    plt.ylim(10**(-2.3),1e4)
    plt.xlim(2,2000)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_Theta_ells():

    ells = [15, 100, 500, 1000, 1500, 2000] 
    k_eta0, file_ells, transfer_ell = load_Theta_ell_data("Theta_ell_of_k.txt")

    
    colors = ["#BF1A2F", "#F0A202", "#F4E285", "#98CE00", "#16E0BD", "#759EB8"]
    plt.figure()

    for i, ell in enumerate(ells):
        #Find the index in the file that matches the ell we want
        try:
            column_idx = np.where(file_ells == ell)[0][0]
        except IndexError:
            print(f"Warning: ell={ell} not found in the data file!")
            continue

        
        normfactor = np.sqrt(ell*(ell+1))
        y_values = normfactor * transfer_ell[:, column_idx]
        
        plt.plot(k_eta0, y_values, label=rf"$\ell={ell}$", color=colors[i])
    
    plt.xlabel(r"$k\eta_0$")
    plt.ylabel(r"$\sqrt{\ell(\ell+1)}\,\Theta_\ell(k)$")
    plt.legend()
    plt.show()

def plot_integrand_Theta_ells():
    ells = [15, 100, 500, 1000, 1500, 2000] 
    k,file_ells, transfer_ell = load_Theta_ell_data("Theta_ell_of_k.txt")


    

    colors = ["#BF1A2F", "#F0A202", "#F4E285", "#98CE00", "#16E0BD", "#759EB8"]

    
    plt.figure()

    for i, ell in enumerate(ells):

        #Find the index in the file that matches the ell we want
        try:
            column_idx = np.where(file_ells == ell)[0][0]
        except IndexError:
            print(f"Warning: ell={ell} not found in the data file!")
            continue

        
        normfactor = ell*(ell+1)
        abs_squared_Theta = np.abs(transfer_ell[:,column_idx])**2

        y_values = normfactor*abs_squared_Theta/k

        
        
        
        plt.plot(k, y_values, label=rf"$\ell={ell}$", color=colors[i])
        
        plt.xlabel(r"$k\eta_0$")
        plt.ylabel(r"$\ell(\ell+1)\,|\Theta_\ell(k)|^2/k$")
    
        plt.legend()
    plt.tight_layout()
    plt.show()


def plot_compare_to_Planck_data():
    
    ell_planck, Dell_planck, err_up_planck, err_down_planck = load_low_ell_TT_data("planck_cell_low.txt")

    Cell_planck = Dell_planck
    
    ell_fiducial, Cell_fiducial = load_cell_data("cells.txt")


    
    
    

    plt.figure()
    plt.semilogx(ell_fiducial, Cell_fiducial, label=r"Fiducial $C_\ell^{TT}$", color="blue")
    plt.errorbar(ell_planck, Cell_planck, yerr=[err_up_planck, err_down_planck], fmt='x', label=r"Planck 2018 Data", color="red", ecolor="gray", capsize=3)

   
    plt.xlabel(r"Multipole $\ell$")
    plt.ylabel(r"$\frac{\ell(\ell+1)}{2\pi}C_\ell^{TT}\left[\mu K^2\right]$")

    plt.xlim(2,10**(3.3))
    plt.legend()
    plt.tight_layout()
    plt.show()

    

    




### calling the plots ###

if __name__ == "__main__":
    plot_style()
    # plot_Cell()
    plot_Cell_contributions()
    # plot_Pk()
    # plot_Theta_ells()
    # plot_integrand_Theta_ells()
    # plot_compare_to_Planck_data()