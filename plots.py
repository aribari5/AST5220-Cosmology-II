import numpy as np
import matplotlib.pyplot as plt

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
        "text.usetex": True
    })


def load_background_data():
    # Load the data from cosmology.txt
    filename = "./cosmology.txt"
    data = np.loadtxt(filename)

    return data


def load_mcmc_results():
    # Load the results from the MCMC fitting

    filename = "results_supernovafitting.txt"
    data = np.loadtxt(filename, skiprows=200)   # Skip the burnin of the first chains

    # Extract parameters
    chi2    = data[:,0]
    h       = data[:,1]
    Omega_M = data[:,2]
    Omega_K = data[:,3]

    return chi2, h, Omega_M, Omega_K
        
#===========================================#
# Now the functions for the different plots 
#===========================================#


def plot_luminosity_distance(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    z           = data[:,0]                            # Redshift
    d_L         = data[:,1]                            # in Gyr
    errorbars   = data[:,2]                            # in Gyr

    cosmo_data = load_background_data()
    x           = cosmo_data[:,0]
    z_model     = np.exp(-x) - 1
    dL_model    = cosmo_data[:,11]   

    # We wish to plot d_L / z (remember the errorbars!):
    dL_over_z = d_L / z
    err_over_z = errorbars / z

    plt.figure()

    plt.errorbar(
        z,
        dL_over_z,
        yerr=err_over_z,
        fmt='o',
        markersize=4,
        capsize=3,
        label="Supernova data",
        color='blue',
    )

    plt.plot(
        z_model,
        dL_model/z_model,
        label="Fiducial model",
        color="red")

    plt.xlabel(r"$z$")
    plt.ylabel(r"$d_L(z)/z$ [Gyr]")
    plt.legend()
    plt.xlim(0, 2.5)


    plt.tight_layout()
    plt.show()

def plot_dHpdx_over_Hp(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    x           = data[:,0]                            # x = -ln(1+z)
    dHpdx_over_Hp = data[:,4]/data[:,3]                       # dHpdx / Hp

    plt.figure()

    plt.plot(
        x,
        dHpdx_over_Hp,
        label=r"$\frac{d\mathcal{H}/dx}{\mathcal{H}}$",
        color='green'
    )
    
    #Vertical lines for matter, radiation and dark energy domination

    plt.axvline(x=-8.0067, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=-0.4054, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   
    #plt.axvline(x=-2, color='gray', linestyle=':', alpha=0.7, label="DE dominated era")      



    plt.xlabel(r"$x$")
    plt.ylabel(r"$\frac{d\mathcal{H}/dx}{\mathcal{H}}$")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.show()

def plot_etaHp_over_c(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    x           = data[:,0]                            # x = -ln(1+z)
    etaHp_over_c = data[:,1]*data[:,3]/(3*10**8)                       # eta*Hp / c

    plt.figure()

    plt.plot(
        x,
        etaHp_over_c,
        label=r"$\frac{\eta \mathcal{H}}{c}$",
        color='purple'
    )

    plt.xlabel(r"$x$")
    plt.ylabel(r"$\frac{\eta \mathcal{H}}{c}$")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.show()

def plot_Hp():
    data        = load_background_data()
    x           = data[:,0]                            # x = -ln(1+z)
    Hp          = data[:,3]                            # Hp

    plt.figure()

    plt.plot(
        x,
        Hp,
        label=r"$\mathcal{H}(x)$",
        color='black'
    )

    plt.xlabel(r"$x$")
    plt.ylabel(r"$\mathcal{H}(x)$")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.show()
    
def plot_densities():
    data        = load_background_data()
    x           = data[:,0]                            # x = -ln(1+z)
    OmegaB      = data[:,5]                            # Baryon density
    OmegaCDM    = data[:,6]                            # CDM density
    OmegaLambda = data[:,7]                            # Dark energy density
    OmegaR      = data[:,8]                            # Radiation density
    OmegaNu     = data[:,9]                            # Neutrino density
    OmegaK      = data[:,10]                           # Curvature density
    Omega_M     = OmegaB + OmegaCDM
    Omega_Rel   = OmegaR + OmegaNu


    plt.figure()

    plt.plot(
        x,
        Omega_M,
        label=r"$\Omega_\mathrm{Matter} = \Omega_b + \Omega_{\mathrm{CDM}}$",
        color='blue'
    )

    plt.plot(
        x,
        Omega_Rel,
        label=r"$\Omega_\mathrm{Relativistic} = \Omega_r + \Omega_{\nu}$",
        color='red'
    )

    plt.plot(
        x,
        OmegaLambda,
        label=r"$\Omega_\Lambda$",
        color='green'
    )

    plt.axvline(x=-8.0067, color='gray', linestyle='--', alpha=0.4, label="Radiation dominated era stops")    
    plt.axvline(x=-0.4054, color='gray', linestyle='-.', alpha=0.4, label="DE dominated era starts")  
    plt.axhline(y=1, color='gray', linestyle=':', alpha=0.4, label="Total density")

    plt.xlabel(r"$x$")
    plt.ylabel(r"Density parameters")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.show()



# Setting the style and calling the plots
if __name__ == "__main__":
    plot_style()
    # plot_luminosity_distance("data/supernovadata.txt") 
    # plot_dHpdx_over_Hp("cosmology.txt") 
    # plot_etaHp_over_c("cosmology.txt") # a bit wrong i think
    # plot_Hp() #wrong
    #plot_densities("cosmology.txt")    #wrong





