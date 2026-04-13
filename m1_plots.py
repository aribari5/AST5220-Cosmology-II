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


def load_background_data():
    # Load the data from cosmology.txt
    filename = "./cosmology.txt"
    data = np.loadtxt(filename)

    return data

def calculate_x_radiation_matter_equality():
    
    data = load_background_data()

    Omega_r     = data[0,16]
    Omega_nu    = data[0,17]

    Omega_b     = data[0,13]
    Omega_CDM   = data[0,14]

    # From analytical calculaton
    x_rad_mat_eq = np.log((Omega_r + Omega_nu)/(Omega_b + Omega_CDM) )

    print(f"Matter-radiation equality occurs at x = {x_rad_mat_eq:.3f}, z = {np.exp(-x_rad_mat_eq)-1:.3f}")
    

    return x_rad_mat_eq



def calculate_x_matter_dark_energy_equality():

    data = load_background_data()

    Omega_b     = data[0,13]
    Omega_CDM   = data[0,14]
    Omega_Lambda = data[0,15]

    # From analytical calculaton
    x_mat_DE_eq = -(1/3)*np.log((Omega_Lambda)/(Omega_b + Omega_CDM) )

    print(f"Matter-dark energy equality occurs at x = {x_mat_DE_eq:.3f}, z = {np.exp(-x_mat_DE_eq)-1:.3f}")

    return x_mat_DE_eq

def caluclate_x_onset_of_acceleration():
    data = load_background_data()

    Omega_b     = data[0,13]
    Omega_CDM   = data[0,14]
    Omega_Lambda = data[0,15]

    # From analytical calculaton
    x_ons_of_acc = (-1/3)*np.log((2*Omega_Lambda)/(Omega_b+Omega_CDM))

    print(f"Onset of acceleration occurs at x = {x_ons_of_acc:.3f}, z = {np.exp(-x_ons_of_acc)-1:.3f}") 

    return x_ons_of_acc


def load_mcmc_results():
    # Load the results from the MCMC fitting

    filename = "results_supernovafitting.txt"
    data = np.loadtxt(filename, skiprows=500)   # Skip the burnin of the first chains

    # Extract parameters
    chi2    = data[:,0]
    h       = data[:,1]
    Omega_M = data[:,2]
    Omega_K = data[:,3]

    return chi2, h, Omega_M, Omega_K
        
#===========================================#
# Now the functions for the different plots 
#===========================================#

def plot_eta_of_x(filename):
    data       = np.loadtxt(filename, skiprows=1)     # Skip the header
    x          = data[:,0]
    eta         = data[:,1]/(1e6*sp.constants.parsec)    

    plt.figure()

    plt.semilogy(
        x,
        eta,
        label=r"$\eta(x)$",
        color='blue',
    )
    
    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$\eta(x)$ [Mpc]")
    plt.legend()
    plt.xlim(-12, 0)
    plt.ylim(1e0, 1e4)

    plt.tight_layout()
    plt.savefig("figures/eta_of_x.pdf")
    plt.show()

def plot_t_of_x(filename):
    data       = np.loadtxt(filename, skiprows=1)     # Skip the header
    x          = data[:,0]
    t          = data[:,2]/(1e9*365.25*24*3600)

    plt.figure()

    plt.semilogy(
        x,
        t,
        label=r"$t(x)$",
        color='green',
    )

    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$t(x)$ [Gyr]")
    plt.legend()
    # plt.xlim(0, 0)
    # plt.ylim(0,0)

    plt.tight_layout()
    plt.savefig("figures/t_of_x.pdf")
    plt.show()




def plot_luminosity_distance(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    z           = data[:,0]                            # Redshift
    d_L         = data[:,1]                            # in Gpc
    errorbars   = data[:,2]                            # in Gpc

    cosmo_data  = load_background_data()

    x           = cosmo_data[:,0]
    z_model     = np.exp(-x) - 1
    dL_model_SI = cosmo_data[:,11]   # in m!
    dL_model_Gpc= dL_model_SI / (1e9*sp.constants.parsec)  #Gpc

    # Best fit model from MCMC
    mcmc_data = load_mcmc_results()

    chi2        = mcmc_data[0]
    h           = mcmc_data[1]
    Omega_M     = mcmc_data[2]
    Omega_K     = mcmc_data[3]

    best_fit_index = np.argmin(chi2)
    best_fit_h = h[best_fit_index]
    best_fit_Omega_M = Omega_M[best_fit_index]
    best_fit_Omega_K = Omega_K[best_fit_index]
    best_fit_Omega_Lambda = 1 - best_fit_Omega_M - best_fit_Omega_K
    
    





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

    plt.semilogx(
        z_model,
        dL_model_Gpc/z_model,
        label="Fiducial model",
        color="red")
    

    plt.xlabel(r"$z$")
    plt.ylabel(r"$d_L(z)/z$ [Gpc]")
    plt.legend()
    plt.xlim(0, 1.5)
    plt.ylim(3.5,8)


    plt.tight_layout()
    plt.savefig("figures/dL_z_SN.pdf")
    plt.show()
    


def plot_dHpdx_over_Hp(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    x           = data[:,0]                            
    dHpdx_over_Hp = data[:,4]/data[:,3]                       # dHpdx / Hp

    plt.figure()

    plt.plot(
        x,
        dHpdx_over_Hp,
        label=r"$\frac{d\mathcal{H}/dx}{\mathcal{H}}$",
        color='green'
    )
    
    #Vertical lines for matter, radiation and dark energy domination

    x_rad_mat_eq = calculate_x_radiation_matter_equality()
    x_mat_DE_eq = calculate_x_matter_dark_energy_equality()
    
    plt.axvline(x=x_rad_mat_eq, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=x_mat_DE_eq, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   

    #Vertical line for onset of acceleration

    x_onset_of_acc =  caluclate_x_onset_of_acceleration()

    plt.axvline(x=x_onset_of_acc, color='gray', linestyle=':', alpha=0.7, label="Onset of acceleration")   


    # Horizontal lines for the analytical convergence in the different regimes  
    plt.axhline(y=-1, color='orange', linestyle='--', alpha=0.7, label="Radiation dominated era convergence")
    plt.axhline(y=-0.5, color='orange', linestyle='-.', alpha=0.7, label="Matter dominated era convergence")
    plt.axhline(y=1, color='orange', linestyle=':', alpha=0.7, label="DE dominated era convergence")

    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$\frac{d\mathcal{H}/dx}{\mathcal{H}}$")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    # plt.savefig("figures/Hp'_over_Hp.pdf")
    plt.show()

def plot_ddHpddx_over_Hp(filename):
    data       = np.loadtxt(filename, skiprows=1)
    x          = data[:,0]
    ddHpddx_over_Hp = data[:,12]/data[:,3]                       # ddHpddx / Hp
    plt.figure()

    plt.plot(
        x,
        ddHpddx_over_Hp,
        label=r"$\frac{d^2\mathcal{H}/dx^2}{\mathcal{H}}$",
        color='green'
    )
    
    #Vertical lines for matter, radiation and dark energy domination

    x_rad_mat_eq = calculate_x_radiation_matter_equality()
    x_mat_DE_eq = calculate_x_matter_dark_energy_equality()
    
    plt.axvline(x=x_rad_mat_eq, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=x_mat_DE_eq, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   

    #Vertical line for onset of acceleration

    x_onset_of_acc =  caluclate_x_onset_of_acceleration()

    plt.axvline(x=x_onset_of_acc, color='gray', linestyle=':', alpha=0.7, label="Onset of acceleration")   


    # Horizontal lines for the analytical convergence in the different regimes  
    plt.axhline(y=1, color='orange', linestyle='--', alpha=0.7, label="Radiation and DE dominated era convergence")
    plt.axhline(y=0.25, color='orange', linestyle='-.', alpha=0.7, label="Matter dominated era convergence")
    


    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$\frac{d^2\mathcal{H}/dx^2}{\mathcal{H}}$")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.savefig("figures/Hp''_over_Hp.pdf")
    plt.show()





def plot_etaHp_over_c(filename):
    data        = np.loadtxt(filename, skiprows=1)     # Skip the header
    x           = data[:,0]                            # x = -ln(1+z)
    etaHp_over_c = data[:,1]*data[:,3]/sp.constants.c                      # eta*Hp / c

    plt.figure()

    plt.plot(
        x,
        etaHp_over_c,
        label=r"$\frac{\eta \mathcal{H}}{c}$",
        color='purple'
    )
    plt.axhline(y=1, color='black', linestyle='--', alpha=0.7, label=r"Convergence to 1 at early times")
    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$\frac{\eta \mathcal{H}}{c}$")
    plt.xlim(-14.0, 0)
    plt.ylim(0.75,3)

    #Vertical lines for matter, radiation and dark energy domination

    x_rad_mat_eq = calculate_x_radiation_matter_equality()
    x_mat_DE_eq = calculate_x_matter_dark_energy_equality()
    
    plt.axvline(x=x_rad_mat_eq, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=x_mat_DE_eq, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   

    #Vertical line for onset of acceleration

    x_onset_of_acc =  caluclate_x_onset_of_acceleration()

    plt.axvline(x=x_onset_of_acc, color='gray', linestyle=':', alpha=0.7, label="Onset of acceleration")   


    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/etaHp_over_c.pdf")
    plt.show()

def plot_Hp():
    data        = load_background_data()

    Mpc_in_km = 1e3*sp.constants.parsec
    x           = data[:,0]
    Hp = data[:,3] * Mpc_in_km / 100.0                                         
   

    plt.figure()

    plt.semilogy(
        x,
        Hp,
        label=r"$\mathcal{H}(x)$",
        color='black'
    )

    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"$\mathcal{H}(x)$")
    plt.xlim(-12,5)
    plt.ylim(1e-1, 1e3)

    #Vertical lines for matter, radiation and dark energy domination

    x_rad_mat_eq = calculate_x_radiation_matter_equality()
    x_mat_DE_eq = calculate_x_matter_dark_energy_equality()
    
    plt.axvline(x=x_rad_mat_eq, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=x_mat_DE_eq, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   

    #Vertical line for onset of acceleration

    x_onset_of_acc =  caluclate_x_onset_of_acceleration()

    plt.axvline(x=x_onset_of_acc, color='gray', linestyle=':', alpha=0.7, label="Onset of acceleration")   


    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/Hp.pdf")
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
        label=r"$\Omega_\mathrm{Relativistic} = \Omega_\gamma + \Omega_{\nu}$",
        color='red'
    )

    plt.plot(
        x,
        OmegaLambda,
        label=r"$\Omega_\Lambda$",
        color='green'
    )

    plt.plot(
        x,
        Omega_M+Omega_Rel+OmegaLambda,
        label=r"$\Sigma_i\Omega_i$",
        color='black',
        linestyle=':'
    )

    #Vertical lines for matter, radiation and dark energy domination

    x_rad_mat_eq = calculate_x_radiation_matter_equality()
    x_mat_DE_eq = calculate_x_matter_dark_energy_equality()
    
    plt.axvline(x=x_rad_mat_eq, color='gray', linestyle='--', alpha=0.7, label="Radiation dominated era stops")    
    plt.axvline(x=x_mat_DE_eq, color='gray', linestyle='-.', alpha=0.7, label="DE dominated era starts")   

    #Vertical line for onset of acceleration

    x_onset_of_acc =  caluclate_x_onset_of_acceleration()

    plt.axvline(x=x_onset_of_acc, color='gray', linestyle=':', alpha=0.7, label="Onset of acceleration")   

         
    # Can also print abs(sum of densities - 1) as a separate sanity check.
    plt.axhline(y=1, color='gray', linestyle=':', alpha=0.4, label="Total expected density = 1")

    plt.xlabel(r"$x=\ln a$")
    plt.ylabel(r"Density parameters")
    plt.legend()
    #plt.xlim(-2.5, 0)

    plt.tight_layout()
    plt.savefig("figures/densities.pdf")
    plt.show()


def plot_mcmc_scatterplot():

    data        = load_mcmc_results()
    chi2        = data[0]
    h           = data[1]
    Omega_M     = data[2]
    Omega_K     = data[3]

    # min chi2 and corresponding parameters
    min_chi2 = np.min(chi2)
    arg_min_chi2 = np.argmin(chi2)

    best_fit_h = h[arg_min_chi2]
    best_fit_Omega_M = Omega_M[arg_min_chi2]
    best_fit_Omega_K = Omega_K[arg_min_chi2]

    # Calculate Omega_Lambda from the flatness condition
    Omega_Lambda = 1 - Omega_M - Omega_K


    # Empty lists get appended with the parameters that are within the 1-sigma
    one_sigma_Omega_M = []
    one_sigma_Omega_Lambda = []
    # Empty lists get appended with the parameters that are within the 2-sigma (i.e. 2sig-1sig region)
    two_sigma_Omega_M = []
    two_sigma_Omega_Lambda = []

    for i in range(len(chi2)):
        if chi2[i] <= min_chi2 + 3.53: 
            one_sigma_Omega_M.append(Omega_M[i])
            one_sigma_Omega_Lambda.append(Omega_Lambda[i])
   
        elif min_chi2 + 3.53 < chi2[i] <= min_chi2 + 8.02: 
            two_sigma_Omega_M.append(Omega_M[i])
            two_sigma_Omega_Lambda.append(Omega_Lambda[i])
    


    

    # Scatterplot in the Omega_M - Omega_Lambda plane
    plt.figure()
    plt.scatter(two_sigma_Omega_M, two_sigma_Omega_Lambda, color='purple', alpha=0.5, label=r"2-$\sigma$ region")
    plt.scatter(one_sigma_Omega_M, one_sigma_Omega_Lambda, color='blue', alpha=0.5, label=r"1-$\sigma$ region")

    # Add the line corresponding to a flat universe
    Omega_M_flat = np.linspace(0, 1, 100)
    Omega_Lambda_flat = 1 - Omega_M_flat
    plt.plot(Omega_M_flat, Omega_Lambda_flat, color='black', linestyle='--', label="Flat Universe")


    # Add the best fit as a separate square
    plt.scatter(best_fit_Omega_M, 1 - best_fit_Omega_M - best_fit_Omega_K, color='red', label="Best fit from Planck", zorder=5)

    plt.xlabel(r"$\Omega_M$")
    plt.ylabel(r"$\Omega_\Lambda$")
    plt.legend()
    plt.xlim(0, 1)
    plt.ylim(0, 1.5)
    plt.tight_layout()
    plt.savefig("figures/MCMC_scatterplot.pdf")
    plt.show()

def plot_mcmc_Omega_Lambda_posterior():

    data        = load_mcmc_results()
    chi2        = data[0]
    h           = data[1]
    Omega_M     = data[2]
    Omega_K     = data[3]

    # min chi2 and corresponding parameters
    min_chi2 = np.min(chi2)
    arg_min_chi2 = np.argmin(chi2)

    best_fit_h = h[arg_min_chi2]
    best_fit_Omega_M = Omega_M[arg_min_chi2]
    best_fit_Omega_K = Omega_K[arg_min_chi2]

    best_fit_Omega_Lambda = 1 - best_fit_Omega_M - best_fit_Omega_K

    # Calculate Omega_Lambda from the flatness condition
    Omega_Lambda = 1 - Omega_M - Omega_K

    # Values for Gaussian approximation of the posterior 
    mean_Omega_Lambda = np.mean(Omega_Lambda)
    std_Omega_Lambda = np.std(Omega_Lambda)

    #Fit a Gaussian to the posterior
    x  = np.linspace(min(Omega_Lambda), max(Omega_Lambda), 200)

    plt.figure()

    _,bins,_ =  plt.hist(Omega_Lambda, bins=45, color='blue', alpha=0.7)
    bin_width = bins[1] - bins[0]

    gaussian_fit = sp.stats.norm.pdf(x, mean_Omega_Lambda, std_Omega_Lambda)*len(Omega_Lambda)*(bin_width)  # Scale the Gaussian to match the histogram


    plt.axvline(mean_Omega_Lambda, color='red', linestyle='--', label=f"Mean")
    plt.axvline(mean_Omega_Lambda - std_Omega_Lambda, color='red', linestyle=':', label=r"$\pm1$-$\sigma$")
    plt.axvline(mean_Omega_Lambda + std_Omega_Lambda, color='red', linestyle=':')
    plt.axvline(best_fit_Omega_Lambda, color='black', linestyle='--', label=f"Best fit from Planck")
    plt.plot(x, gaussian_fit, color='orange', linestyle='-', label="Gaussian fit to posterior")

    plt.xlabel(r"$\Omega_\Lambda$")
    plt.ylabel(r"Frequency")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/MCMC_Omega_Lambda_posterior.pdf")
    plt.show()


def plot_mcmc_H0():
    data        = load_mcmc_results()
    chi2        = data[0]
    h           = data[1]
    Omega_M     = data[2]
    Omega_K     = data[3]

    # min chi2 and corresponding parameters
    min_chi2 = np.min(chi2)
    arg_min_chi2 = np.argmin(chi2)

    best_fit_h = h[arg_min_chi2]
    best_fit_Omega_M = Omega_M[arg_min_chi2]
    best_fit_Omega_K = Omega_K[arg_min_chi2]

    best_fit_Omega_Lambda = 1 - best_fit_Omega_M - best_fit_Omega_K

    # Calculate H0 from h
    H0 = best_fit_h * 100.0

    # Values for Gaussian approximation of the posterior
    mean_h = np.mean(h)
    std_h = np.std(h)

    x  = np.linspace(min(h), max(h), 200)

    plt.figure()

    _,bins,_ =  plt.hist(h, bins=45, color='blue', alpha=0.7)
    bin_width = bins[1] - bins[0]

    gaussian_fit = sp.stats.norm.pdf(x, mean_h, std_h)*len(h)*(bin_width)  # Scale the Gaussian to match the histogram

    plt.axvline(mean_h, color='red', linestyle='--', label=f"Mean={100.0*mean_h:.1f} km/s/Mpc")
    plt.axvline(mean_h - std_h, color='red', linestyle=':', label=r"$\pm1$-$\sigma$")
    plt.axvline(mean_h + std_h, color='red', linestyle=':')
    plt.axvline(best_fit_h, color='black', linestyle='--', label=f"Best fit from Planck={H0:.1f} km/s/Mpc")
    plt.plot(x, gaussian_fit, color='orange', linestyle='-', label="Gaussian fit to posterior")

    plt.xlabel(r"$h$")
    plt.ylabel(r"Frequency")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/MCMC_H0_posterior.pdf")
    plt.show()





# Setting the style and calling the plots
if __name__ == "__main__":
    plot_style()

    calculate_x_radiation_matter_equality()
    calculate_x_matter_dark_energy_equality()
    caluclate_x_onset_of_acceleration()

    # plot_eta_of_x("cosmology.txt") 
    # plot_t_of_x("cosmology.txt")
    # plot_luminosity_distance("data/supernovadata.txt") 
    # plot_dHpdx_over_Hp("cosmology.txt") 
    # plot_ddHpddx_over_Hp("cosmology.txt")  
    # plot_etaHp_over_c("cosmology.txt") 
    # plot_Hp()  
    # plot_densities()
    # plot_mcmc_scatterplot()
    # plot_mcmc_Omega_Lambda_posterior()
    # plot_mcmc_H0()


