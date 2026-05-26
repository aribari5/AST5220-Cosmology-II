#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 0.0;     //3.046;
  double TCMB        = 2.7255;

  // Best fit parameters (From MCMC with background parameters)
  double h_bestfit            = 0.702;
  double Omega_M_bestfit      = 0.255;
  double Omega_K_bestfit      = 0.079;
  double Omega_Lambda_bestfit = 0.666;
  double Omega_CDM_bestfit    = Omega_M_bestfit - OmegaB;




  // Recombination parameters
  double Yp          = 0.0;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();

  // For best fit parameters
  BackgroundCosmology cosmo_bestfit(h_bestfit, OmegaB, Omega_CDM_bestfit, Omega_K_bestfit, Neff, TCMB); 
  // cosmo_bestfit.solve();
  // cosmo_bestfit.info(); 

  
  // // Output background evolution quantities

  // cosmo.output("cosmology.txt");
  // cosmo_bestfit.output("cosmology_bestfit.txt");


  // // Do the supernova fits. Uncomment when you are ready to run this
  // // Make sure you read the comments on the top of src/SupernovaFitting.h
  // // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");  // Done:)
  

  // // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  // rec.output("recombination.txt");
  
  // Remove when module is completed
  // return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double k_small = 0.001 / Constants.Mpc;
  double k_intermdeiate = 0.01 / Constants.Mpc;
  double k_large = 0.1 / Constants.Mpc;

  pert.output(k_small, "perturbations_k0.001.txt");
  pert.output(k_intermdeiate, "perturbations_k0.01.txt");
  pert.output(k_large, "perturbations_k0.1.txt");


  // Run w different terms on off new name files

  // Remove when module is completed
  // return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  power.output_pk("powerspectrum.txt");
  power.output_Theta_ells("Theta_ell_of_k.txt", power.get_thetaT_ell_of_k_spline());
  
  // Remove when module is completed
  // return 0;

  Utils::EndTiming("Everything");
}

