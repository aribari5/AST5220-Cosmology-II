#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================

  Vector x_array      = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr       = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_saha_arr  = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector ne_arr       = Utils::linspace(x_start, x_end, npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    
    // Store the result in Xe_saha_arr...
    // ...

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the Saha result in Xe_arr 
      //=============================================================================
      //...
      //...

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!) 
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result. Store in Xe_arr 
      //=============================================================================
      //...
      //...
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x,
  // functions are working
  //=============================================================================
  //...
  //...

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB0     = cosmo->get_OmegaB(0.0);

  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaM      = cosmo->get_OmegaM();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  
  const double rho_crit0    = 3.0*pow(H0_over_h,2)*pow(Constants.c,2)/(8.0*Constants.pi*G);       // Critical density today in kg/m^3
  const double TCMB0       = cosmo->get_TCMB(0.0);                                                   // CMB temperature today in K


  const double n_b         = (1-Yp)*OmegaB0*rho_crit0/(m_H * pow(a,3));      // Number density of baryons at x in 1/m^3
  const double n_H         = (1-Yp)*n_b;                                     // Number density of hydrogen at x in 1/m^3
  const double T_b         = TCMB0/a;                                        // Temperature of baryons at x in K


  const double C           = (1/n_b)*pow( (m_e*T_b)/(2*Constants.pi),3/2 )*exp(-epsilon_0/T_b); // Constant from Saha eq. Should there be a kb in exp?

  //=============================================================================
  // Computing Xe and ne from the Saha equation
  //=============================================================================
  
  // Electron fraction and number density

  // Must add  a check feks if 4/C << sim. 300? then use 1.0 early Universe approx.
  double Xe = ( -C + sqrt( pow(C,2) - 4*C) )/2;
  double ne = Xe*n_H;
  

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters

  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB0     = cosmo->get_OmegaB(0.0);
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaM      = cosmo->get_OmegaM();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  
  const double H           = cosmo->H_of_x(x);
  
  const double rho_crit0   = 3.0*pow(H0_over_h,2)*pow(Constants.c,2)/(8.0*Constants.pi*G);       // Critical density today in kg/m^3
  const double TCMB0       = cosmo->get_TCMB(0.0);                                                   // CMB temperature today in K
  const double T_b         = TCMB0/a;


  // All of the constants needed for RHS of Peebles

  const double n_b          = (1-Yp)*OmegaB0*rho_crit0/(m_H * pow(a,3));
  const double n_H          = (1-Yp)*n_b; 
  const double n_1s         = (1-X_e)*n_H;
  const double lambda_2s1s  = Constants.lambda_2s1s;
  const double lambda_alpha = H*pow(3*epsilon_0,3)/( pow(8*Constants.pi,2)*n_1s );
  const double alpha        = 1/137.0359992;
  const double phi2_of_Tb   = 0.448*log(epsilon_0/T_b);
  const double alpha2_of_Tb = (64*Constants.pi)/(sqrt(27*Constants.pi)) * pow(alpha/m_e,2) * sqrt(epsilon_0/T_b) * phi2_of_Tb;
  const double beta_of_Tb   = alpha2_of_Tb*pow(m_e*T_b/(2*Constants.pi),3/2)*exp(-epsilon_0/T_b);
  const double beta2_of_Tb  = beta_of_Tb*exp((3*epsilon_0)/(4*T_b));
  const double Cr_of_Tb     = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2_of_Tb);

  //=============================================================================
  // RHS of Peebles ODE for dXedx
  //=============================================================================
  
  dXedx[0] =  Cr_of_Tb/H* ( beta_of_Tb*(1-X_e)- n_H*alpha2_of_Tb*pow(X_e,2) );

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over 
  // Since the IC is at x=0 (tau(0) = 0) the array should go from 0.0 -> x_start 
  const int npts = 1000;
  Vector x_array = Utils::linspace(0.0, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

