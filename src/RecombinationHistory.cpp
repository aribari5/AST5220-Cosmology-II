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

  // Compute and spline sound horizon
  solve_for_sound_horizon();
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
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_saha_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);
  Vector ne_saha_arr(npts_rec_arrays);

  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB0     = cosmo->get_OmegaB(0.0);
  const double H0 = cosmo->get_H0() * Constants.km / Constants.Mpc;  // 1/s
  const double rho_crit0   = 3.0*pow(H0,2)/(8.0*Constants.pi*Constants.G);       // Critical density today in kg/m^3
  const double m_H = Constants.m_H;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    
    // Store the result in Xe_saha_arr...

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    Xe_saha_arr[i] = Xe_current;
    ne_arr[i] = ne_current;
    ne_saha_arr[i] = ne_current;

    // Debugging test
    // std::cout << "---------------------------------\n";
    // std::cout << "x: " << x_array[i] << " Xe_saha: " << Xe_current << " ne_saha: " << ne_current << "\n";
    // isnan(Xe_current) ? std::cout << "Its NaN" << "\n"
    //          : std::cout << "Its a real number" << "\n";
    // isnan(ne_current) ? std::cout << "Its NaN" << "\n"
    //          : std::cout << "Its a real number" << "\n";

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!) 
      // Implement rhs_peebles_ode
      //==============================================================
  
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result. Store in Xe_arr 
      //=============================================================================

      Vector Xe_ic{Xe_current};       // vector with i.c. for Xe

      // x-array from THIS point to end
      Vector x_peebles(x_array.begin()+i, x_array.end()); // x_peebles[i] corresponds to the current value of x, when the Peebles approx. is activated

      // Solve ODE

      ODESolver peebles_Xe_ode;

      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
  
      peebles_Xe_ode.solve(dXedx, x_peebles, Xe_ic, gsl_odeiv2_step_rkf45);

      Vector Xe_peebles = peebles_Xe_ode.get_data_by_component(0);

      // Fill remaining arrays

      for(size_t j = 0; j < Xe_peebles.size(); j++){
        Xe_arr[i + j] = Xe_peebles[j];

        double a = exp(x_peebles[j]);
        double n_b = OmegaB0 * rho_crit0 / (m_H * pow(a,3));
        double n_H = (1 - Yp) * n_b;

        ne_arr[i + j] = Xe_peebles[j] * n_H;
      }

      break; 
      
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x,
  // functions are working
  //=============================================================================
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe of x");
  Xe_saha_of_x_spline.create(x_array, Xe_saha_arr, "Xe Saha of x");

  ne_of_x_spline.create(x_array, ne_arr, "ne of x");
  ne_saha_of_x_spline.create(x_array, ne_saha_arr, "ne Saha of x");

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
  
  const double H0          = cosmo->get_H0() * Constants.km / Constants.Mpc;     // 1/s
  const double rho_crit0   = 3.0*pow(H0,2)/(8.0*Constants.pi*Constants.G);       // Critical density today in kg/m^3
  const double TCMB0       = cosmo->get_TCMB(0.0);                               // CMB temperature today in K


  const double n_b         = OmegaB0*rho_crit0/(m_H * pow(a,3));             // Number density of baryons at x in 1/m^3
  const double n_H         = (1-Yp)*n_b;                                     // Number density of hydrogen at x in 1/m^3
  const double T_b         = TCMB0/a;                                        // Temperature of baryons at x in K


  const double C           = (1/n_b)*pow( (m_e*k_b*T_b)/(2*Constants.pi*hbar*hbar), 1.5 )*exp(-epsilon_0/(k_b*T_b)); 

  //=============================================================================
  // Computing Xe and ne from the Saha equation
  //=============================================================================
  
  // Electron fraction and number density
  

  double Xe;
  if (C > 400.0){
    Xe = 1.0;
  }
  else{
    Xe = (-C + sqrt(C*C + 4.0*C)) / 2.0;
  }

  // debugging
  // list of x test values
//   std::vector<double> test_x_values = {-12.0, -10.0, -8.0, -6.0, -4.0, -2.0, 0.0}; 
//   const double epsilon = 0.05; // tolerance for comparing x values

//  for (const double test_x : test_x_values) {
//   if ((test_x - epsilon < x) && (x < test_x + epsilon)) {
//     std::cout << "---------------------------------\n";
//     std::cout << "Xe is = " << Xe << " at x = " << x << ":\n";
//     std::cout << "n_H is = " << n_H << " at x = " << x << ":\n";
//     std::cout << "---------------------------------\n";
//   }
// }


  double ne = Xe*n_H;
  
  // debugging for loop checking if any values in Xe or ne are NaN
  // std::cout << "---------------------------------\n";
  // std::cout << "---------Debugging NaNs----------\n";
  // std::cout << "---------------------------------\n";
  // std::cout << "x: " << x << " C: " << C << "\n";
  // std::cout << "T_b: " << T_b << "\n";
  // std::cout << "epsilon/(kT): " << epsilon_0/(k_b*T_b) << "\n";
  // std::cout << "n_b: " << n_b << " n_H: " << n_H << "\n";
  // std::cout << "OmegaB: " << OmegaB << " OmegaB0: " << OmegaB0 << "\n";
  // std::cout << "rho_crit0: " << rho_crit0 << "\n";
  // std::cout << "m_H: " << m_H << "\n";

  // std::cout << "---------------------------------\n";
  // for (int i = 0; i < 1; i++) {
  //   if (std::isnan(T_b)) {
  //     std::cout << "T_b is NaN at x = " << x << "\n";
  //     break;
  //   }
  //   else {
  //     std::cout << "T_b is a real number " << T_b << " at x = " << x << "\n";
  //   }
  // }
  // std::cout << "---------------------------------\n";
  // for (int i = 0; i < 1; i++) {
  //   if (std::isnan(C)) {
  //     std::cout << "C is NaN at x = " << x << "\n";
  //     break;
  //   }
  //   else {
  //     std::cout << "C is a real number " << C << " at x = " << x << "\n";
  //   }
  // }
  // std::cout << "---------------------------------\n";

  // for (int i = 0; i < 1; i++) {
  //   if (std::isnan(Xe)) {
  //     std::cout << "Xe is NaN at x = " << x << "\n";
  //     break;
  //   }
  //   else {
  //     std::cout << "Xe is a real number " << Xe << " at x = " << x << "\n";
  //   }
  // }
  // std::cout << "---------------------------------\n";

  // for (int i = 0; i < 1; i++) {
  //   if (std::isnan(ne)) {
  //     std::cout << "ne is NaN at x = " << x << "\n";
  //     break;
  //   }
  //   else {
  //     std::cout << "ne is a real number " << ne << " at x = " << x << "\n";
  //   }

  // }


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

  const double H0_over_h   = Constants.H0_over_h;

  // Cosmological parameters

  const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB0     = cosmo->get_OmegaB(0.0);
  const double OmegaCDM    = cosmo->get_OmegaCDM();
  const double OmegaM      = cosmo->get_OmegaM();
  const double OmegaNu     = cosmo->get_OmegaNu();
  const double OmegaR      = cosmo->get_OmegaR();
  const double OmegaLambda = cosmo->get_OmegaLambda();
  
  const double H           = cosmo->H_of_x(x);
  const double H0          = cosmo->get_H0() * Constants.km / Constants.Mpc;     // 1/s
  const double rho_crit0   = 3.0*pow(H0,2)/(8.0*Constants.pi*Constants.G);       // Critical density today in kg/m^3
  const double TCMB0       = cosmo->get_TCMB(0.0);                                                   // CMB temperature today in K
  const double T_b         = TCMB0/a;


  // All of the constants needed for RHS of Peebles

  const double n_b          = OmegaB0*rho_crit0/(m_H * pow(a,3));
  const double n_H          = (1-Yp)*n_b; 
  const double n_1s         = (1-X_e)*n_H;
  const double lambda_alpha = H*pow(3*epsilon_0,3)/( pow(8*Constants.pi,2)*n_1s );
  const double alpha        = 1/137.0359992;
  const double phi2_of_Tb   = 0.448*log(epsilon_0/((k_b*T_b)));
  const double alpha2_of_Tb = (64*Constants.pi)/(sqrt(27*Constants.pi)) * pow(alpha/m_e,2) * sqrt(epsilon_0/(k_b*T_b)) * phi2_of_Tb;
  const double beta_of_Tb   = alpha2_of_Tb*pow(m_e*k_b*T_b/(2*Constants.pi*hbar*hbar),1.5)*exp(-epsilon_0/(k_b*T_b));
  const double beta2_of_Tb  = beta_of_Tb*exp((3*epsilon_0)/(4*(k_b*T_b)));
  const double Cr_of_Tb     = (lambda_2s1s + lambda_alpha)/(lambda_2s1s + lambda_alpha + beta2_of_Tb);

  //=============================================================================
  // RHS of Peebles ODE for dXedx
  //=============================================================================

  // Course website advice to check if we are fully recombined and if so set derivative to zero to avoid overflow and NaN's.
  const double ratio = epsilon_0 / (k_b * T_b);
    if (ratio > 70.0) {
        dXedx[0] = 0.0;          // fully recombined, implying derivative exactly zero
        return GSL_SUCCESS;
    }

  
  dXedx[0] =  Cr_of_Tb/H* ( beta_of_Tb*(1-X_e)- n_H*alpha2_of_Tb*pow(X_e,2) );

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  
  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    
  
    const double ne = ne_of_x(x);
    const double H  = cosmo->H_of_x(x);

    dtaudx[0] = -(Constants.c*Constants.sigma_T*ne)/H;

    return GSL_SUCCESS;
  };
  
  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  
  ODESolver tau_ode_solver;

  // Set up x-arrays to integrate over 
  // Since the IC is at x=0 (tau(0) = 0) the array should go from x_end -> x_start 

  const int npts = 10000;
  Vector x_array = Utils::linspace(x_end, x_start, npts);

  double tau_initial = 0.0;     
  Vector tau_ic{tau_initial};                   // vector with i.c. for tau

  tau_ode_solver.solve(dtaudx, x_array, tau_ic,gsl_odeiv2_step_rkf45);//    ,gsl_odeiv2_step_rkf45

  auto tau_array = tau_ode_solver.get_data_by_component(0);         // get the 0th component of the sol.




  tau_of_x_spline.create(x_array, tau_array, "tau of x");         // create spline
  
 // Similar for SAHA regime only
 ODESolver tau_saha_ode_solver;
  Vector x_array_saha = Utils::linspace(x_end, x_start, npts);  // same grid

  ODEFunction dtaudx_saha = [&](double x, const double *tau, double *dtaudx){
      const double ne_saha = ne_saha_of_x(x);        // using Saha ne
      const double H  = cosmo->H_of_x(x);
      dtaudx[0] = -(Constants.c*Constants.sigma_T*ne_saha) / H;
      return GSL_SUCCESS;
  };

  tau_saha_ode_solver.solve(dtaudx_saha, x_array_saha, tau_ic, gsl_odeiv2_step_rkf45); // same IC
  auto tau_saha_array = tau_saha_ode_solver.get_data_by_component(0);

  tau_of_x_spline_saha.create(x_array_saha, tau_saha_array, "tau saha");

  //=============================================================================
  // Computing the visibility functions and splining everything
  //=============================================================================
  
  // Calculating "available" values in the visibility function, then splining
  Vector g_tilde_array(npts);

    for(int i = 0; i < npts; i++){
      double x = x_array[i];
      g_tilde_array[i] = -dtaudx_of_x(x) * exp(-tau_of_x(x));
    }

  g_tilde_of_x_spline.create(x_array, g_tilde_array, "g tilde");


  Utils::EndTiming("opticaldepth");
}

void RecombinationHistory::solve_for_sound_horizon(){
  Utils::StartTiming("soundhorizon");
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){
      

      const double Omega_gamma0 = cosmo->get_OmegaR(0.0) - cosmo->get_OmegaNu(0.0);
      const double Omega_b0     = cosmo->get_OmegaB(0.0);

      const double R            = 4.0*Omega_gamma0*exp(-x) / (3.0*Omega_b0);
      const double c_s          = Constants.c * sqrt(R/(3.0*(1.0+R)));

      const double Hp           = cosmo->Hp_of_x(x);

    dsdx[0]= c_s / Hp;

    return GSL_SUCCESS;
  };

    const double Omega_gamma0 = cosmo->get_OmegaR(0.0) - cosmo->get_OmegaNu(0.0);
    const double Omega_b0     = cosmo->get_OmegaB(0.0);

    const double R_initial    = 4.0*Omega_gamma0*exp(-x_start) / (3.0*Omega_b0);

    const double c_s_initial  = Constants.c * sqrt(R_initial/(3.0*(1.0+R_initial))); 

    const double Hp_initial   = cosmo->Hp_of_x(x_start);


  ODESolver sound_horizon_ode_solver;

  const int npts = 10000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  double s_initial = c_s_initial/Hp_initial;
  Vector s_ic{s_initial};                   // vector with i.c. for s

  sound_horizon_ode_solver.solve(dsdx, x_array, s_ic,gsl_odeiv2_step_rkf45);

  auto s_array = sound_horizon_ode_solver.get_data_by_component(0);

  // creating spline

  sound_horizon_of_x_spline.create(x_array, s_array, "sound horizon of x");


  Utils::EndTiming("soundhorizon");
}


//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::tau_of_x_saha(double x) const{
  return tau_of_x_spline_saha(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return ne_of_x_spline(x);
}

double RecombinationHistory::ne_saha_of_x(double x) const{
    return ne_saha_of_x_spline(x);
}

double RecombinationHistory::sound_horizon_of_x(double x) const{
  return sound_horizon_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout <<"====================================================\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout <<"====================================================\n";

  std::cout << "---------------------------------\n"; 
  std::cout << "Yp:          " << Yp          << "\n";

  // std::cout << "---------------------------------\n";
  // std::cout << "Optical depth today:\n";
  // std::cout << "tau(x=0) = "
  //           << tau_of_x(0.0)
  //           << "\n";
  // std::cout << "---------------------------------\n";
  // std::cout << "Optical depth at x=-5:\n";
  // std::cout << "tau(x=-5) = "
  //           << tau_of_x(-5.0)
  //           << "\n";
  // std::cout << "---------------------------------\n";
  // std::cout << "Optical depth at x=-7:\n";
  // std::cout << "tau(x=-7) = "
  //           << tau_of_x(-7.0)
  //           << "\n";
  // std::cout << "---------------------------------\n"; 
  // std::cout << "Optical depth at x=-12:\n";
  // std::cout << "tau(x=-12) = "
  //           << tau_of_x(-12.0)
  //           << "\n";
  // std::cout << "---------------------------------\n"; 


  // Computing tau at x_decoupling, where tau = 1.0
  // Evaluating X_e at decoupling, and finding x where X_e = 0.5:

  const double tau_at_decoupling = 1.0;

  double x_decoupling             = Utils::binary_search_for_value(tau_of_x_spline, tau_at_decoupling);

  double Xe_decoupling            = Xe_of_x(x_decoupling);             
  double x_Xe_half                = Utils::binary_search_for_value(Xe_of_x_spline, 0.5, {-8.0, -6.0});
  double t_decoupling             = cosmo->t_of_x(x_decoupling) / Constants.Gyr; 

  // Computing tau at x_decoupling, where tau = 1.0
  // Evaluating X_e at decoupling, and finding x where X_e = 0.5:
  // SAHA approx
  double x_decoupling_saha        = Utils::binary_search_for_value(tau_of_x_spline_saha, tau_at_decoupling);

  double Xe_decoupling_saha       = Xe_of_x(x_decoupling_saha);             
  double x_Xe_half_saha           = Utils::binary_search_for_value(Xe_saha_of_x_spline, 0.5, {-8.0, -6.0});
  double t_decoupling_saha        = cosmo->t_of_x(x_decoupling_saha) / Constants.Gyr; 


  std::cout << "x_decoupling                         = " << x_decoupling << "\n";
  std::cout << "Redshift of decoupling: z_decoupling = " << exp(-x_decoupling) - 1.0 << "\n";
  std::cout << "Time of decoupling: t_decoupling     = " << t_decoupling << " Gyr\n";

  std::cout << "---------------------------------\n";
  std::cout << "x_decoupling (Saha)                         = " << x_decoupling_saha << "\n";
  std::cout << "Redshift of decoupling (Saha): z_decoupling = " << exp(-x_decoupling_saha) - 1.0 << "\n";
  std::cout << "Time of decoupling (Saha): t_decoupling     = " << t_decoupling_saha << " Gyr\n";

  std::cout << "---------------------------------\n";


  std::cout << "X_e at decoupling: Xe(x_decoupling)             = " << Xe_decoupling << "\n";
  std::cout << "X_e at decoupling (Saha): Xe(x_decoupling_saha) = " << Xe_decoupling_saha << "\n";

  std::cout << "---------------------------------\n";

  std::cout << "x where X_e = 0.5: x_Xe_half                          = " << x_Xe_half << "\n";
  std::cout << "Redshift where X_e = 0.5: z_Xe_half                   = " << exp(-x_Xe_half) - 1.0 << "\n";

  std::cout << "x where X_e = 0.5 (Saha): x_Xe_half_saha              = " << x_Xe_half_saha << "\n";
  std::cout << "Redshift where X_e = 0.5 (Saha): z_Xe_half_saha       = " << exp(-x_Xe_half_saha) - 1.0 << "\n";


  std::cout << "---------------------------------\n";
  std::cout << "Sound horizon at last scattering (x = -7): "
          << sound_horizon_of_x(-7.0) / Constants.Mpc << " Mpc\n";
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
    fp << x                    << " ";  // 0, dimensionless
    fp << Xe_of_x(x)           << " ";  // 1, dimensionless
    fp << ne_of_x(x)           << " ";  // 2, 1/m^3
    fp << tau_of_x(x)          << " ";  // 3, dimensionless
    fp << dtaudx_of_x(x)       << " ";  // 4, dimensionless
    fp << ddtauddx_of_x(x)     << " ";  // 5, dimensionless
    fp << g_tilde_of_x(x)      << " ";  // 6, normalized
    fp << dgdx_tilde_of_x(x)   << " ";  // 7, normalized
    fp << ddgddx_tilde_of_x(x) << " ";  // 8, normalized
    fp << sound_horizon_of_x(x) << " ";  // 9, m
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

