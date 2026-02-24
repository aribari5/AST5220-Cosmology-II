#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  //=============================================================================
  // TODO: Compute OmegaR, OmegaNu, OmegaLambda, H0, ...
  //=============================================================================

  H0 = 100.0*h;                                                                 // Hubble parameter today in km/s/Mpc
  H0_SI = H0 * Constants.km / Constants.Mpc;                                    // Hubble parameter today in 1/s
  OmegaR = ( pow(Constants.pi,3)*pow(Constants.k_b*TCMB,4)*8*Constants.G ) 
  / ( 15.0*pow(Constants.c,5)*pow(Constants.hbar,3)*3.0*pow(H0_SI,2) );         // Radiation density today

  OmegaNu = Neff*(7.0/8.0)*pow(4.0/11.0,4.0/3.0)*OmegaR;                        // Neutrino density today

  OmegaLambda = 1.0 - OmegaB - OmegaCDM - OmegaR - OmegaNu - OmegaK;            // Dark energy density today
}

//====================================================
// Do all the solving. Computing eta(x) and t(x), splining.
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  x_start = -10.0;
  x_end   = 0.0;
  int npts    = 100;

  Vector x_array = Utils::linspace(x_start, x_end, npts);



  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // The rhs of the detadx ODE
    //=============================================================================

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Setting initial condition, solving the ODE and making the spline.
  //=============================================================================
  
  double eta_initial = 0.0;     
  Vector eta_ic{eta_initial};       // vector with i.c. for eta

  ODESolver ode_solver;

  ode_solver.solve(detadx, x_array, eta_ic,gsl_odeiv2_step_rkf45);    // solve. Had to include a stepper since it expected 4 argumetns.
  
  Vector eta_array = ode_solver.get_data_by_component(0);         // get the 0th component of the sol.

  eta_of_x_spline.create(x_array, eta_array, "Eta of x");         // create spline

  std::cout << "---------------------------------\n";
  std::cout << "Conformal time today:\n";
  std::cout << "eta(x=0) = "
            << eta_of_x(0.0)/(Constants.c)/Constants.Gyr
            << " Gyr\n";
  std::cout << "---------------------------------\n";


  Utils::EndTiming("Eta");


  // Similar for t(x)
  Utils::StartTiming("t of x");

  //x_start = -10.0;
  //x_end   = 0.0;
  //int npts    = 100;
  // and
  //Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The (rhs of) ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    double H_SI = H_of_x(x) * Constants.km / Constants.Mpc;             // Hubble parameter today in 1/s 
    dtdx[0] = 1/H_SI;

    return GSL_SUCCESS;
  };


  // Setting initial condition, solving the ODE and making the spline.

  double t_initial = 1/(2*H_of_x(x_start));     // in the r. dom. era, t = 1/2H(x)

  Vector t_ic{t_initial};                   // vector with i.c. for t

  //still using ODESolver ode_solver;

  ode_solver.solve(dtdx, x_array, t_ic,gsl_odeiv2_step_rkf45);    // solve. Had to include a stepper since it expected 4 argumetns.
  
  Vector t_array = ode_solver.get_data_by_component(0);         // get the 0th component of the sol.

  t_of_x_spline.create(x_array, t_array, "t of x");         // create spline

  std::cout << "---------------------------------\n";
  std::cout << "Age of universe today:\n";
  std::cout << "t(x=0) = "
            << t_of_x(0.0)/(Constants.Gyr)
            << " Gyr\n";
  std::cout << "---------------------------------\n";

  Utils::EndTiming("t of x");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  //=============================================================================
  // The Hubble parameter as a function of x = exp(a). 
  //=============================================================================

  double H = H0 * sqrt( 
    (OmegaB + OmegaCDM)*exp(-3.0*x) 
    + (OmegaR + OmegaNu)*exp(-4.0*x) 
    + OmegaLambda 
    + OmegaK*exp(-2.0*x) );

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  //=============================================================================
  // The conformal Hubble parameter as a function of x = exp(a). Using Hp = a*H = exp(x)*H
  //=============================================================================

  double H_SI = H_of_x(x) * Constants.km / Constants.Mpc;  // in 1/s
  double Hp_SI = exp(x) * H_SI;  // conformal Hubble in 1/s
  return Hp_SI;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  //=============================================================================
  // The derivative of the conformal Hubble parameter with respect to x. 
  //=============================================================================

  // Intermediate variables
  double H = H_of_x(x);         
  double Hp = Hp_of_x(x);
  


  double dHpdx = Hp + exp(x)*pow(H0,2.0)/(2*H)
                * ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) );

  return dHpdx;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
  //=============================================================================
  // The double derivative of the conformal Hubble parameter with respect to x. 
  //=============================================================================


  // Intermediate variables
  double H = H_of_x(x);         
  double Hp = Hp_of_x(x);
  double dHpdx = dHpdx_of_x(x);
  double dHdx = pow(H0,2.0)/(2.0*H) *
                ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) );




  double ddHpddx = dHpdx + pow(H0,2.0)/2.0 * 
                  (exp(x)*H-exp(x)*dHdx)/pow(H,2.0) *
                  ( -3.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    -4.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    -2.0*OmegaK*exp(-2.0*x) ) +
                    exp(x)/H * 
                    ( 9.0*(OmegaB + OmegaCDM)*exp(-3.0*x) 
                    +16.0*(OmegaR + OmegaNu)*exp(-4.0*x) 
                    +4.0*OmegaK*exp(-2.0*x));
  return ddHpddx;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  //=============================================================================
  // The baryon density as a function of x = exp(a).
  //=============================================================================

  double OmegaB_x = OmegaB * exp(3.0*x) * pow(H0/H_of_x(x),2.0);

  return OmegaB_x;
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  //=============================================================================
  // The radiation density as a function of x = exp(a).
  //=============================================================================
  
  double OmegaR_x = OmegaR * exp(4.0*x) * pow(H0/H_of_x(x),2.0);

  return OmegaR_x;
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  //=============================================================================
  // The neutrino density as a function of x = exp(a).
  //=============================================================================
  
  double OmegaNu_x = OmegaNu * exp(4.0*x) * pow(H0/H_of_x(x),2.0);

  return OmegaNu_x;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  //=============================================================================
  // The CDM density as a function of x = exp(a).
  //=============================================================================
    
  double OmegaCDM_x = OmegaCDM * exp(3.0*x) * pow(H0/H_of_x(x),2.0);


  return OmegaCDM_x;
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  //=============================================================================
  // The dark energy density as a function of x = exp(a).
  //=============================================================================
  
  double OmegaLambda_x = OmegaLambda * pow(H0/H_of_x(x),2.0);

  return OmegaLambda_x;
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  //=============================================================================
  // The curvature density as a function of x = exp(a).
  //=============================================================================
  
  double OmegaK_x = OmegaK * exp(2.0*x) * pow(H0/H_of_x(x),2.0);

  return OmegaK_x;
}
double BackgroundCosmology::get_OmegaM(double x) const{ 
  //=============================================================================
  // The total matter density as a function of x = exp(a). (i saw this in the header but they were not in the cpp file so i added them here)
  //=============================================================================
  
  return get_OmegaB(x) + get_OmegaCDM(x);
}
 double BackgroundCosmology::get_OmegaRtot(double x) const{ 
  //=============================================================================
  // The total radiation density as a function of x = exp(a). (i saw this in the header but they were not in the cpp file so i added them here)
  //=============================================================================
  
  return get_OmegaR(x) + get_OmegaNu(x);
}

 double BackgroundCosmology::get_OmegaMnu(double x) const{ 
  //=============================================================================
  // The total matter and neutrino density as a function of x = exp(a). (i saw this in the header but they were not in the cpp file so i added them here) (also idk what it is)
  //=============================================================================
  
  return get_OmegaM(x) + get_OmegaNu(x);
 }

double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // The luminosity distance as a function of x = exp(a).
  //=============================================================================
  
  
  double d_L_seconds = get_angular_distance_of_x(x) * exp(-2.0 * x);
  double d_L = d_L_seconds  / Constants.Mpc / 1e3;                    // Gpc

  return d_L;
}

double BackgroundCosmology::get_r_of_x(double x) const{
  //=============================================================================
  // r for different curvature cases as a function of x = exp(a).
  //=============================================================================
  
  double chi = get_comoving_distance_of_x(x);
  
  

  if (std::abs(OmegaK) < 1e-10){                                      // Flat universe

    return chi;
  }
  
  
  double omega_arg = sqrt( abs(OmegaK) )  * (H0*chi)/Constants.c;

  if (OmegaK < 0.0){                                                 // Closed universe

    
    return chi * ( sin(omega_arg)/(omega_arg) );
    }

  else {                                                              // Open universe

    return chi * ( sinh(omega_arg)/(omega_arg) );
  }

  
}

double BackgroundCosmology::get_angular_distance_of_x(double x) const{
  //=============================================================================
  // The angular distance as a function of x = exp(a).
  //=============================================================================
  
  double d_A = exp(x)*get_r_of_x(x);

  return d_A;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // The comoving distance as a function of x = exp(a).
  //=============================================================================
  
  double comoving_distance = eta_of_x(0.0) - eta_of_x(x);

  return comoving_distance;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                               << " ";
    fp << eta_of_x(x)                     << " ";
    fp << t_of_x(x)                       << " ";
    fp << Hp_of_x(x)                      << " ";
    fp << dHpdx_of_x(x)                   << " ";
    fp << get_OmegaB(x)                   << " ";
    fp << get_OmegaCDM(x)                 << " ";
    fp << get_OmegaLambda(x)              << " ";
    fp << get_OmegaR(x)                   << " ";
    fp << get_OmegaNu(x)                  << " ";
    fp << get_OmegaK(x)                   << " ";
    fp << get_luminosity_distance_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

