#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  // the k and x arrays

  k_array = exp(Utils::linspace(log(k_min),log(k_max), n_k));
  x_array = Utils::linspace(-18.0, x_end, n_x);   // putting a x_start = -18 cap, since tau'' drops at the boundary.

  Psi_array       = Vector(n_x * n_k);
  Pi_array        = Vector(n_x * n_k);

  Theta_array     = Vector2D(Constants.n_ell_theta, Vector(n_x*n_k));

  int n_ell_tot_full = Constants.n_ell_tot_full;

  // This is the array which will hold our solutions
  this->y_array = Vector2D(n_ell_tot_full);
  for (int i = 0; i < n_ell_tot_full; i++)
    this->y_array[i] = Vector(n_x * n_k);

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find x_tc end and index of x_tc in the x_array
    std::pair<double,int> tc_time = get_tight_coupling_time(k, x_array);

    double x_end_tight = tc_time.first;
    int idx_end = tc_time.second;

    //===================================================================
    // Tight coupling integration
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight


    // debugging

    if(idx_end <= 1){
    std::cout << "idx_end too small: " << idx_end << std::endl;
    throw "Tight coupling interval too short";  
      }
    else(std::cout << "Tight coupling ends at x = " << x_end_tight << " with index " << idx_end << std::endl);


    Vector x_tc(x_array.begin(), x_array.begin() + idx_end+1);

    ODESolver solver_tc;
    solver_tc.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
    

    // Sol at the end of tc
    Vector y_tight_coupling(Constants.n_ell_tot_tc);
    for(int i = 0; i < Constants.n_ell_tot_tc; ++i){
        y_tight_coupling[i] = solver_tc.get_data_by_component(i).back();
    }

    //====i===============================================================
    //Full equation integration
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    Vector x_full(x_array.begin() + idx_end, x_array.end());
    ODESolver solver_full;
    solver_full.solve(dydx_full, x_full, y_full_ini, gsl_odeiv2_step_rkf45);

    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
  for(int ix = 0; ix < n_x; ++ix){

      const int index = ix + n_x * ik;

      

      // are we in tc or full reigime
      if(ix < idx_end){ 
        // if yes, the solutions are in solver_tc

        // we fetch the ix-th point, which is in tc.

        double a          = exp(x_array[ix]);
        double H0         = cosmo->get_H0();
        double Hp         = cosmo->Hp_of_x(x_array[ix]);
        double OmegaR     = cosmo->get_OmegaR();
        double tau_prime  = rec->dtaudx_of_x(x_array[ix]);

        double ck_over_Hp = (Constants.c * k) / Hp;

        y_array[Constants.ind_deltacdm][index] = solver_tc.get_data_by_component(Constants.ind_deltacdm_tc)[ix];
        y_array[Constants.ind_deltab][index]   = solver_tc.get_data_by_component(Constants.ind_deltab_tc)[ix];
        y_array[Constants.ind_vcdm][index]     = solver_tc.get_data_by_component(Constants.ind_vcdm_tc)[ix];
        y_array[Constants.ind_vb][index]       = solver_tc.get_data_by_component(Constants.ind_vb_tc)[ix];
        y_array[Constants.ind_Phi][index]      = solver_tc.get_data_by_component(Constants.ind_Phi_tc)[ix];

        for(int ell = 0; ell < Constants.n_ell_theta_tc; ++ell){
          y_array[Constants.ind_start_theta+ell][index] = solver_tc.get_data_by_component(Constants.ind_start_theta+ell)[ix];
        }

        // Theta2 is computed using an approx. in tc. #override
        double Theta1 = y_array[Constants.ind_start_theta + 1][index];
        double Theta2 = - (20.0 / 45.0) * ck_over_Hp / tau_prime * Theta1;

        y_array[Constants.ind_start_theta + 2][index] = Theta2;
        

        
      }
      else { 
        // if no, the solutions are in solver_full.
        // we then need a new idx for the full solution since it starts at x_end_tight, not x_start
        int ix_full = ix - idx_end;

        y_array[Constants.ind_deltacdm][index] = solver_full.get_data_by_component(Constants.ind_deltacdm)[ix_full];
        y_array[Constants.ind_deltab][index]   = solver_full.get_data_by_component(Constants.ind_deltab)[ix_full];
        y_array[Constants.ind_vcdm][index]     = solver_full.get_data_by_component(Constants.ind_vcdm)[ix_full];
        y_array[Constants.ind_vb][index]       = solver_full.get_data_by_component(Constants.ind_vb)[ix_full];
        y_array[Constants.ind_Phi][index]      = solver_full.get_data_by_component(Constants.ind_Phi)[ix_full];

        for(int ell = 0; ell < Constants.n_ell_theta; ++ell){
            y_array[Constants.ind_start_theta+ell][index] = solver_full.get_data_by_component(Constants.ind_start_theta+ell)[ix_full];
            }
      }


      double a          = exp(x_array[ix]);
      double H0         = cosmo->get_H0();
      double OmegaR     = cosmo->get_OmegaR();

      double Theta2 = y_array[Constants.ind_start_theta + 2][index];
    
      Psi_array[index] = -y_array[Constants.ind_Phi][index]- (12.0*H0*H0*OmegaR*Theta2)/pow(a*k*Constants.c, 2);
      Pi_array[index]  = Theta2;   // for now
    }
  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_array, k_array, y_array[Constants.ind_deltacdm]);
  delta_b_spline.create(x_array, k_array, y_array[Constants.ind_deltab]);
  v_cdm_spline.create(x_array, k_array, y_array[Constants.ind_vcdm]);
  v_b_spline.create(x_array, k_array, y_array[Constants.ind_vb]);
  Phi_spline.create(x_array, k_array, y_array[Constants.ind_Phi]);

  Psi_spline.create(x_array, k_array, Psi_array);
  Pi_spline.create(x_array, k_array, Pi_array);

  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for(int ell = 0; ell < Constants.n_ell_theta; ++ell){
        Theta_spline[ell].create(x_array, k_array, y_array[Constants.ind_start_theta+ell]);

  }
}


//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  // note that for the tc regime, we only need l = 0,1 for photons
  //=============================================================================
  
  // IC for scalar quantities (Gravitational potential, baryons and CDM)

  double f_nu          = 0.0;               //cosmo->get_OmegaNu() / cosmo->get_OmegaR();   // 0 if no neturinoes included
  double Hp            = cosmo->Hp_of_x(x);


  double Psi_ic        = -1.0/(3.0/2.0 + 2.0/5.0*f_nu);
  double Phi_ic        = -(1.0 + 2.0/5.0*f_nu )* Psi_ic;
  double delta_cdm_ic  = -3.0/2.0 * Psi_ic;
  double delta_b_ic    = delta_cdm_ic;
  double v_cdm_ic      = -(Constants.c*k/(2*Hp)) * Psi_ic;
  double v_b_ic        = v_cdm_ic;

  
  // IC for photon temperature perturbations (Theta_ell)

  double Theta0_ic = -0.5 * Psi_ic;
  double Theta1_ic = (Constants.c*k/(6.0*Hp)) * Psi_ic;

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
    // Work in progress. Since I have chosen to not include Neutrinos (for now),
    // I will instead use my time wisely and focus on the other parts of the code
    // to get it running.
  }

  // Store the IC in the y_tc vector

  // Phi       = Phi_ic;
  // delta_cdm = delta_cdm_ic;
  // delta_b   = delta_b_ic;
  // v_cdm     = v_cdm_ic;
  // v_b       = v_b_ic;

  

  y_tc[Constants.ind_Phi_tc]      = Phi_ic;
  y_tc[Constants.ind_deltacdm_tc] = delta_cdm_ic;
  y_tc[Constants.ind_deltab_tc]   = delta_b_ic;
  y_tc[Constants.ind_vcdm_tc]     = v_cdm_ic;
  y_tc[Constants.ind_vb_tc]       = v_b_ic;
  
  Theta[0]  = Theta0_ic;
  Theta[1]  = Theta1_ic;
  

  // debugging

  std::cout << "ck/Hp         = " << Constants.c*k/Hp << std::endl;

  std::cout << "IC. set at x  = " << x << " for k = " << k << std::endl;
  std::cout << "Phi_ic        = " << Phi_ic << std::endl;
  std::cout << "delta_cdm_ic  = " << delta_cdm_ic << std::endl;
  std::cout << "delta_b_ic    = " << delta_b_ic << std::endl;
  std::cout << "v_cdm_ic      = " << v_cdm_ic << std::endl;
  std::cout << "v_b_ic        = " << v_b_ic << std::endl;

  std::cout << "Theta0_ic     = " << Theta0_ic << std::endl;
  std::cout << "Theta1_ic     = " << Theta1_ic << std::endl;


  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];
  double *Theta_p         = &y[Constants.ind_start_thetap];
  double *Nu              = &y[Constants.ind_start_nu];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // IC for salar quantities (Gravitational potental, baryons and CDM), taken from the tight coupling regime

  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;
  Phi       = Phi_tc;



  // IC for photon temperature perturbations (Theta_ell), including higher order now

  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];


  double Hp        = cosmo->Hp_of_x(x);
  double tau_prime = rec->dtaudx_of_x(x);


  for(int ell = 2; ell < n_ell_theta; ell++){
  Theta[ell] = - (double(ell)/(2.0*ell + 1.0)) * (Constants.c*k)/(Hp*tau_prime) * Theta[ell-1];
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // Not including polarization for now.
    // Maybe later:)
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // ...
    // ...
    // Not including neutrinos for now.
    // Maybe later:)
  }

  // Store the IC in the y vector


  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double,int> Perturbations::get_tight_coupling_time(const double k,const Vector& x_array) const{
  double x_tight_coupling_end = 0.0;

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================


  double x_onset_of_rec = -8.3;

  // defining grid we search over

  const double x_start_search = x_start;
  const double x_end_search   = -6.0;

  const int    n_steps = 8000;                                        // such that dx is roughyl 0.001
  const double dx      = (x_end_search - x_start_search) / n_steps;  

  for(int i=0; i<= n_steps; i++){

    double x_i = x_start_search + i*dx;

      // condition 3) from course appendix
      if(x_i > x_onset_of_rec){
          x_tight_coupling_end = x_i; 
      }

      // relevant quantities for conditions
      double Hp = cosmo->Hp_of_x(x_i);
      double ck_over_Hp = (Constants.c * k) / Hp;
      double dtaudx = rec->dtaudx_of_x(x_i);
      double abs_dtaudx = std::abs(dtaudx);

      double min_val = 10.0*std::min(1.0,ck_over_Hp);

      // condition 2) from course appendix
      if (abs_dtaudx <= min_val) {
          x_tight_coupling_end = x_i;
          break;
      }

      // condition 1) from course appendix. is it even needed?
      // if (abs_dtaudx <= 10.0) {
      //     x_tight_coupling_end = x_i;
      //     break;
      // }

    // x_tight_coupling_end = x_i; // why did i put this here?

  }

  // final check
  if (x_tight_coupling_end > x_onset_of_rec) {
    std::cout << "x_tight_coupling_end = "<< x_tight_coupling_end << " is greater than x_onset_of_rec = " << x_onset_of_rec << ". Setting x_tight_coupling = x_onset_of_rec." << std::endl;
    x_tight_coupling_end = x_onset_of_rec;
  }

  // Find the index in the x_array corresponding to x_tight_coupling_end. useful to have.
  auto it = std::lower_bound(x_array.begin(), x_array.end(), x_tight_coupling_end);
  int idx_end = std::distance(x_array.begin(), it);


  return {x_tight_coupling_end, idx_end};
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================

  Vector k_array      = this-> k_array;
  Vector x_array      = this-> x_array;
  const auto& y_array = this-> y_array;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){

      
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================

      // const double Hp       = cosmo->Hp_of_x(x);

      const double tau       = rec->tau_of_x(x);
      const double g_tilde   = rec->g_tilde_of_x(x);
    
      double dPhidx = Phi_spline.deriv_x(x, k);
      double dPsidx = Psi_spline.deriv_x(x, k);


      double S1_term = g_tilde*(y_array[Constants.ind_start_theta][index]
                     + Psi_array[index]+(1.0/4.0) * Pi_array[index]);
      double S2_term = exp(-tau) * (dPsidx- dPhidx );
      double S3_term = 0.0; // 0 for now, no polarization
      double S4_term = 0.0; // ------------"------------- 

      // Temperature source
      ST_array[index] = S1_term + S2_term + S3_term + S4_term;

      // Polarization source
      if(Constants.polarization){
        SE_array[index] = 0.0; // 0 for now, no polarization
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");
  if(Constants.polarization){
    SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];



  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];


  double H0            = cosmo->get_H0();
  double Hp            = cosmo->Hp_of_x(x);
  double Hp_prime      = cosmo->dHpdx_of_x(x);
  double Omega_gamma0  = cosmo->get_OmegaR();
  double Omega_b0      = cosmo->get_OmegaB();
  double Omega_CDM0    = cosmo->get_OmegaCDM();

  double tau           = rec->tau_of_x(x);
  double tau_prime     = rec->dtaudx_of_x(x);
  double tau_2prime    = rec->ddtauddx_of_x(x);

  double R            = 4.0*Omega_gamma0*exp(-x) / (3.0*Omega_b0);
  double ck_over_Hp   = (Constants.c*k)/Hp;

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  
  double Theta2       = - (20.0/45.0) * ck_over_Hp/tau_prime * Theta[1];      // No polarization (for now)

  double Psi          = -Phi-12.0*H0*H0/(Constants.c*Constants.c*k*k)*(Omega_gamma0*Theta2)*exp(-2.0*x);


 
  

  // SET: Scalar quantities (Phi, delta, v, ...)
  
  dPhidx        = Psi - (1.0/3.0)*ck_over_Hp*ck_over_Hp*Phi + pow(H0/(2.0*Hp),2)*(Omega_CDM0*delta_cdm*exp(-x) 
                      + Omega_b0*delta_b*exp(-x) + 4.0*Omega_gamma0*Theta[0]*exp(-2.0*x));                           // No neutrinos

  //---------------------------------------

  double Theta0_prime = -ck_over_Hp*Theta[1] - dPhidx;
  
  double q_numerator   = -(((1.0-R)*tau_prime + (1.0+R)*tau_2prime)*(3.0*Theta[1]+v_b) - (ck_over_Hp)*Psi
                         + (1.0-Hp_prime/Hp)*(ck_over_Hp)*(-Theta[0]+2.0*Theta2) - (ck_over_Hp)*Theta0_prime);
  double q_denominator = (1.0+R)*tau_prime + (Hp_prime/Hp) -1.0;

  //---------------------------------------

  double q             = q_numerator / q_denominator;

  double v_b_prime     = (1.0/(1.0+R))*(-v_b - ck_over_Hp*Psi + R*(q + ck_over_Hp*(-Theta[0]+2.0*Theta2) -ck_over_Hp*Psi) );
  double Theta1_prime  = (1.0/3.0)*(q-v_b_prime);
  
  //---------------------------------------

  ddelta_cdmdx  = ck_over_Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx    = ck_over_Hp * v_b - 3.0*dPhidx;
  
  
  dv_cdmdx      = -v_cdm+ck_over_Hp*Psi;
  dv_bdx        = -v_b-ck_over_Hp*Psi+tau_prime*R*(3.0*Theta1_prime + v_b_prime);

  // debugging
  
  std::cout << "x            = " << x << std::endl;
  std::cout << "k            = " << k << std::endl;
  std::cout << "H0           = " << H0 << std::endl;
  std::cout << "Hp           = " << Hp << std::endl;
  std::cout << "Hp'          = " << Hp_prime << std::endl;
  std::cout << "Omega_gamma0 = " <<  Omega_gamma0<< std::endl;
  std::cout << "Omega_b0     = " << Omega_b0 << std::endl;
  std::cout << "Omega_CDM0   = " << Omega_CDM0 << std::endl;
  std::cout << "tau          = " << tau << std::endl;       // this one is just for debugging
  std::cout << "tau'         = " << tau_prime << std::endl;
  std::cout << "tau''        = " << tau_2prime << std::endl;

  std::cout << "Theta1 = " << Theta[1]  << std::endl;
  std::cout << "Theta2 = " << Theta2 << std::endl;
  std::cout << "Psi    = " << Psi << std::endl;
  

  std::cout << "dPhidx  = " << dPhidx << std::endl;
  std::cout << "Theta0' = " << Theta0_prime << std::endl;
  std::cout << "Theta1' = " << Theta1_prime << std::endl;
  std::cout << "q_numer = " << q_numerator << std::endl;
  std::cout << "q_denom = " <<  q_denominator << std::endl;
  std::cout << "q       = " << q << std::endl;
  std::cout << "v_b'    = " <<  v_b_prime<< std::endl;
  std::cout << "delta_cdm' = " << ddelta_cdmdx << std::endl;
  std::cout << "delta_b'   = " << ddelta_bdx << std::endl;
  std::cout << "v_cdm'     = " << dv_cdmdx << std::endl;
  std::cout << "v_b'       = " << dv_bdx << std::endl;


  // Photon multipoles (Theta_ell) (in tc regime only Th0 and Th1 evolved)
  dThetadx[0] = Theta0_prime;
  dThetadx[1] = Theta1_prime;

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // No neutrinos for now
    
  }
      // debugging dydx for tc regime, NaNs?
//   for(int i = 0; i < Constants.n_ell_tot_tc; i++){
//     if(!std::isfinite(dydx[i])){
//         std::cout << "NaN in dydx_"<<i<< "at x=" << x << " k=" << k << std::endl;
//         throw "NaN detected";
//     }
//   else{
//     std::cout << "dydx_"<<i<<"="<<dydx[i]<< "at x=" << x << " k=" << k << " is finite."<< std::endl;
//   }
// }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];
  
  double H0            = cosmo->get_H0();
  double Hp            = cosmo->Hp_of_x(x);
  double Hp_prime      = cosmo->dHpdx_of_x(x);
  double Omega_gamma0  = cosmo->get_OmegaR();
  double Omega_b0      = cosmo->get_OmegaB();
  double Omega_CDM0    = cosmo->get_OmegaCDM();
  double eta_of_x      = cosmo->eta_of_x(x);

  double R             = 4.0*Omega_gamma0*exp(-x) / (3.0*Omega_b0); 

  double tau_prime     = rec->dtaudx_of_x(x);
  double tau_2prime    = rec->ddtauddx_of_x(x);

  double ck_over_Hp    = (Constants.c*k)/Hp;

  double Pi            = Theta[2];                                             //+ Theta_p[0] + Theta_p[2], however these are 0 (for now)  
  double Theta2        = Theta[2];


  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // Scalar quantities (Phi, delta, v, ...)
  
  double Psi           = -Phi -12.0*H0*H0/(Constants.c*Constants.c*k*k)*(Omega_gamma0*Theta2)*exp(-2.0*x);             // No neutrinos
  dPhidx               = Psi - (1.0/3.0)*ck_over_Hp*ck_over_Hp*Phi + 0.5*pow(H0/(Hp),2)*(Omega_CDM0*delta_cdm*exp(-x) 
                       + Omega_b0*delta_b*exp(-x) + 4.0*Omega_gamma0*Theta[0]*exp(-2.0*x));                           // No neutrinos

  ddelta_cdmdx = ck_over_Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx   = ck_over_Hp * v_b - 3.0*dPhidx;

  dv_cdmdx     = -v_cdm + ck_over_Hp * Psi;
  dv_bdx       = -v_b - ck_over_Hp*Psi + tau_prime*R*(3.0*Theta[1] + v_b); 




  // Photon multipoles (Theta_ell)

  dThetadx[0] = -ck_over_Hp*Theta[1] - dPhidx;
  dThetadx[1] = (Constants.c*k)/(3.0*Hp)*(Theta[0]-2*Theta[2]+Psi) + tau_prime*(Theta[1] + (1.0/3.0)*v_b);
  
  for(int ell = 2; ell < n_ell_theta-1; ell++){

    dThetadx[ell] = ck_over_Hp*(1.0/(2.0*ell+1))*(ell*Theta[ell-1] - (ell+1.0)*Theta[ell+1]) + tau_prime*(Theta[ell]-0.1*Pi*(ell==2));
  }

    // for the last multipole 
    int ell = n_ell_theta-1;  // since we start counting from 0, the last multipole is n_ell_theta-1 
    dThetadx[ell] = ck_over_Hp*(Theta[ell-1] - (ell+1)/(k*cosmo->eta_of_x(x))*Theta[ell]) + tau_prime*Theta[ell];
    

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // ...
    // ...
    // No polarization for now
  }

  // SET: Neutrino mutlipoles (Nu_ell)
  if(neutrinos){
    // ...
    // ...
    // No neutrinos for now
  }
  if(!std::isfinite(dPhidx) || !std::isfinite(dThetadx[0])) {
    std::cout << "NaN detected at x=" << x << " k=" << k << std::endl;
    throw "NaN in RHS";
}

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

