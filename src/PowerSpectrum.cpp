#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration();

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
    
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  double z_max    = 40000.0;                          // as mentioned^

  const double delta_z = M_PI / 8.0;                  // idk where this came from tbh
  const int n_z        = int(z_max / delta_z);

  Vector z_array  = Utils::linspace(0.0,z_max, n_z);

  for(size_t i = 0; i < ells.size(); i++){
    int ell = ells[i];

    
    Vector j_ell = Vector(z_array.size(), 0.0);

    for(size_t i = 0; i < z_array.size(); i++){
      j_ell[i] = Utils::j_ell(ell, z_array[i]);
    }

    j_ell_splines[i].create(z_array, j_ell, "j_"+ std::to_string(ell) +"_spline" );

  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array,
    Vector & x_array,
    std::function<double(double,double)> &source_function){

  Utils::StartTiming("lineofsight");
    
  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double eta0     = cosmo->eta_of_x(0.0);

  // Vector x_array  = pert-> x_array; 
  double dx       = x_array[1] - x_array[0]; 

  for(size_t ik = 0; ik < k_array.size(); ik++){

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    // Trapezoidal rule
    //=============================================================================
    double k = k_array[ik];

    Vector LoS_ints(ells.size(), 0.0);
    for(int ix = 0; ix < x_array.size(); ix++){

      double x          = x_array[ix];
      double eta        = cosmo->eta_of_x(x);
      double S_tilde    = source_function(x, k);

      for(size_t ell = 0; ell < ells.size(); ell++){

        double j_ell = j_ell_splines[ell](k*(eta - eta0));
        
        // Trapezoidal rule:  multiply 0.5 at the boundaries
        if (ix == 0 || ix == x_array.size() - 1)
          LoS_ints[ell] += 0.5 * S_tilde * j_ell * dx;

        else
        LoS_ints[ell] += S_tilde * j_ell * dx;

        // Store the result for Source_ell(k) in results[ell][ik]
        result[ell][ik] = LoS_ints[ell];
        
      }
    }

  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(){
  const int nells      = ells.size();
  
  //============================================================================
  // TODO: Make linear spaced k-array with dk ~ 2pi/eta0/N where N >~ 6
  //============================================================================
  
  const double eta0 = cosmo->eta_of_x(0.0);
  const int N       = 6;
  const double dk   = 2.0*M_PI/(eta0/N);
  const int n_k_i   = int((k_max - k_min) / dk);

  Vector k_array = Utils::linspace(k_min, k_max, n_k_i);
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, n_x);



  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array,x_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  for (int i = 0; i < ells.size(); i++) {
    thetaT_ell_of_k_spline[i].create(k_array, thetaT_ell_of_k[i], "thetaT_"+ std::to_string(ells[i])+"ell_of_k_spline");
  }

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================
  if(Constants.polarization){

    // No polarization for now
    // ...
    // ...
    // ...

  }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
    
  //============================================================================
  // TODO: Make log-spaced k-array with n = (kmax-kmin)/ dk points where
  // dk ~ 2.0 * M_PI / eta0 / N with N >~ 32 
  //============================================================================

  const double eta0   = cosmo->eta_of_x(0.0);
  const int N         = 32;
  const double dk     = 2.0*M_PI/(eta0/N);
  int n_k_i           = int((k_max - k_min) / dk);

  Vector log_k_array  = Utils::linspace(log(k_min), log(k_max), n_k_i);
  

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dlog_k
  // Trapezoidal rule
  //============================================================================

  double dlog_k     = log_k_array[1] - log_k_array[0];
  
  Vector result(nells,0.0);

  for(int i = 0; i < nells; i++) {

    int ell         = ells[i];
    double LoS_int  = 0.0;

    for(int ik = 0; ik < n_k_i; ik++) {

      double k      = exp(log_k_array[ik]);
      double P_k    = primordial_power_spectrum(k);

      double f_ell  = f_ell_spline[i](k);
      double g_ell  = g_ell_spline[i](k);

      // Trapezoidal rule
      if (ik == 0 || ik == n_k_i - 1){
        LoS_int += 2.0*M_PI*P_k*f_ell*g_ell * dk;
      }
      else{
        LoS_int += 4.0*M_PI*P_k*f_ell*g_ell * dk;
      }
    }

    result[i] = LoS_int;
  }

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k) const{
  double k_mpc = k / Constants.Mpc; 

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  double Hp             = cosmo->Hp_of_x(x);
  double Phi            = pert->get_Phi(x,k);
  double ck_over_Hp     = (Constants.c * k) / Hp;

  double Delta_M        = 2.0/3.0 * ck_over_Hp * ck_over_Hp * Phi;
  double abs_Delta_M    = std::abs(Delta_M);

  double P_k            = primordial_power_spectrum(k);

  double pofk           = abs_Delta_M * abs_Delta_M * P_k;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

