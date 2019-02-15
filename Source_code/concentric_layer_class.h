//SJL: 10/15
/*
Class structure for a layer for the HERCULES code
*/

#ifndef LAYER_H
#define LAYER_H

class concentric_layer 
{
 public:
  //number of mu points and max legendre degree
  int Nmu, kmax;
  //the rotation rate of layer, density (only of concentric layer), mass, volume
  double omega, rho, M, Vol;
  //equitorial and polar radius
  double a, b;
  //Moment of inertia
  double I;
  //Vector of the angles
  std::vector<double> mu;
  
  //Surface points and updated points
  std::vector<double> xi;
  
  //precomputed integrals
  double Mbar;
  std::vector<double> Js;
  std::vector< std::vector<double> > Nbar_mu;
  std::vector< std::vector<double> > Kbar_mu;
  std::vector< std::vector<double> > Nbar_integrand;
  std::vector< std::vector<double> > Kbar_integrand;
  
  //functions
  void calc_I(std::vector< std::vector<double> >& Plegendre);
  void calc_Mbar(std::vector< std::vector<double> >& Plegendre);
  void calc_M();
  void calc_Js(std::vector< std::vector<double> >& Plegendre);
  void precalc_Nbar(std::vector< std::vector<double> >& Plegendre);
  void precalc_Kbar(std::vector< std::vector<double> >& Plegendre);
  double calc_Pbar(int k, double mu_r);

  double calc_VRegI(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu);
  double calc_VRegII(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu);
  double calc_VRegIV(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu);
  double calc_mu_r(double xi_calc, std::vector<double>::iterator& mu_r_pos, int& mu_r_ind);
  std::vector<double> calc_N(std::vector< std::vector<double> >& Plegendre, double mu_r, std::vector<double>::iterator& mu_r_pos, int& mu_r_ind, double xi_cal);
  std::vector<double> calc_K(std::vector< std::vector<double> >& Plegendre, double mu_r, std::vector<double>::iterator& mu_r_pos, int& mu_r_ind, double xi_cal);

  void write_binary(std::ofstream& file);
  void read_binary(std::ifstream& file);

};


#endif
