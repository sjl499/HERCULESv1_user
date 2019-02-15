//SJL: 10/15
/*
Class structure for a planet for the HERCULES code
*/

#ifndef PLANET_H
#define PLANET_H

class planet 
{
 public:
  //Number of mu points, number of layers, number of different materials
  int Nmu, Nlayer, Nmaterial;
  //max legendre degree to calculate to
  int kmax;
  //mass, angular momentum of planet
  double Mtot, Ltot;
  //rotation rate of corotating region, pressure of outer layer
  double omega_rot, pmin;
  //equatorial radius, aspect ratio (b/a)
  double amax, aspect;
  //Reference density for constant density layers
  double ref_rho;
  //Target mass and AM
  double Mtot_tar, Ltot_tar;
  //core properties
  double Ucore, pcore;

  //Vector of the real density, pressure, dpdr of the layers, angles and potential layers
  std::vector<double> real_rho, press, dpdr, mu, Ulayers;
  //vector of flags for which layers are which materials
  std::vector<int> flag_material;

  //Array of the masses of each of the materials
  std::vector<double> M_materials;
  std::vector<int> material_lay;

  //Vector containing all the layers
  std::vector<concentric_layer> layers;
  //Array of the EOS for the different materials
  std::vector<EOS> materials;

  //Array of precomputed Legendre polynomials
  std::vector< std::vector<double> > Plegendre;
  
  //Array of the Mass and AM  outside the equitorial radii of each layer
  std::vector<double> Mout, Lout;

  //Array of the whole planet Jn's
  std::vector<double> Js;
  
  /////////////////////////////////////////
  //functions
  void precalc_Plegendre();
  void initialise_ellipsoids(parameters& params);
  double calc_V(double calc_r, int ind_mu);
  double calc_Q(double calc_r, int ind_mu, int ind_layer);
  double calc_X(double calc_r, int ind_mu, int ind_layer);
  void calc_Mtot();
  void calc_Ltot();
  void calc_Mout();
  void calc_Lout();
  void calc_Js();
  
  void write_binary(std::ofstream& file);
  void read_binary(std::ifstream& file);
  
};


#endif
