//SJL: 11/15
/*
Functions for the class structure for a layer for the HERCULES code
*/
#include "header.h"

/////////////////////////////////////////////////////
//Function to calculate the moment of inertia of a layer
void concentric_layer::calc_I(std::vector< std::vector<double> >& Plegendre)
{
std::vector<double> integrand;
integrand.resize(Nmu);
for (int i=0; i<Nmu-1; ++i){
  integrand[i]=(8.0*M_PI/15.0)*rho*pow((xi[i]*a),5.0)*(1-Plegendre[1][i]);
  }
integrand[Nmu-1]=0.0; //in case of rounding error do last term seperately

 I=integrate(integrand, 1.0/(Nmu-1), 0, Nmu-1);

}


/////////////////////////////////////////////////////
//Function to calculate the moment of inertia of a layer
void concentric_layer::calc_Mbar(std::vector< std::vector<double> >& Plegendre)
{
  std::vector<double> integrand;
  integrand.resize(Nmu);
  for (int i=0; i<Nmu; ++i){
    integrand[i]=pow(xi[i],3.0);
  }
  
  Mbar=integrate(integrand, 1.0/(Nmu-1), 0, Nmu-1);
  
}

/////////////////////////////////////////////////////
//Function to calculate the mass and volume
void concentric_layer::calc_M()
{
  Vol=Mbar*(4.0*M_PI*pow(a,3)/3.0);
  M=Vol*rho;

}


/////////////////////////////////////////////////////
//Function to precalcualte the Nbar integrand
void concentric_layer::precalc_Nbar(std::vector< std::vector<double> >& Plegendre)
{
  // for k=0
  for (int i=0; i<Nmu; ++i){
    Nbar_integrand[0][i]=pow(xi[i],2.0);
  }

  // for k=1
  for (int i=0; i<Nmu; ++i){
    Nbar_integrand[1][i]=log(xi[i])*Plegendre[1][i];
  }

  // for k>1
  for (int k=2; k<(kmax+1); ++k){
    for (int i=0; i<Nmu; ++i){
      Nbar_integrand[k][i]=pow(xi[i],(2.0-2.0*k))*Plegendre[k][i];
  }
  }

  //Calcualte the Nbar_mu values (i.e. the integral to each mu)
  int ind_mid=(int)(Nmu/2.0);
  for (int k=0; k<(kmax+1); ++k){
    //set the zero value and set the integral from 0 to 1
    Nbar_mu[k][0]=0.0;
    Nbar_mu[k][Nmu-1]=integrate(Nbar_integrand[k], 1.0/(Nmu-1), 0, Nmu-1);
    //for the first half of points do the reverse integral
    for (int i=1; i<ind_mid; ++i){
      Nbar_mu[k][i]=Nbar_mu[k][Nmu-1]-integrate(Nbar_integrand[k], 1.0/(Nmu-1), i, Nmu-1);
    }
    for (int i=ind_mid; i<Nmu-1; ++i){
      Nbar_mu[k][i]=integrate(Nbar_integrand[k], 1.0/(Nmu-1), 0, i);
    }
  }
 
}

/////////////////////////////////////////////////////
//Function to precalcualte the Kbar integrand
void concentric_layer::precalc_Kbar(std::vector< std::vector<double> >& Plegendre)
{
  // for k=0
  for (int i=0; i<Nmu; ++i){
    Kbar_integrand[0][i]=pow(xi[i],3.0);
  }

  // for k>1
  for (int k=1; k<(kmax+1); ++k){
    for (int i=0; i<Nmu; ++i){
      Kbar_integrand[k][i]=pow(xi[i],(3.0+2.0*k))*Plegendre[k][i];
    }
  }  

  //Calcualte the Kbar_mu values (i.e. the integral to each mu)
  int ind_mid=(int)(Nmu/2.0);
  for (int k=0; k<(kmax+1); ++k){
    //set the zero value and set the integral from 0 to 1
    Kbar_mu[k][0]=0.0;
    Kbar_mu[k][Nmu-1]=integrate(Kbar_integrand[k], 1.0/(Nmu-1), 0, Nmu-1);
    //for the first half of points do the reverse integral
    for (int i=1; i<ind_mid; ++i){
      Kbar_mu[k][i]=Kbar_mu[k][Nmu-1]-integrate(Kbar_integrand[k], 1.0/(Nmu-1), i, Nmu-1);
    }
    for (int i=ind_mid; i<Nmu-1; ++i){
      Kbar_mu[k][i]=integrate(Kbar_integrand[k], 1.0/(Nmu-1), 0, i);
    }
  }


}


//////////////////////////////////////////////////////
//Function to calculate Pbar 
double concentric_layer::calc_Pbar(int k, double mu_r)
{
  double Pbar;
  if (k==0){
    Pbar=1.0;
  }
  else {
  //Calculate Pbar as the sum of two Legendre polynomials
    Pbar=(boost::math::legendre_p(2*k+1,mu_r)-boost::math::legendre_p(2*k-1,mu_r))/(4.0*k+1.0);
  }

  return Pbar;
}


/////////////////////////////////////////////////////
//Function to calculate the moment of inertia of a layer
//Mbar must be up to date to call this function
void concentric_layer::calc_Js(std::vector< std::vector<double> >& Plegendre)
{
  std::vector<double> integrand;
  integrand.resize(Nmu);

  //loop over all the 2k values to find the corresponding J's
  for (int k=0; k<Js.size(); ++k)
    {
      //Find the integrand at each mu
      for (int i=0; i<Nmu; ++i)
	{
	  integrand[i]=pow(xi[i],(2*k+3.0))*Plegendre[k][i];
	}
      Js[k]=(-3.0/(Mbar*(2*k+3.0)))*integrate(integrand, 1.0/(Nmu-1), 0, Nmu-1);
    }

  /*
  //////TESTING///////////////
  //Analytical for Maclauren spheroid
  double testing;
  double l;
  l=pow(a,2.0)/pow(b,2.0)-1.0;
  for (int k=0; k<Js.size(); ++k)
    {
      testing=3.0*pow(-1.0,1+k)*pow((l/(1.0+l)),k)/(2.0*k+1.0)/(2.0*k+3.0);
      std::cout << 2*k << "\t" << (Js[k]-testing)/testing << std::endl; 
    }
  */

}

/////////////////////////////////////////////////////
//Function to calculate the contribution from the layer in the case that r<=b
double concentric_layer::calc_VRegI(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu)
{
  double xi_calc=calc_r/a;
  double V;
  double N_temp;
  
  //Start from the zeroeth order term
  N_temp=(3.0/2.0)*(Nbar_mu[0][Nmu-1]-pow(xi_calc,2.0))/Mbar;
  V=(N_temp)/a;

  //k=1
  N_temp=-3.0*Nbar_mu[1][Nmu-1]/Mbar;
  V+=-pow(xi_calc,2.0)*(N_temp)*Plegendre[1][ind_mu]/a;
  
  //Loop over all other terms
  for (int k=2; k<(kmax+1); ++k){
    N_temp=(3.0/(2.0*k-2.0))*Nbar_mu[k][Nmu-1]/Mbar;
    V+=-pow(xi_calc,2.0*k)*N_temp*Plegendre[k][ind_mu]/a;
    }
  V=V*M;
  //Add the internal gravity term
  V+=(4.0*M_PI/3.0)*rho*pow(calc_r,2.0);

  return V*CONST::G;

}

/////////////////////////////////////////////////////
//Function to calculate the gravity in the regime b<r<a
double concentric_layer::calc_VRegII(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu)
{
  double xi_calc=calc_r/a;
  double V;
  double mu_r;
  std::vector<double>::iterator mu_r_pos;
  int mu_r_ind;

  //Calculate the mu_r and get the position of the point below it
  mu_r=calc_mu_r(xi_calc, mu_r_pos, mu_r_ind);

  //Calculate the N and K
  std::vector<double> Ni, Ki;
  Ni=calc_N(Plegendre, mu_r, mu_r_pos, mu_r_ind, xi_calc);
  Ki=calc_K(Plegendre, mu_r, mu_r_pos, mu_r_ind, xi_calc);

  //Add the first terms
  V=1.0-Ki[0]+xi_calc*Ni[0];

  //Loop over all the other terms
  for (int k=1; k<(kmax+1); ++k){
    V+=-(pow(xi_calc, 2*k+1.0)*Ni[k]+pow(xi_calc, -2.0*k)*(Js[k]-Ki[k]))*Plegendre[k][ind_mu];
  }
  
  return V*CONST::G*M/calc_r;
}


/////////////////////////////////////////////////////
//Function to calculate the gravity in the regime r>=a
double concentric_layer::calc_VRegIV(std::vector< std::vector<double> >& Plegendre, double calc_r, int ind_mu)
{
  double xi_calc=calc_r/a;
  double V=1.0;

  //Loop over each terms
  for (int k=1; k<(kmax+1); ++k){
    V+=-pow(xi_calc, -2.0*k)*Js[k]*Plegendre[k][ind_mu];
  }

  return V*CONST::G*M/calc_r;

}

/////////////////////////////////////////////////////
//Function to calc mu_r of the surface
double concentric_layer::calc_mu_r(double xi_calc, std::vector<double>::iterator& pos, int& ind)
{
  double mu_r=0.0;

  //SJL 11/15
  //Initially find mu_r by linear interpolation

  //find the elements below mu_r
  pos=std::find_if(xi.begin(), xi.end(), std::bind2nd(std::less<double>(),xi_calc));
  pos=pos-1;
  ind=std::distance(xi.begin(), pos);
  //and linearly interp
  mu_r=mu[ind]+(xi_calc-xi[ind])*(mu[ind+1]-mu[ind])/(xi[ind+1] - xi[ind]);

  return mu_r;
}

/////////////////////////////////////////////////////
//Function to calc the N needed for gravity calculations
std::vector<double> concentric_layer::calc_N(std::vector< std::vector<double> >& Plegendre, double mu_r, std::vector<double>::iterator& mu_r_pos, int& mu_r_ind, double xi_calc)
{
  std::vector<double> Ni;
  Ni.resize(kmax+1);
  std::vector<double> integrand_extra(2);

  ////////////
  // for k=0
  //Calculate the integrand for mu_r_ind to mu_r
  integrand_extra[0]=Nbar_integrand[0][mu_r_ind];
  integrand_extra[1]=pow(xi_calc,2.0);

  Ni[0]=Nbar_mu[0][mu_r_ind]+integrate(integrand_extra, (mu_r-mu[mu_r_ind]), 0, 1);
  Ni[0]=1.5*(Ni[0]-mu_r*integrand_extra[1])/Mbar;

  // for k=1
  integrand_extra[0]=Nbar_integrand[1][mu_r_ind];
  integrand_extra[1]=log(xi_calc)*boost::math::legendre_p(2,mu_r);

  Ni[1]=Nbar_mu[1][mu_r_ind]+integrate(integrand_extra, (mu_r-mu[mu_r_ind]), 0, 1);
  Ni[1]=-3.0*(Ni[1]-calc_Pbar(1, mu_r)*log(xi_calc))/Mbar;

  // for k>1
  for (int k=2; k<(kmax+1); ++k){
    integrand_extra[0]=Nbar_integrand[k][mu_r_ind];
    integrand_extra[1]=pow(xi_calc, 2.0-2.0*k)*boost::math::legendre_p(2*k,mu_r);

    Ni[k]=Nbar_mu[k][mu_r_ind]+integrate(integrand_extra, (mu_r-mu[mu_r_ind]), 0, 1);
    Ni[k]=(-3.0/(2.0*k-2.0))*(calc_Pbar(k, mu_r)*pow(xi_calc, 2.0-2.0*k)-Ni[k])/Mbar;
  }

  return Ni;
}


/////////////////////////////////////////////////////
//Function to calc the K needed for gravity calculations
std::vector<double> concentric_layer::calc_K(std::vector< std::vector<double> >& Plegendre, double mu_r, std::vector<double>::iterator& mu_r_pos, int& mu_r_ind, double xi_calc)
{
  std::vector<double> Ki;
  Ki.resize(kmax+1);
  std::vector<double> integrand_extra(2);

  ////////////
  // for k=0
  //Calculate the integrand for mu_r_ind to mu_r
  integrand_extra[0]=Kbar_integrand[0][mu_r_ind];
  integrand_extra[1]=pow(xi_calc,3.0);

  Ki[0]=Kbar_mu[0][mu_r_ind]+integrate(integrand_extra, (mu_r-mu[mu_r_ind]), 0, 1);
  Ki[0]=(Ki[0]-mu_r*integrand_extra[1])/Mbar;

  // for k>0
  for (int k=1; k<(kmax+1); ++k){
    integrand_extra[0]=Kbar_integrand[k][mu_r_ind];
    integrand_extra[1]=pow(xi_calc, 2*k+3)*boost::math::legendre_p(2*k,mu_r);

    Ki[k]=Kbar_mu[k][mu_r_ind]+integrate(integrand_extra, (mu_r-mu[mu_r_ind]), 0, 1);
    Ki[k]=(-3.0/(2.0*k+3.0))*(Ki[k]-calc_Pbar(k, mu_r)*pow(xi_calc, 2.0*k+3.0))/Mbar;
  }
  

  return Ki;
}



///////////////////////////////////////////////////////////////
//Function to write planer to binary
void concentric_layer::write_binary(std::ofstream& file)
{
  //write out the size constraining parameters
  file.write ((char*)&Nmu, sizeof (int));
  file.write ((char*)&kmax, sizeof (int));

  //write single parameter values
  file.write ((char*)&omega, sizeof(double));
  file.write ((char*)&rho, sizeof(double));
  file.write ((char*)&M, sizeof(double));
  file.write ((char*)&Vol, sizeof(double));
  file.write ((char*)&a, sizeof(double));
  file.write ((char*)&b, sizeof(double));
  file.write ((char*)&I, sizeof(double));

  //write the vector of mu's and xi's
  file.write (reinterpret_cast<const char*>(&mu[0]), Nmu*sizeof(double));
  file.write (reinterpret_cast<const char*>(&xi[0]), Nmu*sizeof (double));

  //write the Js
  file.write (reinterpret_cast<const char*>(&Js[0]), (kmax+1)*sizeof(double));
  
}



///////////////////////////////////////////////////////////////
//Function to read planet from binary
void concentric_layer::read_binary(std::ifstream& file)
{
  //read in the size defining parameters
  file.read ((char*)&Nmu, sizeof (int));
  file.read ((char*)&kmax, sizeof (int));

  //write single parameter values
  file.read ((char*)&omega, sizeof(double));
  file.read ((char*)&rho, sizeof(double));
  file.read ((char*)&M, sizeof(double));
  file.read ((char*)&Vol, sizeof(double));
  file.read ((char*)&a, sizeof(double));
  file.read ((char*)&b, sizeof(double));
  file.read ((char*)&I, sizeof(double));

  //resize and read the vector of mu's and xi's
  mu.resize(Nmu);
  file.read (reinterpret_cast<char*>(&mu[0]), Nmu*sizeof(double));
  xi.resize(Nmu);
  file.read (reinterpret_cast<char*>(&xi[0]), Nmu*sizeof(double));

  //read in the Js
  Js.resize(kmax+1);
  file.read (reinterpret_cast<char*>(&Js[0]), (kmax+1)*sizeof(double));

  //define the size of the other arrays
  Nbar_mu.resize(kmax+1);
  Kbar_mu.resize(kmax+1);
  Nbar_integrand.resize(kmax+1);
  Kbar_integrand.resize(kmax+1);
  for (int k=0; k<(kmax+1); ++k){
    Nbar_integrand[k].resize(Nmu);
    Kbar_integrand[k].resize(Nmu);
    Nbar_mu[k].resize(Nmu);
    Kbar_mu[k].resize(Nmu);
  }	

}
