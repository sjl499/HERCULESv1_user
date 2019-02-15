//SJL: 11/15 
/*
Function for the planet class in HERCULES
 */

#include "header.h"

////////////////////////////////////////////////////
//Function to precompute the array of Legendre polynomials
void planet::precalc_Plegendre()
{
  for (int k=0; k<(kmax+1); ++k){
    for (int i=0; i<Nmu-1; ++i){
      Plegendre[k][i]=boost::math::legendre_p(2*k,mu[i]);
    }						
    Plegendre[k][Nmu-1]=1.0;
  }
}

////////////////////////////////////////////////////
//Function to initialise a planet with ellipsoidal shells
void planet::initialise_ellipsoids(parameters& params)
{

  std::cout<< "\t Initialising with ellipsoids"  << std::endl;

  //define the mu points to use
  //For now assume linearly spaced in mu
  std::cout<< "\t \t Using evenly spaced mu points"  << std::endl;
  double Deltamu=1.0/(Nmu-1);
  mu[0]=0.0;
  for (int count=1; count<Nmu; ++count){
    mu[count]=mu[count-1]+Deltamu;
  }

  //Now precompute the Legendre polynomials at the mu points
  precalc_Plegendre();

  //Loop over each layer to set raddi, omega
  double Vtot=0.0;
  int counter=Nlayer-1;
  double temp_r=0.0; //Temp array to work outwards in radii
  double Deltar;
  //loop over each material
  for (int i=Nmaterial-1; i>=0; --i){
    Deltar=(layers[counter-material_lay[i]+1].a-temp_r)/material_lay[i]; //radius gap between layers
    //Then loop over each layer in a material
    for (int count=counter; count>=counter-material_lay[i]+1; --count){
      //calcuate the equatorial and polar radii
      temp_r+=Deltar;
      layers[count].a=temp_r;
      layers[count].b=temp_r*aspect;
      
      //set up the raay of mu and xi's
      layers[count].mu=mu;
      for (int i=0; i<Nmu; ++i){
	layers[count].xi[i]=layers[count].b/(sqrt(pow(layers[count].a*mu[i],2.0)+pow(layers[count].b,2.0)-pow(layers[count].b*mu[i],2.0)));
      }
      
      //Calculate the volume
      layers[count].Vol=(4.0*M_PI/3.0)*layers[count].b*pow(layers[count].a,2.0); 
      Vtot+=layers[count].Vol;
      
      //set constant rotation
      layers[count].omega=omega_rot;
      
      //set the Nmu, kmax
      layers[count].Nmu=Nmu;
      layers[count].kmax=kmax;

    }
    counter-=material_lay[i];
  }



  //Assign the densities depending on the Mconc flag
  double temp_M;
  if (params.flag_Mconc==0){
    //For constant density assign the density and reassign the masses

    //Loop over to assign the density and mass to layers
    temp_M=0.0;
    for (int i=0; i<Nlayer; ++i){
      layers[i].rho=ref_rho;
      layers[i].M=layers[i].Vol*ref_rho;
      temp_M+=layers[i].M;
    }
    Mtot=temp_M;
  }
  else if (params.flag_Mconc==1 || params.flag_Mconc==2){
    //Assign the density of layers depending on the mass of each material layer
    //Note in this section M is a proxy for volume
    std::vector<double> rho_materials;
    rho_materials.resize(Nmaterial);
    std::vector<int>::iterator pos;
    double temp_V;
    int inside=0, outside=0;
    for (int count=0; count<Nmaterial; ++count){
      //Find the index of the layer on the outside of the material
      pos=std::find(flag_material.begin(), flag_material.end(), count+1);
      inside=std::distance(flag_material.begin(), pos);
      
      //Now find the mass of all the material in the range for that material
      temp_M=M_materials[count];
      //loop over all layers and add mass if nesc.
      for (int i=0; i<Nlayer; ++i){
	if (flag_material[i]<count){
	  if (inside==layers.size()) {
	    temp_M+=-rho_materials[flag_material[i]]*(layers[outside].Vol);
	  }
	  else {
	    temp_M+=-rho_materials[flag_material[i]]*(layers[outside].Vol-layers[inside].Vol);
	  }
	}
      }
      
      //find the volume within a material
      temp_V=0.0;
      for (int i=outside; i<inside; ++i){
	if (inside==layers.size()) {
	  temp_V+=layers[i].Vol;
	}
	else {
	  temp_V+=layers[i].Vol-layers[inside].Vol;
	}
      }
      
      //Then find the density
      rho_materials[count]=temp_M/temp_V;
      
      //Loop over to assign the density and mass to layers
      for (int i=outside; i<inside; ++i){
	layers[i].rho=rho_materials[count];
	layers[i].M=layers[i].Vol*rho_materials[count];
      }
      
      //Move inside to outside
      outside=inside;
    }
  }
  else {
    std::cerr << "Error in input. Unrecognised Mconc flag." << std::endl;
    std::exit(1);
  }
  
  
  //Loop over the layers to assign real density
  double temp=0.0;
  for (int count=0; count<Nlayer; ++count){
    //sum over the layers to find the 'real' density
    temp+=layers[count].rho;
    real_rho[count]=temp;
  }


  std::cout<< "\t \t Min density " << real_rho[0] << "(SI)" << std::endl;
  std::cout<< "\t \t Max density " << real_rho[Nlayer-1] << "(SI)" << std::endl;

  //Calculate the precalc values and J's, I's etc.
  for (int count=0; count<Nlayer; ++count){
    //precompute values
    layers[count].calc_Mbar(Plegendre);
    layers[count].precalc_Nbar(Plegendre);
    layers[count].precalc_Kbar(Plegendre); 

    //body constants
    layers[count].calc_I(Plegendre);    
    layers[count].calc_Js(Plegendre);

    /*
    //////TEST///////
    //Test the Js
    double temp_l;
    std::cout << "testing Js" <<std::endl;
    for (int k=0; k<(kmax+1); ++k){
      temp_l=pow(aspect, -2)-1.0;
      std::cout << (layers[count].Js[k]-3.0*pow(-1, 1+k)*pow((temp_l/(1.0+temp_l)),k)/((2.0*k+1.0)*(2.0*k+3.0)))/(3.0*pow(-1, 1+k)*pow((temp_l/(1.0+temp_l)),k)/((2.0*k+1.0)*(2.0*k+3.0))) << std::endl;
      //std::cout << k*2 << " " << layers[count].Js[k] << std::endl;
      //std::cout << pow(-1, 1+k) << std::endl;
    }
    std::cout << "testing Js over" <<std::endl;
    */
  }

  //Calculate the angular momentum and mass distribution
  calc_Mout();
  calc_Lout();
  calc_Js();

}

//////////////////////////////////////////////////////
//Function to calculate the gravitational potential at a single point
double planet::calc_V(double calc_r, int ind_mu)
{
  //Initialise the potential
  double V=0.0;
  //loop over all the layers and test which regime 
  for (int i=0; i<Nlayer; ++i){
    if (calc_r>=layers[i].a){
      V+=layers[i].calc_VRegIV(Plegendre, calc_r, ind_mu );
      //std::cout << "RegIV" << std::endl;
    }
    else if (calc_r<=layers[i].b){
      V+=layers[i].calc_VRegI(Plegendre, calc_r, ind_mu );
      //std::cout << "RegI" << std::endl;
    }
    else {
      V+=layers[i].calc_VRegII(Plegendre, calc_r, ind_mu);
      //std::cout << "RegII" << std::endl;
    }
  }

  return V;
}

//////////////////////////////////////////////////////
//Function to calculate the rotational potential at a single point
double planet::calc_Q(double calc_r, int ind_mu, int ind_layer)
{
  //Initialise the potential
  double Q=0.0;

  Q=(pow(layers[ind_layer].xi[ind_mu]*layers[ind_layer].a,2.0)*pow(layers[ind_layer].omega,2)*(1.0-Plegendre[1][ind_mu]))/3.0;

  return Q;
  
}

//////////////////////////////////////////////////////
//Function to calculate the value of the itterative equation
 double planet::calc_X(double calc_r, int ind_mu, int ind_layer)
 {
   double X;
   double temp;
   //gravitational potential
   X=calc_V(calc_r, ind_mu);
   //rotational potential
   temp=calc_Q(calc_r, ind_mu, ind_layer);
   X+=temp;
   //subtract potential needed
   X-=Ulayers[ind_layer];
   
   return X;
 }


///////////////////////////////////////////////////////
//Function to calculate the mass of the planet from the layers
void planet::calc_Mtot()
{
  double temp_Mtot=0.0;

  for (int count=0; count<Nlayer; ++count){
    temp_Mtot+=layers[count].M;
    }
  Mtot=temp_Mtot;
}

///////////////////////////////////////////////////////
//Function to calculate the AM of the planet
void planet::calc_Ltot()
{
  double temp_Ltot=0.0;

  //Loop over all the layers only considering the moment inertia of mass with that layer
  for (int count=0; count<Nlayer-1; ++count){
    for (int j=0; j<count+1; ++j){
      temp_Ltot+=((layers[j].rho/layers[count].rho)*layers[count].I-
		  (layers[j].rho/layers[count+1].rho)*layers[count+1].I)*layers[count].omega;
    }
  }
  //Need to deal with inner layer seperately
  for (int j=0; j<Nlayer; ++j){
    temp_Ltot+=(layers[j].rho/layers[Nlayer-1].rho)*layers[Nlayer-1].I*layers[Nlayer-1].omega;
  }
  
  Ltot=temp_Ltot;
}


///////////////////////////////////////////////////////
//Function to calculate the mass beyond the equitorial radii of each layer
void planet::calc_Mout()
{
  //loop over all the layers
  //Outermost layer is easy
  Mout[0]=0.0;
  for (int count=1; count<Nlayer; ++count){
    Mout[count]=0.0;
    //loop over all layers outside the current layer
    for (int i=0; i<count; ++i){
      //Find the value of mu_c by looping over all mu to find mu above it
      double mu_c=0.0;
      int ind_mu_c=0;
      while ((layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))>(layers[count].a/layers[i].a)){
	++ind_mu_c;
      }
      ind_mu_c-=1;
      //now linearly interpolate to find mu_c and xi_c
      if (ind_mu_c==Nmu-2){
	mu_c=mu[ind_mu_c]+((layers[count].a/layers[i].a)-layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))
	  *(mu[ind_mu_c+1]-mu[ind_mu_c])/
	  (- layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)));
      }
      else {
	mu_c=mu[ind_mu_c]+((layers[count].a/layers[i].a)-layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))
	  *(mu[ind_mu_c+1]-mu[ind_mu_c])/
	  (layers[i].xi[ind_mu_c+1]*sqrt(1-pow(mu[ind_mu_c+1],2.0)) - layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)));
      }
      
      //Integrate Kbar to find the first term
      std::vector<double> integrand_extra(2);
      double Mout_temp=0.0;
      integrand_extra[0]=layers[i].Kbar_integrand[0][ind_mu_c];
      integrand_extra[1]=pow((layers[count].a/layers[i].a)/(sqrt(1-pow(mu_c,2.0))),3.0);
      
      //add the second term and multiply by constants
      Mout_temp=pow(layers[i].a,3.0)*(layers[i].Kbar_mu[0][ind_mu_c]+integrate(integrand_extra, (mu_c-mu[ind_mu_c]), 0, 1));
      Mout_temp-=pow(layers[count].a,3.0)*mu_c/(sqrt(1-pow(mu_c,2.0)));
      Mout[count]+=Mout_temp*layers[i].rho*4.0*M_PI/3.0;
    }

    //std::cout << "Mass outside " << layers[count].a << "\t" << Mout[count] << std::endl;
  }
}


///////////////////////////////////////////////////////
//Function to calculate the AM beyond the equitorial radii of each layer
void planet::calc_Lout()
{
  //loop over all the layers
  //Outermost layer is easy
  Lout[0]=0.0;
  for (int count=1; count<Nlayer; ++count){
    Lout[count]=0.0;
    //loop over all layers outside the current layer
    for (int i=0; i<count; ++i){
      //Find the value of mu_c by looping over all mu to find mu above it
      //Also calculate the integrand for each point 
      double mu_c=0.0;
      int ind_mu_c=0;
      std::vector<double> integrand;
      integrand.resize(Nmu);
      while ((layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))>(layers[count].a/layers[i].a)){
	integrand[ind_mu_c]=(1-pow(mu[ind_mu_c],2.0))*
	  (pow(layers[i].xi[ind_mu_c]*layers[i].a,5.0)-
	   pow(std::max(layers[i+1].xi[ind_mu_c]*layers[i+1].a,(layers[count].a/sqrt(1-pow(mu[ind_mu_c],2.0)))),5.0));
	++ind_mu_c;
      }
      ind_mu_c-=1;
      //now linearly interpolate to find mu_c
      if (ind_mu_c==Nmu-2){
	mu_c=mu[ind_mu_c]+((layers[count].a/layers[i].a)-layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))
	  *(mu[ind_mu_c+1]-mu[ind_mu_c])/
	  (-layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)));
      }
      else {
	mu_c=mu[ind_mu_c]+((layers[count].a/layers[i].a)-layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)))
	  *(mu[ind_mu_c+1]-mu[ind_mu_c])/
	  (layers[i].xi[ind_mu_c+1]*sqrt(1-pow(mu[ind_mu_c+1],2.0)) - layers[i].xi[ind_mu_c]*sqrt(1-pow(mu[ind_mu_c],2.0)));
      }
      
      //Integrate to find the first term
      std::vector<double> integrand_extra(2);
      double Lout_temp=0.0;
      integrand_extra[0]=integrand[ind_mu_c];
      integrand_extra[1]=0.0;
      
      //Integrate and multiply by constants
      Lout_temp=(integrate(integrand, 1.0/(Nmu-1), 0, ind_mu_c)+integrate(integrand_extra, (mu_c-mu[ind_mu_c]), 0, 1));
      Lout[count]+=Lout_temp*real_rho[i]*layers[i].omega*4.0*M_PI/5.0;

    }

  }
  
}


///////////////////////////////////////////////////////
//Function to calculate the gravitational moments of the whole body
void planet::calc_Js()
{
  for (int k=0; k<kmax+1; ++k){
    Js[k]=0.0;
    for (int i=0; i<Nlayer; ++i){
      Js[k]+=layers[i].M*pow(layers[i].a,2.0*k)*layers[i].Js[k];
    }
    Js[k]/=(Mtot*pow(amax, 2.0*k));
  }
}

///////////////////////////////////////////////////////////////
//Function to write planer to binary
void planet::write_binary(std::ofstream& file)
{
  //write out the size constraining parameters
  file.write ((char*)&Nmu, sizeof (int));
  file.write ((char*)&Nlayer, sizeof (int));
  file.write ((char*)&Nmaterial, sizeof (int));
  file.write ((char*)&kmax, sizeof (int));

  //write single parameter values
  file.write ((char*)&Mtot, sizeof(double));
  file.write ((char*)&Ltot, sizeof(double));
  file.write ((char*)&omega_rot, sizeof(double));
  file.write ((char*)&pmin, sizeof(double));
  file.write ((char*)&amax, sizeof(double));
  file.write ((char*)&aspect, sizeof(double));
  file.write ((char*)&ref_rho, sizeof(double));
  file.write ((char*)&Mtot_tar, sizeof(double));
  file.write ((char*)&Ltot_tar, sizeof(double));
  file.write ((char*)&Ucore, sizeof(double));
  file.write ((char*)&pcore, sizeof(double));

  //write the vectors
  file.write (reinterpret_cast<const char*>(&real_rho[0]), Nlayer*sizeof (real_rho[0]));
  file.write (reinterpret_cast<const char*>(&press[0]), Nlayer*sizeof (press[0]));
  file.write (reinterpret_cast<const char*>(&dpdr[0]), Nlayer*sizeof(dpdr[0]));
  file.write (reinterpret_cast<const char*>(&mu[0]), Nmu*sizeof(mu[0]));
  file.write (reinterpret_cast<const char*>(&Ulayers[0]), Nlayer*sizeof(Ulayers[0]));
  file.write (reinterpret_cast<const char*>(&flag_material[0]), Nlayer*sizeof(flag_material[0]));
  file.write (reinterpret_cast<const char*>(&M_materials[0]), Nmaterial*sizeof(M_materials[0]));
  file.write (reinterpret_cast<const char*>(&Mout[0]), Nlayer*sizeof(Mout[0]));
  file.write (reinterpret_cast<const char*>(&Lout[0]), Nlayer*sizeof(Lout[0]));
  file.write (reinterpret_cast<const char*>(&Js[0]), (kmax+1)*sizeof(Js[0]));

  //loop over all the layers and print
  for (int count=0; count<Nlayer; ++count){
    layers[count].write_binary(file);
  }

  //loop over all the material and print the input files
  for (int count=0; count<Nmaterial; ++count){
    materials[count].write_binary(file);
  }
}



///////////////////////////////////////////////////////////////
//Function to read planet from binary
void planet::read_binary(std::ifstream& file)
{
  //read in the size defining parameters
  file.read ((char*)&Nmu, sizeof (int));
  file.read ((char*)&Nlayer, sizeof (int));
  file.read ((char*)&Nmaterial, sizeof (int));
  file.read ((char*)&kmax, sizeof (int));

  //read single parameter values
  file.read ((char*)&Mtot, sizeof(double));
  file.read ((char*)&Ltot, sizeof(double));
  file.read ((char*)&omega_rot, sizeof(double));
  file.read ((char*)&pmin, sizeof(double));
  file.read ((char*)&amax, sizeof(double));
  file.read ((char*)&aspect, sizeof(double));
  file.read ((char*)&ref_rho, sizeof(double));
  file.read ((char*)&Mtot_tar, sizeof(double));
  file.read ((char*)&Ltot_tar, sizeof(double));
  file.read ((char*)&Ucore, sizeof(double));
  file.read ((char*)&pcore, sizeof(double));

  
 //resize and read the vectors
  real_rho.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&real_rho[0]), Nlayer*sizeof(double)); 
  press.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&press[0]), Nlayer*sizeof(double));
  dpdr.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&dpdr[0]), Nlayer*sizeof(double));
  mu.resize(Nmu);
  file.read (reinterpret_cast<char*>(&mu[0]), Nmu*sizeof(double));
  Ulayers.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&Ulayers[0]), Nlayer*sizeof(double));
  flag_material.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&flag_material[0]), Nlayer*sizeof(int));
  M_materials.resize(Nmaterial);
  file.read (reinterpret_cast<char*>(&M_materials[0]), Nmaterial*sizeof(double));
  Mout.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&Mout[0]), Nlayer*sizeof(double));
  Lout.resize(Nlayer);
  file.read (reinterpret_cast<char*>(&Lout[0]), Nlayer*sizeof(double));
  Js.resize(kmax+1);
  file.read (reinterpret_cast<char*>(&Js[0]), (kmax+1)*sizeof(double));


  //loop over all the layers and read
  layers.resize(Nlayer);
  for (int count=0; count<Nlayer; ++count){
    layers[count].read_binary(file);
  }

  //loop over all the material and read the input files
  materials.resize(Nmaterial);
  for (int count=0; count<Nmaterial; ++count){
    materials[count].read_binary(file);
  }

  //define the length of the Plegendre
  Plegendre.resize(kmax+1);
  for (int k=0; k<(kmax+1); ++k){
    Plegendre[k].resize(Nmu);
  }	
  
}
