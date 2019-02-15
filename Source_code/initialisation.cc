//SJL: 10/15
/*
Function that sets up the HERCULES model before itteration
*/

#include "header.h"

void initialisation(planet& p, parameters& param)
{
  std::cout << "Intialising " <<  std::endl;
  //Read in the input data from the input file
  //This also reads in the EOS data
  read_input(p, param);
  
  //Initialise the planet depending on start flag
  ////////////////////////////////////////////////////////////////////////
  //1) Start from concentric spheroids
  ////////////////////////////////////////////////////////////////////////
  if (param.flag_start==0){
    //Initialise the major properties with constant rotation
    p.initialise_ellipsoids(param);
    //Test as to whether to conserve AM
    if (param.flag_Lconc==0)
      {
	
	for (int count=p.Nlayer-1; count>=0; --count){
	  p.layers[count].omega=p.omega_rot;
	}
	std::cout << "\t \t Constant rotation with omega = " << p.omega_rot << std::endl;
      }
    else if (param.flag_Lconc==1)
      {
	//Here we assume that the whole body is corotating
	double temp=0.0;
	for (int count=p.Nlayer-1; count>=0; --count){
	  temp+=p.layers[count].I;
	}
	p.omega_rot=p.Ltot/temp;
	for (int count=p.Nlayer-1; count>=0; --count){
	  p.layers[count].omega=p.omega_rot;
	}
	std::cout << "\t \t AM conservation selected"  << std::endl;
      }
    else {
      // Print an error and exit
      std::cerr << "Error in input. Unrecognised L conservation flag" << std::endl;
      std::exit(1);
    }
    
  }
  ////////////////////////////////////////////////////////////////////////
  //1) Start from an output file
  ////////////////////////////////////////////////////////////////////////
  else if (param.flag_start==1){
    // Print the initialisation
    std::cout << "\t Initialisation by reading start file " << param.start_file <<  std::endl;

    //Create new planet and paramter file to read in file
    planet read_p;
    parameters read_param;
    read_binary_structure(read_p, read_param, param.start_file);

    //Test to see if the number of layers and mu are the same
    if ((read_p.Nlayer==p.Nlayer) && (read_p.Nmu==p.Nmu)){
      //attribute the necessary components to the structure
      p.omega_rot=read_p.omega_rot;
      p.amax=read_p.amax;
      p.aspect=read_p.aspect;
      p.real_rho=read_p.real_rho;
      p.press=read_p.press;
      p.dpdr=read_p.dpdr;
      p.mu=read_p.mu;


      //for the layers assign the whole layer and correct kmax
      p.layers=read_p.layers;
      for (int count=0; count<p.Nlayer; ++count){
	p.layers[count].kmax=p.kmax;
      }

      //now precompute the various parameters
      p.precalc_Plegendre();
      for (int count=0; count<p.Nlayer; ++count){
	//precompute values
	p.layers[count].calc_Mbar(p.Plegendre);
	p.layers[count].calc_M();
	p.layers[count].precalc_Nbar(p.Plegendre);
	p.layers[count].precalc_Kbar(p.Plegendre); 
	
	//body constants
	p.layers[count].calc_I(p.Plegendre);    
	p.layers[count].calc_Js(p.Plegendre);
      }

      p.calc_Mtot();
      p.calc_Ltot();
      

    }

    else {
      std::cerr << "Error in input. Start file not compatable. Needs to be coded." << std::endl;
      std::exit(1);
    }
  }
  else {
    // Print an error and exit
    std::cerr << "Error in input. Unrecognised start flag" << std::endl;
    std::exit(1);
  }

  for (int count=0; count<p.Nlayer; ++count){
    p.Ulayers[count]=p.calc_V(p.layers[count].a, 0);
    p.Ulayers[count]+=p.calc_Q(p.layers[count].a, 0, count);
  }
  
  std::cout << "\t Initialisation complete" <<std::endl;
}
