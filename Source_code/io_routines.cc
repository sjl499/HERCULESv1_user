//SJL 11/15
/*
Functions to write data from HERCULES code
 */

#include "header.h"

////////////////////////////////////////////////////////////////////
/*
Function that reads in the input from the HERCULES input file.
This function also defines the length of the vectors
*/
void read_input(planet& p, parameters& param)
{
  std::cout << "\t Reading input file " <<  std::endl;
  //define the input file
  std::ifstream input("HERCULES_input.txt");

  //Exit if couldn't open the input file stream for reading
  if (!input){
    // Print an error and exit
    std::cerr << "Input file 'HERCULES_input.txt' unreadable" << std::endl;
    std::exit(1);
  }
  
  //Define a throw away string and the delimiter dividing input sections
  std::string temp("temp");
  std::string delimiter("#");
  
  
  //read header lines
  while (delimiter.compare(temp.substr(0,1)) != 0)
    {
      input >> temp;
    }
  
  //////////////////////////////////////////////////
  //Read in run properties
  input >> temp >> temp >> param.run_name;
  getline(input, temp);
  input >> temp >> temp >> param.nint_max;
  getline(input, temp);
  input >> temp >> temp >> param.toll;
  getline(input, temp);
  input >> temp >> temp >> param.xi_nint_max;
  getline(input, temp);
  input >> temp >> temp >> param.xi_toll;
  getline(input, temp);
  input >> temp >> temp >> param.dxi;
  getline(input, temp);
  input >> temp >> temp >> param.flag_start;
  getline(input, temp);
  input >> temp >> temp >> param.start_file;
  getline(input, temp);
  input >> temp >> temp >> param.flag_Mconc;
  getline(input, temp);
  input >> temp >> temp >> param.flag_Lconc;
  getline(input, temp);

  //Read the vector of rotational profile paramters
  param.omega_param.resize(3);
  double temp_doub; //define temporary arrays to read vectors 
  input >> temp >> temp >> temp; //read input definition
  //Loop to read in each item of vector
  for (int count=0; count < 3; ++count) 
    {
      input >> temp_doub >> temp;
      param.omega_param.at(count)=temp_doub;
    }
  getline(input, temp);

  //Continue reading run properties
  input >> temp >> temp >> param.flag_iter_print;
  getline(input, temp);

  //std::cout << param.xi_nint_max << "\t" << param.xi_toll << "\t" <<param.dxi << "\t" << param.flag_start << std::endl;
  
  //more header files
  while (delimiter.compare(temp.substr(0,1)) != 0){
    input >> temp;
  }
  /////////////////////////////////////////////////
  //Read in planet properties
  input >> temp >> temp >> p.aspect;
  getline(input, temp);
  input >> temp >> temp >> p.Nlayer;
  getline(input, temp);
  input >> temp >> temp >> p.Nmaterial;
  getline(input, temp);
  
  //Read the vector of layer numbers
  int temp_int; //define temporary arrays to read vectors 
  p.material_lay.resize(p.Nmaterial);
  input >> temp >> temp >> temp; //read input definition
  //Loop to read in each item of vector
  for (int count=0; count < p.Nmaterial; ++count) 
    {
      input >> temp_int >> temp;
      p.material_lay.at(count)=temp_int;
    }
  getline(input, temp);
  
  //Test whether the number of layers defined correctly
  temp_int = std::accumulate(p.material_lay.begin(), p.material_lay.end(), 0);
  if (temp_int!=p.Nlayer){
    // Print an error and exit
    std::cerr << "Error in input. Number of layers inconsistent" << std::endl;
    std::cerr << "Nlayers given " << p.Nlayer << ". Nmaterial_lay total " << temp_int << std::endl;
    std::exit(1);
  }
  
  //Define the array of material flags by looping over materials
  p.flag_material.resize(p.Nlayer);
  int ind=0;
  for (int i=0; i < p.Nmaterial; ++i) {
    for (int j=0; j < p.material_lay.at(i); ++j) {
      p.flag_material.at(ind)=i;
      ++ind;
    }
  }

  //Read the vector of max layer radii
  std::vector<double> material_amax;
  material_amax.resize(p.Nmaterial);
  input >> temp >> temp >> temp; //read input definition
  //Loop to read in each item of vector
  for (int count=0; count < p.Nmaterial; ++count) 
    {
      input >> temp_doub >> temp;
      material_amax.at(count)=temp_doub;
    }
  getline(input, temp);
  p.amax=material_amax.at(0);

  
  //Read in the EOS files for materials and intialise EOS classes
  std::string file_name; //define temporary arrays to read vectors 
  p.materials.resize(p.Nmaterial);
  input >> temp >> temp >> temp; //read input definition
  //Loop to read in each item of vector
  for (int count=0; count < p.Nmaterial; ++count) 
    {
      input >> file_name  >> temp;
      p.materials.at(count).read_EOSfile(file_name);
    }
  getline(input, temp);
  
  //Read the vector of material masses
  p.M_materials.resize(p.Nmaterial);
  input >> temp >> temp >> temp; //read input definition
  //Loop to read in each item of vector
  for (int count=0; count < p.Nmaterial; ++count) 
    {
      input >> temp_doub >> temp;
      p.M_materials.at(count)=temp_doub;
    }
  getline(input, temp);

  //Read in the refference density
  input >> temp >> temp >> p.ref_rho;
  getline(input, temp);

  //sum over the material layers to find the total mass
  //Note this is only provisional
  p.Mtot=VecDoub::VecSum(p.M_materials);
  p.Mtot_tar=p.Mtot;
  
  
  //Continue reading planet properties
  input >> temp >> temp >> p.kmax;
  getline(input, temp);
  input >> temp >> temp >> p.Nmu;
  getline(input, temp);
  
  
  //test that the number of angular points is sufficient
  if (p.Nmu<(2.0*p.kmax)){
      // Print an error and exit
    std::cerr << "Error in input. Too few angular points" << std::endl;
    std::cerr << "Nmu given " << p.Nmu << ". Min Nmu needed " << 2.0*p.kmax << std::endl;
    std::exit(1);
  }
  
  //continue reading planet properties
  input >> temp >> temp >> p.omega_rot;
  getline(input, temp);
  input >> temp >> temp >> p.Ltot;
  getline(input, temp);
  input >> temp >> temp >> p.pmin;

  //set the target Ltot
  p.Ltot_tar=p.Ltot;


  ////////////////////////////////////////////////
  //Set up the full run name
  std::ostringstream name;
  name << param.run_name << "_L" << p.Ltot_tar/CONST::LEM <<  "_N"  << p.Nlayer << "_Nm" << p.Nmu << "_k" << p.kmax << "_f" << param.flag_start << param.flag_Mconc << param.flag_Lconc << "_p" << p.pmin/1E5 << "_l" << param.omega_param[0] << "_" << param.omega_param[1] << "_" << param.omega_param[2]  ;
  param.run_name=name.str();
  std::cout << "\t Run name: " << param.run_name << std::endl;
  //std::cout << param.run_name.length() << std::endl;
  
  ////////////////////////////////////////////////
  //Set up the rest of the planet etc.
  //Redefine the length of the vectors we haven't set yet
  p.layers.resize(p.Nlayer);
  p.real_rho.resize(p.Nlayer);
  p.press.resize(p.Nlayer);
  p.dpdr.resize(p.Nlayer);
  p.mu.resize(p.Nmu);
  p.Ulayers.resize(p.Nlayer);
  p.Mout.resize(p.Nlayer);
  p.Lout.resize(p.Nlayer);
  p.Js.resize(p.kmax+1);
  
  //Loop over the layers to define the size of layer vectors
  for (int count=0; count < p.Nlayer; ++count){
    p.layers[count].xi.resize(p.Nmu);
    p.layers[count].Js.resize(p.kmax+1);
    p.layers[count].Nbar_mu.resize(p.kmax+1);
    p.layers[count].Kbar_mu.resize(p.kmax+1);
    p.layers[count].Nbar_integrand.resize(p.kmax+1);
    p.layers[count].Kbar_integrand.resize(p.kmax+1);
    for (int k=0; k<(p.kmax+1); ++k){
      p.layers[count].Nbar_integrand[k].resize(p.Nmu);
      p.layers[count].Kbar_integrand[k].resize(p.Nmu);
      p.layers[count].Nbar_mu[k].resize(p.Nmu);
      p.layers[count].Kbar_mu[k].resize(p.Nmu);
    }	
  }

  //Assign the radii of the outermost layers in each material 
  temp_int=0;
  for (int count=0; count<p.Nmaterial; ++count){
    p.layers[temp_int].a=material_amax[count];
    temp_int+=p.material_lay[count];
  }

  //Set the size of the precompute Plegendre vector
  p.Plegendre.resize(p.kmax+1);
  for (int i=0; i<(p.kmax+1); ++i){
    p.Plegendre[i].resize(p.Nmu);
  }	
  
  //Close file
  input.close();
  std::cout << "\t Input file read successful" <<std::endl;
  
}


////////////////////////////////////////////////////////////////////
/* Function to write the structure to a binary file that can be read in to plotting routines or in as a restart file */
void write_binary_structure(planet& p, parameters& params, std::string fileID )
{
  //Define standard file name
  std::string dir="Output/";
  
  //Open the file
  std::ofstream file;
  file.open (dir+params.run_name+"_"+fileID, std::ios::out | std::ios::binary);
  
  //Write the 'heaader' information
  params.write_binary(file);

  //Write the structure (note some things not printed)
  p.write_binary(file);
  
  file.close();
  
}

////////////////////////////////////////////////////////////////////
/* Function to read a structure from a binary file */
void read_binary_structure(planet& p, parameters& params, std::string file_name )
{
  
  //Open the file
  std::ifstream file;
  file.open (file_name, std::ios::in | std::ios::binary);
  
  //read the parameters
  params.read_binary(file);

  //read the planet
  p.read_binary(file);
  
  file.close();
  
}
