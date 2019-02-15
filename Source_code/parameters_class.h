//SJL: 10/15
/*
Class structure for the parameters for the HERCULES code
*/

#ifndef PARAMETERS_H
#define PARAMETERS_H

class parameters
{
 public:
  //run name and name of start file if needed
  std::string run_name;
  //max itterations of whole model itteration
  int nint_max;
  //tollerance for whole model itteration
  double toll;

  //Max number of itterations for each xi point
  int xi_nint_max;
  //tollerance as to when to stop itteration for xi
  double xi_toll;
  //Initial stepsize for dQdxi
  double dxi;

  //flags for different parts of the code
  int flag_start;
  std::string start_file;
  int flag_Mconc, flag_Lconc, flag_iter_print;

  //Rotational profile properties
  std::vector<double> omega_param;

  //functions
  void write_binary(std::ofstream& file);
  void read_binary(std::ifstream& file);
  

};


#endif
