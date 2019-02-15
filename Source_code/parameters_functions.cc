//SJL: 11/15
/*
Functions for the class structure for the parameters for the HERCULES code
*/

#include "header.h"

///////////////////////////////////////////////////////////////
//Function to write parameters to binary
void parameters::write_binary(std::ofstream& file)
{
  //Set standard string length in characters
  int string_length=200;
  char temp_char[string_length];

  //run name
  std::string temp_string=run_name;
  temp_string.insert(temp_string.end(), string_length - temp_string.size(), ' ');
  std::strcpy(temp_char, temp_string.c_str());
  file.write (temp_char, string_length*sizeof(char));

  //various param values
  file.write ((char*)&nint_max, sizeof (nint_max));
  file.write ((char*)&toll, sizeof (toll));
  file.write ((char*)&xi_nint_max, sizeof (xi_nint_max));
  file.write ((char*)&xi_toll, sizeof (xi_toll));
  file.write ((char*)&dxi, sizeof (dxi));

  //flag_start
  file.write ((char*)&flag_start, sizeof (flag_start));
  
  //start_file
  //need to convert string to standard length and char
  temp_string=start_file;
  temp_string.insert(temp_string.end(), string_length - temp_string.size(), ' ');
  std::strcpy(temp_char, temp_string.c_str());
  file.write (temp_char, string_length*sizeof(char));

  //other flags
  file.write ((char*)&flag_Mconc, sizeof (flag_Mconc));
  file.write ((char*)&flag_Lconc, sizeof (flag_Lconc));
  file.write ((char*)&flag_iter_print, sizeof (flag_iter_print));

  //Rotational profile properties
  file.write (reinterpret_cast<const char*>(&omega_param[0]), 3*sizeof(omega_param[0]));

}



///////////////////////////////////////////////////////////////
//Function to read parameters from binary
void parameters::read_binary(std::ifstream& file)
{
  //set string length
  int string_length = 200;
  char temp_char[string_length+1];
  temp_char[string_length]=0; //this stops the string from reading beyond char

  //Run name
  file.read(temp_char, string_length*sizeof(char));
  std::string temp_string=temp_char;
  temp_string.erase(remove_if(temp_string.begin(), temp_string.end(), isspace), temp_string.end()); //remove white space
  run_name=temp_string;
  
  //various parameters
  file.read((char*)&nint_max, sizeof(int));
  file.read((char*)&toll, sizeof(double));
  file.read((char*)&xi_nint_max, sizeof(int));
  file.read((char*)&xi_toll, sizeof(double));
  file.read((char*)&dxi, sizeof(double));
  
  //flags
  file.read((char*)&flag_start, sizeof(int));
  
  //start file
  file.read(temp_char, string_length*sizeof(char));
  temp_string=temp_char;
  temp_string.erase(remove_if(temp_string.begin(), temp_string.end(), isspace), temp_string.end()); //remove white space
  start_file=temp_string;
  
  //other flags
  file.read((char*)&flag_Mconc, sizeof(int));
  file.read((char*)&flag_Lconc, sizeof(int));
  file.read((char*)&flag_iter_print, sizeof(int));

  //rotational profile parameters
  omega_param.resize(3);
  file.read (reinterpret_cast<char*>(&omega_param[0]), 3*sizeof(double));
  
}
