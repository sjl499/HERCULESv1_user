//SJL: 10/15
/*
Class structure for a planet for the HERCULES code
*/

#ifndef EOS_H
#define EOS_H

class EOS
{
public:
  //Name of input file
  std::string fname;
  //vectors for the pressure, density, temperature in EOS
  std::vector<double> p, rho, T, S;
  //length of vectors
  int Ndata;

  //type of EOS (0=standard linear, 1=standard log, 2=Hubbard) 
  int EOS_type;


  //Functions
  void read_EOSfile(std::string file_name);

  void write_binary(std::ofstream& file);
  void read_binary(std::ifstream& file);

  double calc_rho(double press);
};

#endif
