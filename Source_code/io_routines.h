//SJL: 11/15
/*
IO functions for the HERCULES code
*/

#ifndef HERCULES_IO_H
#define HERCULES_IO_H

void read_input(planet& p, parameters& param);
void write_binary_structure(planet& p, parameters& params, std::string fileID );
void read_binary_structure(planet& p, parameters& params, std::string file_name );


#endif
