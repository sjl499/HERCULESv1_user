//SJL: 10/15
/*
The header file for the HERCULES code
*/

//USe a header gaurd to avoid recompiling
#ifndef HERCULES_H
#define HERCULES_H

//Added for compilation on oddysey
#include <cstdlib>
#include <sstream>

//include io
#include <iostream>
#include <fstream>
//cmath include maths functions
#include <cmath>
//Vector
#include <vector>
//Allows some numeric functions on vectors
#include <numeric>
//Strings
#include <string>



//Legendre polynomials
#include <boost/math/special_functions/legendre.hpp>

//Header for my self defined vector namespaces
#include "vector_operations.h"

//Header for general constants
#include "constants.h"

//Header for other general functions
#include "general_functions.h"

//include the planet and layer classes
#include "parameters_class.h"
#include "concentric_layer_class.h"
#include "EOS_class.h"
#include "planet_class.h"


//include the io functions
#include "io_routines.h"

//For developmetn only
#include <ctime>

//functions defined in other files
void initialisation(planet& p, parameters& params);
void itterate(planet& p, parameters& params);




//end protected segment
#endif
