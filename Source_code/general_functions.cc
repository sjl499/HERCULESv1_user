//SJL 11/15
/*
File containing various functions used in the HERCULES code
*/
#include <vector>
//include io
#include <iostream>

//Added for compilation on oddysey
#include <cstdlib>

////////////////////////////////////////////////
//Function to calculate the integral of a vector with evenly spaced points between two given points
double integrate(std::vector<double>& integrand, double h, int ind_start, int ind_end)
{
  double result;
  //test to see what order integration to use
  if (ind_end-ind_start>=5){
    //Using extended simpson's (eqn 4.1.14 in numerical methods)
    result=(3.0/8.0)*(integrand[ind_start]+integrand[ind_end]);
    result+=(7.0/6.0)*(integrand[ind_start+1]+integrand[ind_end-1]);
    result+=(23.0/24.0)*(integrand[ind_start+2]+integrand[ind_end-2]);
    for (int i=ind_start+3; i<ind_end-2; ++i){
      result+=integrand[i];
    }
  }
  else if (ind_end-ind_start==1){
result=0.5*(integrand[ind_start]+integrand[ind_end]);
}
  else if (ind_end==ind_start){
result=0.0;
  }
  else if(ind_end-ind_start==4)
    {
//Using Bode's rule (eqn 4.1.6 in numerical methods)
result=(14.0)*(integrand[ind_start]+integrand[ind_end]);
result+=(64.0)*(integrand[ind_start+1]+integrand[ind_end-1]);
result+=(24.0)*integrand[ind_start+2];
result=result/45.0;
}
  else if (ind_end-ind_start==3){
//Using simpson's 3/8th rule (eqn 4.1.5 in numerical methods)
result=(3.0)*(integrand[ind_start]+integrand[ind_end]);
result+=(9.0)*(integrand[ind_start+1]+integrand[ind_end-1]);
result=result/8.0;
}
  else if (ind_end-ind_start==2){
//Using simpson's rule (eqn 4.1.4 in numerical methods)
result=integrand[ind_start]+integrand[ind_end];
result+=(4.0)*integrand[ind_start+1];
result=result/3.0;
}
  else
    {
      // Print an error and exit
      std::cerr << "Error in integrate. Non integer ind difference" << std::endl;
      std::exit(2);
}
return (result*h);
}
