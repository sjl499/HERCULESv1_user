//SJL: 11/15
/*
Function that itterates to find the equilibrium structure in the HERCULES code
*/

#include "header.h"


void itterate(planet& p, parameters& params) 
{
  std::cout << "Itterating" << std::endl;

  //print the pre-itteration structure
  //EDIT FOR ODYSSEY
  std::ostringstream temp_ostring;
  temp_ostring << 0;
  write_binary_structure(p, params, temp_ostring.str());

  //set up the column headings for main itteration
  std::cout << "\t Main Itteration # \t convergence \t a \t aspect \t ellipticity \t l \t omega_rot \t Mass error \t AM error " << std::endl;
  
  //set up the counters and convergence test
  int counter=0;
  double converge=1E99;
  std::vector<double> Ulayers_old;

  //Initialise a vector of new xi 
  std::vector< std::vector<double> > xi_new;
  xi_new.resize(p.Nlayer);
  for (int i=0; i<p.Nlayer; ++i){
    xi_new[i].resize(p.Nmu);
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  //Itteration to find the final stucture
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
  //itterate until reach max itterations or convergence
  double X, X_temp, dXdxi;
  double xi_step;
  double xi_old;
  double xi_converge=1.0E99;
  int xi_counter=0;
  while ((converge>params.toll) && (counter<params.nint_max))
    {
    ++counter;
    Ulayers_old=p.Ulayers;
    
    ////////////////////////////////////////////////////////////////////////////////////
    //1) Loop to find the new equipotential surfaces
    ///////////////////////////////////////////////////////////////////////////////////
    //loop over all the layers
    for (int l=0; l<p.Nlayer; ++l){
      //loop over the point is a layer
      for (int i=0; i<p.Nmu; ++i){
	//reset criteria
	xi_counter=0;
	xi_converge=1.0E99;

	//Find dXdxi using default step size
	X=p.calc_X(p.layers[l].xi[i]*p.layers[l].a, i, l);
	X_temp=p.calc_X((p.layers[l].xi[i]+params.dxi)*p.layers[l].a, i, l);
	dXdxi=(X_temp-X)/params.dxi;
	xi_step=-X/dXdxi;

	//initialise xi value
	xi_new[l][i]=p.layers[l].xi[i];
	  
	//loop until convergence
	while ((xi_converge>params.xi_toll) && (xi_counter<params.xi_nint_max))
	  {
	    //use N-R to calculate new xi value
	    xi_old=xi_new[l][i];
	    xi_new[l][i]+=xi_step;

	    //Calculate new X and dXdxi
	    X=p.calc_X(xi_new[l][i]*p.layers[l].a, i, l);
	    X_temp=p.calc_X((xi_new[l][i]+params.dxi)*p.layers[l].a, i, l);
	    dXdxi=(X_temp-X)/params.dxi;
	    xi_step=-X/dXdxi;
	    
	    //calculate convergence criteria and increment counter
	    xi_converge=std::abs((xi_new[l][i]-xi_old)*(p.layers[l].a/p.amax));
	    ++xi_counter;

	  }
      }
    } //out of xi while loop

    ////////////////////////////////////////////////////////////////////////////////////
    //2) Set equipotential surfaces as new density surfaces and calculate planet properties
    ///////////////////////////////////////////////////////////////////////////////////
    //Assign the new xi to the xi vectors in each layer
    for (int count=0; count<p.Nlayer; ++count){
      p.layers[count].xi=xi_new[count];
      //assign the new b values
      p.layers[count].b=xi_new[count][p.Nmu-1]*p.layers[count].a;
    }

    //find the new planetary aspect
    p.aspect=xi_new[0][p.Nmu-1];

    //Recalculate the mass and volume of each of the layers and total mass
    for (int count=0; count<p.Nlayer; ++count){
      p.layers[count].calc_Mbar(p.Plegendre);
      //std::cout << "Mbar " << p.layers[count].Mbar << std::endl;
      p.layers[count].calc_M();
    }
    p.calc_Mtot();

    //Recalculate the precalc values and J's, etc.
    for (int count=0; count<p.Nlayer; ++count){
      //precompute values
      p.layers[count].precalc_Nbar(p.Plegendre);
      p.layers[count].precalc_Kbar(p.Plegendre); 
      
      //body constants 
      p.layers[count].calc_Js(p.Plegendre);
    }

    //calculate the new equitorial gravity
    for (int count=0; count<p.Nlayer; ++count){
      p.Ulayers[count]=p.calc_V(p.layers[count].a, 0);
      p.Ulayers[count]+=p.calc_Q(p.layers[count].a, 0, count);
    }


    ////////////////////////////////////////////////////////////////////////////////////
    //3) Second major section to itterate to conserve mass and AM
    ///////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////
    //Select mass conservation option
    if (params.flag_Mconc==0){
      //do nothing, your job here is done
      if (params.flag_Lconc!=0){
	std::cerr << "Lconc!=0 and Mconc=0 combination not yet programmed"  << std::endl;
	std::exit(1);
      }
    }

    ////////////////////////////////////////////////
    //More advanced Mconc options
    else if (params.flag_Mconc==1 || params.flag_Mconc==2){
      //for these flags we scale either the density or the radius to conserve mass
      
      ////////////////////////////////////////////////////////////////////////////////////
      //4) Loop until the pressure and rotation profile converges
      ///////////////////////////////////////////////////////////////////////////////////

       //set up the column headings for mass itteration
      std::cout << "\t\t Mass/AM Itteration # \t convergence \t p_core \t a \t omega_rot \t edge corot \t AM error " << std::endl;

      //Need to itterate until the change in central pressure is small
      int counter_M=0;
      double converge_M=1E99;
      double pcore_old=0.0, omega_rot_old=0.0;
      while ((converge_M>params.toll) && (counter_M<params.xi_nint_max))
	{
	  ++counter_M;

	  ////////////////////////////////////////////////////////////////////////////////////
	  //4a) Find the pressure in the planet and the corresponding density

	  //In these flags we first must calculate the pressure to find the density
	  p.press[0]=p.pmin;
	  for (int count=1; count<p.Nlayer; ++count){
	    p.press[count]=p.press[count-1]+p.real_rho[count-1]*(p.Ulayers[count]-p.Ulayers[count-1]);
	  }

	  //and find the core potential and pressure
	  p.Ucore=p.calc_V(0.0, 0);
	  p.pcore=p.press[p.Nlayer-1]+p.real_rho[p.Nlayer-1]*(p.Ucore-p.Ulayers[p.Nlayer-1]);

	  //Now use the pressure to recalculate the density 
	  //Pressure is calculated as the average of the two neighbouring zone
	  for (int count=0; count<p.Nlayer-1; ++count){
	    p.real_rho[count]=p.materials[p.flag_material[count]].calc_rho(0.5*(p.press[count]+p.press[count+1]));
	  }
	  p.real_rho[p.Nlayer-1]=p.materials[p.flag_material[p.Nlayer-1]].calc_rho(0.5*(p.press[p.Nlayer-1]+p.pcore));
	  
	  //Distribute the density over the layers
	  p.layers[0].rho=p.real_rho[0];
	  for (int count=1; count<p.Nlayer; ++count){
	    p.layers[count].rho=p.real_rho[count]-p.real_rho[count-1];
	  }
	  
	  ////////////////////////////////////////////////////////////////////////////////////
	  //4a) Mass conservation by scaling density
	  
	  //Now these conservation options split
	  if (params.flag_Mconc==1){
	    //In this case we now scale the densities to give the total mass
	    //There is a seperate scaling factor for each material layer
	    std::vector<double> rho_scale;
	    std::vector<int>::iterator pos;
	    double temp_denominator, temp_numerator;
	    int inside=0, outside=0;
	    rho_scale.resize(p.Nmaterial);
	    //loop over all materials
	    for (int count=0; count<p.Nmaterial; ++count)
	      {
		//find the zones inside of the material
		pos=std::find(p.flag_material.begin(), p.flag_material.end(), count+1);
		inside=std::distance(p.flag_material.begin(), pos);
		
		//Now find the mass of all the material in the range for that material
		temp_numerator=p.M_materials[count];
		//loop over all layers and calculate the numerator
		for (int i=0; i<p.Nlayer; ++i){
		  if (p.flag_material[i]<count){
		    if (inside==p.layers.size()) {
		      temp_numerator+=-rho_scale[p.flag_material[i]]*p.layers[i].rho*
			(p.layers[outside].Vol);
		    }
		    else {
		      temp_numerator+=-rho_scale[p.flag_material[i]]*p.layers[i].rho*
			(p.layers[outside].Vol-p.layers[inside].Vol);
		    }
		  }
		}
		
		//loop to find the denominator
		temp_denominator=0.0;
		for (int i=outside; i<inside; ++i){
		  if (inside==p.layers.size()) {
		    temp_denominator+=p.layers[i].rho*p.layers[i].Vol;
		  }
		  else {
		    temp_denominator+=p.layers[i].rho*(p.layers[i].Vol-p.layers[inside].Vol);
		  }
		}
		//calculate the scaling factor
		rho_scale[count]=temp_numerator/temp_denominator;
		
		//Move inside to outside
		outside=inside;
	      }

	    //Now use the scales to scale the density of each layer
	    for (int count=0; count<p.Nlayer; ++count){
	      p.layers[count].rho=p.layers[count].rho*rho_scale[p.flag_material[count]];
	    }
	    
	    //recalculate the real_rho
	    double temp=0.0;
	    for (int count=0; count<p.Nlayer; ++count){
	      //sum over the layers to find the 'real' density
	      temp+=p.layers[count].rho;
	      p.real_rho[count]=temp;
	    }
	    
	  }

	  ////////////////////////////////////////////////////////////////////////////////////
	  //4b) Mass conservation by scaling radius
	  
	  else {
	    //In this case we change the radii to conserve mass
	    //There is a seperate scaling factor for each material layer
	    std::vector<double> a_scale;
	    std::vector<int>::iterator pos;
	    double temp_denominator, temp_numerator;
	    int inside=0, outside=0;
	    a_scale.resize(p.Nmaterial);
	    //loop over all materials
	    for (int count=p.Nmaterial-1; count>=0; --count)
	      {
		//find the zones inside of the material
		pos=std::find(p.flag_material.begin(), p.flag_material.end(), count+1);
		inside=std::distance(p.flag_material.begin(), pos);
		pos=std::find(p.flag_material.begin(), p.flag_material.end(), count);
		outside=std::distance(p.flag_material.begin(), pos);
		
		
		//Now find the mass of all the material in the range for that material
		temp_numerator=p.M_materials[count];
		//loop over all layers and calculate the numerator
		if (inside==p.layers.size()) {
		  //Do nothing the mass is the entire numerator
		}
		else {
		  for (int i=0; i<inside; ++i){
		    temp_numerator+=a_scale[count+1]*p.layers[i].rho*p.layers[inside].Vol;
		  }  
		}
		
		//loop to find the denominator
		temp_denominator=0.0;
		for (int i=0; i<outside; ++i){
		  temp_denominator+=p.layers[i].rho*p.layers[outside].Vol;
		}
		for (int i=outside; i<inside; ++i){
		  temp_denominator+=p.layers[i].rho*p.layers[i].Vol;
		}
		//calculate the scaling factor
		a_scale[count]=temp_numerator/temp_denominator;
	      }
	    
	    //Now use the scales to scale the radius of each layer
	    for (int count=0; count<p.Nlayer; ++count){
	      p.layers[count].a=p.layers[count].a*pow(a_scale[p.flag_material[count]], (1.0/3.0));
	      p.layers[count].b=p.layers[count].xi[p.Nmu-1]*p.layers[count].a;
	    }
	    p.amax=p.layers[0].a;
	    

	    //check to see if the radii are crossing and if they are correct
	    for (int count=1; count<p.Nlayer; ++count){
	      if (p.layers[count].a >= p.layers[count-1].a)
		{
		  //std::cerr << "radii crossing " << count << " " << p.layers[count].a << " " << p.layers[count-1].a << std::endl;
		  for (int i=0; i<count; ++i){
		    p.layers[i].a=p.layers[i].a+1.005*(p.layers[count].a-p.layers[count-1].a);
		    //std::cerr << "radii crossing " << count << " " << p.layers[count].a << " " << p.layers[count-1].a << std::endl;
		    //std::exit(5);
		  }
		}
	    }

	    //Ensure there isn't a large gap between layers
	    int temp_int=p.Nlayer-1;
	    //for all layer boundaries
	    for (int count=p.Nmaterial-1; count>=1; --count){
	      //test to see if the gap over the material interface is too large compared 
	      temp_int-=p.material_lay[count];
	      double dif_bnd, dif_upper;
	      dif_bnd=p.layers[temp_int].a-p.layers[temp_int+1].a; //boundary gap
	      dif_upper=p.layers[temp_int-1].a-p.layers[temp_int].a; //normal gap for next layer
	      if (dif_bnd>1.05*dif_upper)
		{
		  //If the gap is large then move material in outer layers inwards
		  for (int i=0; i<=temp_int; ++i){
		    p.layers[i].a=p.layers[i].a-(dif_bnd-dif_upper);
		  }
		}
	    }

	    //reassign the b and amax
	    for (int count=0; count<p.Nlayer; ++count){
	      p.layers[count].b=p.layers[count].xi[p.Nmu-1]*p.layers[count].a;
	    }
	    p.amax=p.layers[0].a;
	    
	      
	  } //out of split for different options for mass conservation
	  
	  ////////////////////////////////////////////////////////////////////////////////////
	  //4c) Redefine the properties of the planet
	  
	  //now recalculate the planetary properties
	  for (int count=0; count<p.Nlayer; ++count){
	    p.layers[count].calc_Mbar(p.Plegendre);
	    //std::cout << "Mbar " << p.layers[count].Mbar << std::endl;
	    p.layers[count].calc_M();
	  }
	  //total mass
	  p.calc_Mtot();
	  
	  //Recalculate the precalc values and J's, I's etc.
	  for (int count=0; count<p.Nlayer; ++count){
	    //precompute values
	    p.layers[count].precalc_Nbar(p.Plegendre);
	    p.layers[count].precalc_Kbar(p.Plegendre); 
	    //body constants
	    p.layers[count].calc_I(p.Plegendre);    
	    p.layers[count].calc_Js(p.Plegendre);
	  }

	  //calculate the new equitorial gravity
	  for (int count=0; count<p.Nlayer; ++count){
	    p.Ulayers[count]=p.calc_V(p.layers[count].a, 0);
	    p.Ulayers[count]+=p.calc_Q(p.layers[count].a, 0, count);
	  }

	  ////////////////////////////////////////////////////////////////////////////////////////
	  //4d) Do the AM conservation

	  //index of last corotating layer defined here for printing purposes
	  int ind_orb=0;

	  //////////////////////////////
	  //if no AM conserve just do nothing
	  if (params.flag_Lconc==0){
	    //Do nothing the omega will not change 

	  }

	  //////////////////////////////
	  //else find the omega_rot that conserves AM
	  else if (params.flag_Lconc==1){
	    //We need to itterate to conserve AM
	 
	    //Sum the contributions to the moment of inertia
	    double temp_I1=0.0;
	    for (int count=0; count<p.Nlayer; ++count){
	      temp_I1+=p.layers[count].I;
	    }

	    //find the new omega_rot and assign to all the corotating shells
	    p.omega_rot=(p.Ltot_tar)/temp_I1;
	    for (int count=0; count<p.Nlayer; ++count){
	      p.layers[count].omega=p.omega_rot;
	    }

	    //calculate the new AM of the planet
	    p.calc_Ltot();
	    
	    //Finally then update the equitorial potentials
	    for (int count=0; count<p.Nlayer; ++count){
	      p.Ulayers[count]=p.calc_V(p.layers[count].a, 0);
	      p.Ulayers[count]+=p.calc_Q(p.layers[count].a, 0, count);
	    }
	    
	  }
	  else {
	    std::cerr << "Error in input. Unrecognised Lconc flag" << std::endl;
	    std::exit(1);
	  }
	 
	  
	  ////////////////////////////////////////////////////////////////////////////////////////
	  //4e) Calculate the convergence and loop

	  converge_M=std::abs((pcore_old-p.pcore)/p.pcore);

	  pcore_old=p.pcore;
	  omega_rot_old=p.omega_rot;

	  std::cout << "\t \t \t" << counter_M << "\t" << converge_M << "\t" << pcore_old << "\t" << p.amax << "\t" << omega_rot_old << "\t" << ind_orb << "\t" << (p.Ltot-p.Ltot_tar)/p.Ltot_tar << std::endl; 
	}

    }
    else {
      std::cerr << "Error in input. Unrecognised Mconc flag" << std::endl;
      std::exit(1);
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    //Calculate the convergence
    std::vector<double> temp_U;
    temp_U=VecDoub::VecSubtract(p.Ulayers,Ulayers_old);
    converge=VecDoub::VecSum(temp_U)/(p.Nlayer*p.Ulayers[0]);
    converge=std::abs(converge);
    std::cout << "\t "<< counter <<"\t" << converge << "\t" << p.amax << "\t" << p.aspect << "\t" << sqrt(1.0-pow(p.aspect,2.0)) << "\t" << sqrt(pow(p.aspect, -2.0)-1) << "\t" << p.omega_rot << "\t" << (p.Mtot-p.Mtot_tar)/p.Mtot_tar << "\t" << (p.Ltot-p.Ltot_tar)/p.Ltot_tar << std::endl;

    //If wanted print the structure after each itteration
    if (params.flag_iter_print==1){
      //calculate the mass and AM distribution
      p.calc_Mout();
      p.calc_Lout();
      p.calc_Js();
      
      //Print out file
      std::ostringstream temp_ostring;
      temp_ostring<<counter; 
      write_binary_structure(p, params, temp_ostring.str()); 
    }
    
    }
  
  //If not selected to print each itteration print the final structure
  if (params.flag_iter_print!=1){
    p.calc_Mout();
    p.calc_Lout();
    p.calc_Js();
    std::string temp_name="final";
    write_binary_structure(p, params, temp_name);
  }
}
