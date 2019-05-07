// simulSPC_.cpp : Ce fichier contient la fonction 'main'. L'exécution du programme commence et se termine à cet endroit.
//


/*! \mainpage SimulSPC
 *
 * \section intro_sec Introduction
 *
 * SimulSPC is a c++ comand line programm intended to generate a list of photon time arrival via the numerical simulation of a single particle spectroscopy experiment.
 *
 * \section Dependencies Dependencies
 *
 * Gnu Scientif Library

 * \subsection step1 Step 1
 *
 * \section Source_code Source code
 *
 * https://github.com/MLoum/simulSPC
 */

#include "pch.h"
#include <iostream>
#include "gsl/include/gsl/gsl_integration.h"
#include "Experiment.h"


int main()
{

//TODO read a ini file
Experiment exp;

for(int i=0 ; i < 20; i++)
	exp.tick();
	
	exp.write_photon_vector();

	return 0;
	

}


