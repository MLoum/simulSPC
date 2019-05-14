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
using namespace std;

int main()
{

	//TODO read a ini file
	cout << "Initialisation" << endl;
	Experiment exp;

	unsigned __int64 user_clock = 0;

	// every %
	double progression = 0;
	double progression_step = 0.01;
	unsigned __int64 notification_time = exp.nb_of_ticks_*progression_step;

	cout << "Starting simulation" << endl;
	for (int i = 0; i < exp.nb_of_ticks_; i++)
	{
		exp.tick();
		user_clock++;
		if (user_clock == notification_time)
		{
			progression += progression_step;
			cout << progression * 100 << "%" << endl;
			cout << "nb of detected photons : " << exp.photon_vector_.size() << endl;
			user_clock = 0;
		}
	}

	exp.write_photon_vector();

	cout << "Simulation end" << endl;

	return 0;


}


