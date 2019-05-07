//
// Created by MatthieuL on 07/11/2018.
//

#include "Experiment.h"
#include "OpticalSetup.h"
#include "MDF_gaussian.h"
#include <time.h>
#include <fstream>
using namespace std;

void Experiment::tick() {

	bool is_onePhotonDetected = false;

    for(int i=0; i<nb_of_particle_; i++)
    {
        particle_list_[i].move();
//        particle_list_[i].inter_particle_interraction();
		if(!is_onePhotonDetected)
		{ 
			// Photon probability -> branching
			if (particle_list_[i].light_matter_interraction())
			{
				Photon ph(macroClock_, 0, 0);
				photon_vector_.push_back(ph);
				// Only one photon can be measured for one tick of the experiment. Hence 
				// if one photon is detected, there is no need to compute the light_matter_interraction for the others particles.
				is_onePhotonDetected = true;
			}
		}
		macroClock_++;
    }

}

Experiment::~Experiment() {
delete[] particle_list_;
delete[] photon_array_;
}

void Experiment::convertPhoton_vectorToList()
{
	if (photon_array_ != NULL)
		delete[] photon_array_;
	photon_array_ = new __int64[photon_vector_.size()];

	for (int i = 0; i < photon_vector_.size(); i++)
		photon_array_[i] = photon_vector_[i].compress_data();
}

Experiment::Experiment() : photon_vector_(0) 
{

	// choosing and random number generator... at "random" for now e.g. gsl_rng_default gsl_rng_taus	 gsl_rng_mt19937
	rngGenerator_ = gsl_rng_alloc(gsl_rng_taus);
	// rng seed
	gsl_rng_set(rngGenerator_, time(NULL));

	read_ini_file();

	solvent_.box_size_radial_ = init_parameter("box_size_radial=", '=', ';');
	solvent_.box_size_axial_ = init_parameter("box_size_axial=", '=', ';');

	opticalSetup_.laser_exc_.wl = init_parameter("laser_wl=", '=', ';');
	opticalSetup_.laser_exc_.intensity = init_parameter("laser_intensity=", '=', ';');
	opticalSetup_.objective_.NA = init_parameter("objective_NA=", '=', ';');

    particle_list_ = new Particle[nb_of_particle_];
    for(int i=0; i<nb_of_particle_; i++)
        particle_list_[i].initParticleParam(this);


    // create PSF, CEF and MDF
	// TODO read from file the type of MDF.
	mdf_ = new MDF_gaussian(this);
    
}

void Experiment::write_photon_vector()
{
	int array_size = photon_vector_.size();
	convertPhoton_vectorToList();

	ofstream export_flux("photon_list.txt", ios::out | ios::binary);

	//TODO header of the binary file with the time_stamp_unit.
	for (int i = 0; i < array_size; i++)
	{
		export_flux << photon_array_[i] << endl;
	}

}

void Experiment::read_ini_file()
{
	ifstream data("parameters_data.txt");
	std::string line;
	while (getline(data, line))
	{
		iniFilevector_.push_back(line);
	}
}

double Experiment::init_parameter(std::string parameter, char start, char end)
{
	int found;
	int pos_start;
	int pos_end;
	std::string value;
	double numerical_value = 0;

	for (int i = 0; i < iniFilevector_.size(); i++)
	{
		if (iniFilevector_[i].find(parameter) >= 0)
		{
			pos_start = iniFilevector_[i].find(start);
			pos_end = iniFilevector_[i].find(end);
			for (int j = pos_start + 1; j < pos_end; j++)
			{
				value += iniFilevector_[i][j];
			}
			numerical_value = stod(value);
		}
	}
	return numerical_value;
}
