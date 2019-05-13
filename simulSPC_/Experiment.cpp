//
// Created by MatthieuL on 07/11/2018.
//

#include "Experiment.h"
#include "OpticalSetup.h"
#include "MDF_gaussian.h"
#include <time.h>
#include <fstream>
#include <thread>
using namespace std;

void Experiment::tick() {
	/*for (unsigned i = 0; i < num_cpus; ++i) {
    threads[i] = std::thread([&iomutex, i] {
      {
        // Use a lexical scope and lock_guard to safely lock the mutex only for
        // the duration of std::cout usage.
        std::lock_guard<std::mutex> iolock(iomutex);
        std::cout << "Thread #" << i << " is running\n";
      }

      // Simulate important work done by the tread by sleeping for a bit...
      std::this_thread::sleep_for(std::chrono::milliseconds(200));

    });
  }

  for (auto& t : threads) {
    t.join();
  }
  
 thread : 
 https://thispointer.com/c11-start-thread-by-member-function-with-arguments/
  https://stackoverflow.com/questions/29754438/c11-multithreading-with-class-member-function
  https://www.acodersjourney.com/top-20-cplusplus-multithreading-mistakes/
  https://hackernoon.com/learn-c-multi-threading-in-5-minutes-8b881c92941f


  async :
  https://eli.thegreenplace.net/2016/the-promises-and-challenges-of-stdasync-task-based-parallelism-in-c11/
  https://stackoverflow.com/questions/11758414/class-and-stdasync-on-class-member-in-c
  /*
  Agree. Here is a real-world example explaining your mistake. Imagine a pregnant woman : she needs 9 monthes to "make" a baby. Now if you add 7 more women, the baby won't arrive faster because pregnancy is not meant to be shared amongst several women... Same for your program, 8 core is not faster than 1 core in your case. You must have a way to split the work into parts, assign them to workers and then reassemble the results.
  */

	/*As pointed out in other comments, you cannot create a thread "on a specific core", as C++ has no knowledge of such architectural details. Moreover, in the majority of cases, the operating system will be able to manage the distribution of threads among cores/processors well enough

  std::thread th(&Task::move, particle_list_[i]);
  std::thread th(&Task::light_matter_interraction, particle_list_[i]);
  */

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

	// Multi-threading preparation
	num_cpus_ = std::thread::hardware_concurrency();
	std::vector<std::thread> threads_(num_cpus_);

	read_ini_file();

	time_step_ = init_parameter("time_step=", 0, iniFilevector_.size());

	space_step_ = init_parameter("space_step=", 0, iniFilevector_.size());

	solvent_.box_size_radial_ = init_parameter("box_size_radial=", 0, iniFilevector_.size());
	solvent_.box_size_axial_ = init_parameter("box_size_axial=", 0, iniFilevector_.size());
	solvent_.viscosity_ = init_parameter("viscosity=", 0, iniFilevector_.size());
	solvent_.temperature_ = init_parameter("temperature=", 0, iniFilevector_.size());

	opticalSetup_.laser_exc_.wl = init_parameter("laser_wl=", 0, iniFilevector_.size());
	opticalSetup_.laser_exc_.intensity = init_parameter("laser_intensity=", 0, iniFilevector_.size());
	opticalSetup_.objective_.NA = init_parameter("objective_NA=", 0, iniFilevector_.size());

    // create PSF, CEF and MDF
	// TODO read from file the type of MDF.
	mdf_ = new MDF_gaussian(this);

	vector<int> pos_particle_type(0);
	int found;

	for (int i = 0; i < iniFilevector_.size(); i++)
	{
		found = iniFilevector_[i].find("<Particle");
		if (found > 0)
		{
			pos_particle_type.push_back(i);
		}
	}

	pos_particle_type.push_back(iniFilevector_.size());

	vector<Particle> particle_vector(0);
	Particle* particle_inter;
	int pos_vector(0);

	for (int i = 0; i < (pos_particle_type.size() - 1); i++)
	{
		int nb_of_particle_type = (int)init_parameter("nb_of_particle=", pos_particle_type[i], pos_particle_type[i + 1]);
		for (int j = 0; j < nb_of_particle_type; j++)
		{
			particle_inter = new Particle;
			particle_vector.push_back(*particle_inter);
			particle_vector[pos_vector].initParticleParam(this, pos_particle_type[i], pos_particle_type[i + 1]);
			pos_vector += 1;
		}
	}

	particle_list_ = new Particle[particle_vector.size()];


	for (int i = 0; i < particle_vector.size(); i++)
	{
		particle_list_[i] = particle_vector[i];
	}

	nb_of_particle_ = particle_vector.size();
    
}

void Experiment::write_photon_vector()
{
	int array_size = photon_vector_.size();
	convertPhoton_vectorToList();

	ofstream export_flux("photon_list.txt", std::ios::out | std::ios::binary);
	//TODO header of the binary file with the time_stamp_unit.	
	export_flux.write((char*) photon_array_, array_size*sizeof(__int64));
	/*
	for (int i = 0; i < array_size; i++)			
	{
		{
			//export_flux.write((char*)& photon_array_[i], sizeof(__int64));
			export_flux << (char)photon_array_[i];
		}
	}
	*/
	export_flux.close();
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

double Experiment::init_parameter(string parameter, int ligne_start, int ligne_end)
{
	int found;
	int pos_start;
	int pos_end;
	string string_average, string_sigma;
	double value;
	double average(0), sigma(0);

	for (int i = ligne_start; i < ligne_end; i++)
	{
		found = iniFilevector_[i].find(parameter);
		if (found > 0)
		{
			pos_start = iniFilevector_[i].find('{');
			pos_end = iniFilevector_[i].find(';');
			for (int j = pos_start + 1; j < pos_end; j++)
			{
				string_average += iniFilevector_[i][j];
			}
			average = stod(string_average);

			pos_start = iniFilevector_[i].find(';');
			pos_end = iniFilevector_[i].find('}');
			for (int k = pos_start + 1; k < pos_end; k++)
			{
				string_sigma += iniFilevector_[i][k];
			}
			sigma = stod(string_sigma);
		}
	}
	if (sigma = 0)
	{
		value = average;
	}
	else
	{
		value = gsl_ran_gaussian(rngGenerator_, sigma * average / 100) + average;
	}

	return value;
}
