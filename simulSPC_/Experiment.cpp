//
// Created by MatthieuL on 07/11/2018.
//

#include "Experiment.h"
#include "OpticalSetup.h"
#include "MDF_gaussian.h"
#include "ThreadPool.h"
#include "PSF_planewave.h"
#include "PSF_gaussian.h"
#include "CEF_geometric.h"
#include "CEF_planewave.h"
#include <time.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <future>
#include <functional>
#include <array>

constexpr double N_AVOGADRO = 6.02E23;
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
	is_dont_comput_flag = false;

	if (is_multithreaded_)
	{

		for (int i = 0; i < nb_of_particle_; i++)
		{
			function_array_particle_[i] = std::bind(&Particle::move_and_light_matter_interraction, particle_list_[i]);
		}

		// Creates a copy of classInstace, use std::ref(classInstance) to force pass by ref
		
		for (int i = 0; i < nb_of_particle_; i++)
		{
			//async_futures_particle_tab[i] = pool_.submit(function_array_particle_[i]);
			auto future = pool_.submit(function_array_particle_[i]);
		}
		

		// Get result, should be equals to 4 (3+1)
		

		/*
		// launch threads
		for (int i = 0; i < nb_of_particle_; i++)
		{
			//async_futures_particle_tab[i] = std::async(launch::async, &Particle::move_and_light_matter_interraction, particle_list_[i]);					
			
			//auto result = pool_.enqueue(&Particle::move_and_light_matter_interraction, particle_list_[i]);

			async_futures_particle_tab[i].get();

		}
		*/
		for (int i = 0; i < nb_of_particle_; i++)
		{
			// if one of the particle has emmited a photon that was detected, 
			// we inform the exp_ object so that we don't computethe light_matter_interraction 
			// part of the particles object since only one photon can be aquired 
			// per main_cmock tick
			// NB : It should be different if there are more than one detector.
			if (async_futures_particle_tab[i].get())
			{
				Photon ph(macroClock_, 0, 0);
				photon_vector_.push_back(ph);
				// Only one photon can be measured for one tick of the experiment. Hence 
				// if one photon is detected, there is no need to compute the light_matter_interraction for the others particles.
				is_dont_comput_flag = true;
			}				
		}
	}
	else
	{ 	
		for (int i = 0; i < nb_of_particle_; i++)
		{
			
			particle_list_[i].move();
			//        particle_list_[i].inter_particle_interraction();
			if (!is_onePhotonDetected)
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
		}
	}
	//particle_list_[0].monitor_particle();
	macroClock_++;
	/*polarizationX.push_back(gsl_vector_get(particle_list_[2].inter_vector, 0));
	polarizationY.push_back(gsl_vector_get(particle_list_[2].inter_vector, 1));
	polarizationZ.push_back(gsl_vector_get(particle_list_[2].inter_vector, 2));
	*/
	
}

Experiment::~Experiment() {
	delete[] particle_list_;
	delete[] photon_array_;
	delete[] async_futures_particle_tab;
}

void Experiment::convertPhoton_vectorToList()
{
	if (photon_array_ != NULL)
		delete[] photon_array_;
	photon_array_ = new __int64[photon_vector_.size()];

	for (int i = 0; i < photon_vector_.size(); i++)
		photon_array_[i] = photon_vector_[i].compress_data();
}

//TODO list instead of vector for photon.
Experiment::Experiment() : photon_vector_(0), solvent_(this), pool_(4)
{

	// choosing and random number generator... at "random" for now e.g. gsl_rng_default gsl_rng_taus	 gsl_rng_mt19937
	rngGenerator_ = gsl_rng_alloc(gsl_rng_taus);
	// rng seed
	gsl_rng_set(rngGenerator_, time(NULL));

	// Multi-threading preparation
	num_cpus_ = std::thread::hardware_concurrency();

	pool_.init();
	//std::vector<std::thread> threads_(num_cpus_);

	read_ini_file();

	double simulation_time_s = init_parameter("simulation_time_s=", 0, iniFilevector_.size());

	time_step_ = init_parameter("time_step=", 0, iniFilevector_.size());
	nb_of_ticks_ = simulation_time_s / (time_step_*1E-9);
	space_step_ = init_parameter("space_step=", 0, iniFilevector_.size());
	is_multithreaded_ = init_parameter("is_multithreaded=", 0, iniFilevector_.size());
	double false_negative_proba = init_parameter("false_negative_proba=", 0, iniFilevector_.size());
	// We are searching the mean µ of the a poisson probability distribution P 
	// so that P(1;µ) >= false_negative_proba
	probabilityPhotonThreshold_ = 1E-6;

	while (gsl_ran_poisson_pdf(1, probabilityPhotonThreshold_) < false_negative_proba)
		probabilityPhotonThreshold_ += probabilityPhotonThreshold_ * 0.1;

	angular_step_ = init_parameter("angular_step=", 0, iniFilevector_.size());


	solvent_.box_size_radial_ = init_parameter("box_size_radial=", 0, iniFilevector_.size());
	solvent_.box_size_axial_ = init_parameter("box_size_axial=", 0, iniFilevector_.size());
	solvent_.viscosity_ = init_parameter("viscosity=", 0, iniFilevector_.size());
	solvent_.temperature_ = init_parameter("temperature=", 0, iniFilevector_.size());
	solvent_.refractive_index_ = init_parameter("refractive index=", 0, iniFilevector_.size());

	double w_microchanel_;
	double h_microchanel_;
	double flow_;

	solvent_.w_microchanel_ = init_parameter("w_microchanel_=", 0, iniFilevector_.size());
	solvent_.h_microchanel_ = init_parameter("h_microchanel_=", 0, iniFilevector_.size());
	solvent_.meanFlow_ = init_parameter("meanFlow_=", 0, iniFilevector_.size());

	opticalSetup_.laser_exc_.wl = init_parameter("laser_wl=", 0, iniFilevector_.size());
	opticalSetup_.laser_exc_.intensity = init_parameter("laser_intensity=", 0, iniFilevector_.size());
	opticalSetup_.objective_.NA = init_parameter("objective_NA=", 0, iniFilevector_.size());
	gsl_vector_set(opticalSetup_.laser_exc_.polarization_, 0, init_parameter("laser_polarization_X=", 0, iniFilevector_.size()));
	gsl_vector_set(opticalSetup_.laser_exc_.polarization_, 1, init_parameter("laser_polarization_Y=", 0, iniFilevector_.size()));
	gsl_vector_set(opticalSetup_.laser_exc_.polarization_, 2, init_parameter("laser_polarization_Z=", 0, iniFilevector_.size()));
	//TODO normalize vector polarisation
	//gsl_vector_scale(gsl_vector_set, 1 / normalization);
	//gsl_vector_scale(gsl_vector_set, 1 / normalization);

	int mdf_type = init_parameter("MDF_type=", 0, iniFilevector_.size());
	// TODO enum or #define 
	if (mdf_type == 0)
	{
		//Phenomelogic Gaussian-gaussian MDF
		mdf_ = new MDF_gaussian(this);
		// TODO abstract class flor MDF and create custom MDF class
		//mdf_->r_ = init_parameter("MDF_r_to_z_param=", 0, iniFilevector_.size());
	}
	else
	{
		// custom MDF
		mdf_ = new MDF(this);
		int psf_type = init_parameter("PSF_type=", 0, iniFilevector_.size());
		if (psf_type == 0)
		{
			mdf_->set_psf(new PSF_gaussian(this));
		}
		else if (psf_type == 1)
		{
			mdf_->set_psf(new PSF_planewave(this));
		}

		int cef_type = init_parameter("CEF_type=", 0, iniFilevector_.size());
		if (cef_type == 0)
		{
			mdf_->set_cef(new CEF_geometric(this));
		}
		else if (cef_type == 1)
		{

			mdf_->set_cef(new CEF_planewave(this));
		}
	}

	mdf_->export_mdf();

	solvent_.create_velocity_map();


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
		double concentration_nm = init_parameter("concentration_nm=", pos_particle_type[i], pos_particle_type[i + 1]);
		// The dimension of the box are in nm. 1µm^3 equals 1 femtoliter i.e. 1E-15 liter, so we multiply by 1E-15 x 1E-9 = 1E-24
		double volume_nm_3 = solvent_.box_size_axial_ * solvent_.box_size_radial_ * solvent_.box_size_radial_;
		double volume_liter = volume_nm_3 * 1E-24;

		int nb_of_particle_type = (int)(concentration_nm * 1E-9 * volume_liter * N_AVOGADRO);
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
	//TODO is nb_of_particle_ is zero
	nb_of_particle_ = particle_vector.size();

	//////////////////////////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	
	/*particle_test = new Particle;
	particle_test->y_ = 0;
	int nb_step = (int)(solvent_.box_size_radial_ / (space_step_)) + 1;
	double valeur_inter;
	Test_psf_ = new double* [nb_step];
	for (int i = 0; i < nb_step; ++i)
		Test_psf_[i] = new double[nb_step];
	for (int z = 0; z < nb_step; z++)
	{
		particle_test->z_ = z * space_step_ * 5 - (solvent_.box_size_axial_ / 2);
		for (int x = 0; x < nb_step; x++)
		{
			particle_test->x_ = x * space_step_ - (solvent_.box_size_radial_ / 2);
			delete mdf_->psf_;
			mdf_->psf_ = new PSF_planewave();
			mdf_->psf_->get_laser_intensity(this,*particle_test);
			valeur_inter = mdf_->psf_->I_x;
			Test_psf_[z][x]= valeur_inter;
		}
	}
	ofstream flux("Test_psf_zx.txt");
	for (int z = 0; z < nb_step; z++)
	{
		for (int x = 0; x < nb_step; x++)
		{
			flux << z << "; " << x << "; " << Test_psf_[z][x] << endl;
		}
	}*/

	//Multi-thread
	if (is_multithreaded_)
	{
		async_futures_particle_tab = new std::future<bool>[nb_of_particle_];
		//std::array< std::function<bool()>, 0> function_array;
		std::function<bool()> *function_array_particle_ = new std::function<bool()>[nb_of_particle_];
		for (int i = 0; i < nb_of_particle_; i++)
		{
			function_array_particle_[i] = std::bind(&Particle::move_and_light_matter_interraction, particle_list_[i]);
		}
	}

}

void Experiment::write_photon_vector()
{
	int array_size = photon_vector_.size();
	convertPhoton_vectorToList();

	ofstream export_flux("photon_list.ptn", std::ios::out | std::ios::binary);
	//TODO header of the binary file with the time_stamp_unit.	
	export_flux.write((char*)photon_array_, array_size * sizeof(__int64));
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

void Experiment::algorithm_benchmark()
{
	// choosing and random number generator... at "random" for now e.g. gsl_rng_default gsl_rng_taus	 gsl_rng_mt19937
	rngGenerator_ = gsl_rng_alloc(gsl_rng_taus2);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_mt19937);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_ranlxs0);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_ranlxd1);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_cmrg);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_gfsr4);
	//rngGenerator_ = gsl_rng_alloc(gsl_rng_ranf);		
	//
	

	// rng seed
	gsl_rng_set(rngGenerator_, time(NULL));

	int nb_of_trial = 1E8;

	std::cout << "Starting benchmark with " << nb_of_trial <<"nb of trial" << std::endl;

	double a = 5;
	double b;
	auto start = std::chrono::steady_clock::now();
	for (int i = 0; i < nb_of_trial; i++)
	{
		b = a * i;
	}
	auto end = std::chrono::steady_clock::now();
	auto diff = end - start;
	std::cout << "multiplication cost : " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;

	start = std::chrono::steady_clock::now();
	for (int i = 0; i < nb_of_trial; i++)
	{
		gsl_ran_gaussian(rngGenerator_, 1);
	}
	end = std::chrono::steady_clock::now();
	diff = end - start;
	std::cout << "gsl_ran_gaussian cost : " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;

	start = std::chrono::steady_clock::now();
	for (int i = 0; i < nb_of_trial; i++)
	{
		gsl_ran_gaussian_ziggurat(rngGenerator_, 1);
	}
	end = std::chrono::steady_clock::now();
	diff = end - start;
	std::cout << "gsl_ran_gaussian_ziggurat cost : " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;


	start = std::chrono::steady_clock::now();
	for (int i = 0; i < nb_of_trial; i++)
	{
		gsl_ran_gaussian_ratio_method(rngGenerator_, 1);
	}
	end = std::chrono::steady_clock::now();
	diff = end - start;
	std::cout << "gsl_ran_gaussian_ratio_method cost : " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;

	start = std::chrono::steady_clock::now();
	for (int i = 0; i < nb_of_trial; i++)
	{
		gsl_ran_poisson(rngGenerator_, 1);
	}
	end = std::chrono::steady_clock::now();
	diff = end - start;
	std::cout << "gsl_ran_poisson cost : " << std::chrono::duration <double, std::milli>(diff).count() << " ms" << std::endl;


}
