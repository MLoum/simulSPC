//
// Created by MatthieuL on 07/11/2018.
//

#include "Particle.h"
#include <stdio.h>
#include <gsl/gsl_randist.h>
constexpr double PI = 3.1415159265;
constexpr double kb = 1.38E-23;



Particle::Particle(Experiment &exp) {

	exp_ = &exp;
	
	/*
	x_ = gsl_rng_uniform(exp_->rngGenerator_);
	y_ = gsl_rng_uniform(exp_->rngGenerator_);
	z_ = gsl_rng_uniform(exp_->rngGenerator_);
	*/

}

void Particle::initParticleParam(Experiment *exp, int ligne_start, int ligne_end)
{
	exp_ = exp;

	double r_hydro_max, r_hydro_sigma;

	r_hydro_max = exp_->init_parameter("r_hydro=", ligne_start, ligne_end);
	r_hydro_sigma = exp_->init_parameter("r_hydro=", ligne_start, ligne_end);

	r_hydro_ = gsl_ran_gaussian(exp_->rngGenerator_, r_hydro_sigma * r_hydro_max / 100) + r_hydro_max;

	abs_cross_section_ = exp_->init_parameter("cross_section=", ligne_start, ligne_end);
	quantum_yield_ = exp_->init_parameter("quantum_yield=", ligne_start, ligne_end);
	is_point_ = (bool)(exp_->init_parameter("is_point=", ligne_start, ligne_end));
	is_isotropic_ = (bool)(exp_->init_parameter("is_isotropic=", ligne_start, ligne_end));

	Dx_ = 1E9 * (kb * exp_->solvent_.temperature_) / (6 * PI * exp_->solvent_.viscosity_ * r_hydro_ * 1E-9);
	//FIXME non isotropic particle
	Dy_ = Dx_;
	Dz_ = Dx_;
	mean_delta_x_squared_ = 2 * Dx_ * exp_->time_step_;
	mean_delta_y_squared_ = 2 * Dy_ * exp_->time_step_;
	mean_delta_z_squared_ = 2 * Dz_ * exp_->time_step_;

	// For a sphere
	Dphi1_ = 1E-9 * (kb * exp_->solvent_.temperature_) / (8 * PI * exp_->solvent_.viscosity_ * r_hydro_*r_hydro_*r_hydro_ * 1E-27);
	Dphi2_ = Dphi1_;
	Dphi3_ = Dphi1_;
	mean_delta_phi1_squared_ = 2 * Dphi1_ * exp_->time_step_;
	mean_delta_phi2_squared_ = 2 * Dphi2_ * exp_->time_step_;
	mean_delta_phi3_squared_ = 2 * Dphi3_ * exp_->time_step_;

	//TODO precomputed value of abs_cross_section_/(h*nu)*quantum_yield_
	brightness_ = abs_cross_section_ * quantum_yield_;

	x_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
	y_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_radial_;
	z_ = gsl_rng_uniform(exp_->rngGenerator_) * exp_->solvent_.box_size_axial_;

	x_ = exp_->solvent_.box_size_radial_*1.2 / 2;
	y_ = exp_->solvent_.box_size_radial_*1.2 / 2;
	z_ = exp_->solvent_.box_size_axial_*1.2 / 2;

	//TODO set identity -> one function ?
	gsl_matrix_set_identity(mat_transf_);

	// Absorption dipole
	double polarization1, polarization2, polarization3, normalization;
	polarization1 = exp_->init_parameter("polarization1=", ligne_start, ligne_end);
	polarization2 = exp_->init_parameter("polarization2=", ligne_start, ligne_end);
	polarization3 = exp_->init_parameter("polarization3=", ligne_start, ligne_end);
	normalization = sqrt(pow(polarization1, 2) + pow(polarization2, 2) + pow(polarization3, 2));
	polarization1 /= normalization;
	polarization2 /= normalization;
	polarization3 /= normalization;

	gsl_vector_set(abs_dipole_vector, 0, polarization1);
	gsl_vector_set(abs_dipole_vector, 1, polarization2);
	gsl_vector_set(abs_dipole_vector, 2, polarization3);

}

gsl_vector* Particle::get_rotated_abs_dipole_vector()
{
	gsl_blas_dgemv(CblasNoTrans, 1.0, mat_transf_, abs_dipole_vector, 0.0, abs_dipole_vector_lab);
	return abs_dipole_vector_lab;
}

void Particle::brownianMotionTest()
{
	//TODO Copy the first version
}

void Particle::move() {

	//TODO Dx * exp_->time_step_ change unit to spare a multiplciation
	/*
	x_ += Dx_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
	y_ += Dy_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
	z_ += Dz_ * gsl_ran_gaussian(exp_->rngGenerator_, 1) * exp_->time_step_;
	*/
	// Ziggurat is way faster than the Box-Muller algorithm (i.e. gsl_ran_gaussian)
	x_ += gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_x_squared_);
	y_ += gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_y_squared_);
	z_ += gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_z_squared_);

	// Periodic boundary condition
	if (y_ < 0)
	{
		y_ += exp_->solvent_.box_size_radial_;
	}
	else if (y_ >= exp_->solvent_.box_size_radial_)
	{
		y_ -= exp_->solvent_.box_size_radial_;
	}

	if (z_ < 0)
	{
		z_ += exp_->solvent_.box_size_axial_;
	}
	else if (z_ >= exp_->solvent_.box_size_axial_)
	{
		z_ -= exp_->solvent_.box_size_axial_;
	}


	// convection
	//x_ += exp_->solvent_.get_convection(this);
	if (x_ < 0)
	{
		x_ += exp_->solvent_.box_size_radial_;
	}
	else if (x_ >= exp_->solvent_.box_size_radial_)
	{
		x_ -= exp_->solvent_.box_size_radial_;
	}

	if (!is_isotropic_)
	{
		//polarization and rotation
		double phi1, phi2, phi3;
		phi1 = gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_phi1_squared_);
		phi2 = gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_phi2_squared_);
		phi3 = gsl_ran_gaussian_ziggurat(exp_->rngGenerator_, mean_delta_phi3_squared_);

		double c1 = cos(phi1), s1 = sin(phi1), c2 = cos(phi2), s2 = sin(phi2), c3 = cos(phi3), s3 = sin(phi3);
		
		gsl_matrix_set(mat_rot1, 0, 0, 1);
		gsl_matrix_set(mat_rot1, 1, 1, c1);
		gsl_matrix_set(mat_rot1, 2, 2, c1);
		gsl_matrix_set(mat_rot1, 2, 1, -s1);
		gsl_matrix_set(mat_rot1, 1, 2, s1);

		gsl_matrix_set(mat_rot2, 0, 0, c2);
		gsl_matrix_set(mat_rot2, 0, 2, s2);
		gsl_matrix_set(mat_rot2, 1, 1, 1);
		gsl_matrix_set(mat_rot2, 2, 2, c2);
		gsl_matrix_set(mat_rot2, 2, 0, -s2);

		gsl_matrix_set(mat_rot3, 0, 0, c3);
		gsl_matrix_set(mat_rot3, 0, 1, -s3);
		gsl_matrix_set(mat_rot3, 1, 0, s3);
		gsl_matrix_set(mat_rot3, 1, 1, c3);
		gsl_matrix_set(mat_rot3, 2, 2, 1);
		

		/*
		double a_pl_c = gsl_matrix_get(mat_transf_, 0, 0) + gsl_matrix_get(mat_transf_, 0, 2);
		double a_pl_c_ti_c2 = a_pl_c * c2;
		double a_pl_c_ti_s2 = a_pl_c * s2;
		double c_ti_c1 = gsl_matrix_get(mat_transf_, 0, 2)*c1;
		double c_ti_s1 = gsl_matrix_get(mat_transf_, 0, 2)*s1;
		double b_ti_c1 = gsl_matrix_get(mat_transf_, 0, 1)*c1;
		double b_ti_s1 = gsl_matrix_get(mat_transf_, 0, 1)*s1;


		
		//gsl_matrix_set(mat_transf_, 0, 0, c3*((a + c)*c2 - (c*c1 - b*s1)*s2) + (b*c1 + c*s1)*s3);
		//gsl_matrix_set(mat_transf_, 0, 1, c3*(b*c1 + c*s1) - ((a + c)*c2 - (c*c1 - b * s1)*s2)*s3);
		//gsl_matrix_set(mat_transf_, 0, 2, c2*(c*c1 - b*s1) + (a + c)*s2);
		
		gsl_matrix_set(mat_transf_, 0, 0, c3*(a_pl_c_ti_c2 - (c_ti_c1 - b_ti_s1)*s2) + (b_ti_c1 + c_ti_s1)*s3);
		gsl_matrix_set(mat_transf_, 0, 1, c3*(b_ti_c1 + c_ti_s1) - (a_pl_c_ti_c2 - (c_ti_s1 - b_ti_s1)*s2)*s3);
		gsl_matrix_set(mat_transf_, 0, 2, c2*(c_ti_c1 - b_ti_s1) + a_pl_c_ti_s2);
		
		double d_pl_f = gsl_matrix_get(mat_transf_, 1, 0) + gsl_matrix_get(mat_transf_, 1, 2);
		double d_pl_f_ti_c2 = d_pl_f * c2;
		double d_pl_f_ti_s2 = d_pl_f * s2;
		double f_ti_c1 = gsl_matrix_get(mat_transf_, 1, 2)*c1;
		double f_ti_s1 = gsl_matrix_get(mat_transf_, 1, 2)*s1;
		double e_ti_c1 = gsl_matrix_get(mat_transf_, 1, 1)*c1;
		double e_ti_s1 = gsl_matrix_get(mat_transf_, 1, 1)*s1;

		gsl_matrix_set(mat_transf_, 1, 0, c3*(d_pl_f_ti_c2 - (f_ti_c1 - e_ti_s1)*s2) + (e_ti_c1 + f_ti_s1)*s3);
		gsl_matrix_set(mat_transf_, 1, 1, c3*(e_ti_c1 + f_ti_s1) - (d_pl_f_ti_c2 - (f_ti_c1 - e_ti_s1)*s2)*s3);
		gsl_matrix_set(mat_transf_, 1, 2, c2*(f_ti_c1 - e_ti_s1) + d_pl_f_ti_s2);

		
		//gsl_matrix_set(mat_transf_, 1, 0, c3*(c2*(d + f) - (c1*f - e*s1)*s2) + (c1*e + f*s1)*s3);
		//gsl_matrix_set(mat_transf_, 1, 1, c3*(c1*e + f*s1) - (c2*(d + f) - (c1*f - e*s1)*s2)*s3);
		//gsl_matrix_set(mat_transf_, 1, 2, c2*(c1*f - e*s1) + (d + f)*s2);

		double g_pl_i = gsl_matrix_get(mat_transf_, 2, 0) + gsl_matrix_get(mat_transf_, 2, 2);
		double g_pl_i_ti_c2 = g_pl_i * c2;
		double g_pl_i_ti_s2 = g_pl_i * s2;
		double i_ti_c1 = gsl_matrix_get(mat_transf_, 2, 2)*c1;
		double i_ti_s1 = gsl_matrix_get(mat_transf_, 2, 2)*s1;
		double h_ti_c1 = gsl_matrix_get(mat_transf_, 2, 1)*c1;
		double h_ti_s1 = gsl_matrix_get(mat_transf_, 2, 1)*s1;

		
		//c2*(g + i)
		//gsl_matrix_set(mat_transf_, 2, 0, c3*(c2*(g + i) - (c1*i - h*s1)*s2) + (c1*h + i*s1)*s3);
		//gsl_matrix_set(mat_transf_, 2, 1, c3*(c1*h + i*s1) - (c2*(g + i) - (c1*i - h*s1)*s2)*s3);
		//gsl_matrix_set(mat_transf_, 2, 2, c2*(c1*i - h*s1) + (g + i)*s2);

		gsl_matrix_set(mat_transf_, 1, 0, c3*(g_pl_i_ti_c2 - (i_ti_c1 - h_ti_s1)*s2) + (h_ti_c1 + i_ti_s1)*s3);
		gsl_matrix_set(mat_transf_, 1, 1, c3*(h_ti_c1 + i_ti_s1) - (g_pl_i_ti_c2 - (i_ti_c1 - h_ti_s1)*s2)*s3);
		gsl_matrix_set(mat_transf_, 1, 2, c2*(i_ti_c1 - h_ti_s1) + g_pl_i_ti_s2);
		*/
		//TODO !


		
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat_transf_, mat_rot3, 0.0, mat_inter);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat_inter, mat_rot2, 0.0, mat_transf_);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mat_transf_, mat_rot1, 0.0, mat_inter);
		gsl_matrix_memcpy(mat_transf_, mat_inter);
		
		/*
		gsl_blas_dgemv(CblasNoTrans, 1.0, mat_tran, abs_dipole_vector, 0.0, inter_vector);
		double normalization = sqrt(pow(gsl_vector_get(inter_vector, 0), 2) + pow(gsl_vector_get(inter_vector, 1), 2) + pow(gsl_vector_get(inter_vector, 2), 2));
		gsl_vector_scale(inter_vector, 1 / normalization);
		*/
	}

}

Particle::~Particle() {
	//TODO Something is wrong here
	/*
	if (mat_transf_ != NULL)
		free(mat_transf_);
		*/


	/*
	if (mat_rot1 != NULL)
		free(mat_rot1);
	if (mat_rot2 != NULL)
		free(mat_rot2);
	if (mat_rot3 != NULL)
		free(mat_rot3);
	if (mat_inter != NULL)
		free(mat_inter);
		*/

	/*
	if (abs_dipole_vector != NULL)
		free(abs_dipole_vector);
	if (inter_vector != NULL)
		free(inter_vector);
	*/

}

bool Particle::light_matter_interraction() {

	// 
	if(exp_->is_dont_comput_flag)
	{
		return false;
	}

	// 1) Based on particle position, angle, absorption cross-section and so on, calculate 
	//via the Molecular Detection Efficiency Function (MDF), the probability to detect a photon from the molecule.	
	double photon_probability = exp_->mdf_->get_detection_efficiency(this);

	// 2) Knowing the mean probability of counting a photon, use poissonian statistics to see if it is the case.

	// invoking gsl_ran_poisson is time consuming. Consequently, we first set 
	// a probability threshold and only if the current photon probability is 
	// bigger than the threshold we compute the gsl_ran_poisson
	if (photon_probability > exp_->probabilityPhotonThreshold_)
	{
		if (gsl_ran_poisson(exp_->rngGenerator_, photon_probability) > 0)
		{
			// A photon was detected
			//TODO microtime.
			return true;
		}
	}
	else
		return false;
}

void Particle::inter_particle_interraction() {

}


Particle::Particle()
{

}

/* Debug function that should disappear at some point*/
void Particle::display_particle_position()
{
	cout << x_ << "; " << y_ << "; " << z_ << "; " << endl;
}

void Particle::monitor_particle(int mode)
{
	//TODO enum
	if (mode == 1)
	{
		Monitor_variable monitor_variable;
		monitor_variable.r_ = sqrt(x_*x_ + y_*y_);
		//monitor_variable.x_ = x_;
		//monitor_variable.y_ = y_;
		monitor_variable.z_ = z_;
		
		gsl_blas_ddot(exp_->opticalSetup_.laser_exc_.polarization_, get_rotated_abs_dipole_vector(), &monitor_variable.e_dot_mu);
		monitor_variable.e_dot_mu = monitor_variable.e_dot_mu* monitor_variable.e_dot_mu;
		
		//get_rotated_abs_dipole_vector();
		//monitor_variable.phi_x = gsl_vector_get(abs_dipole_vector_lab, 0);
		//monitor_variable.phi_y = gsl_vector_get(abs_dipole_vector_lab, 1);
		//monitor_variable.phi_z = gsl_vector_get(abs_dipole_vector_lab, 2);
		
		monitor_variable_vector.push_back(monitor_variable);
	}
	else if (mode == 2)
	{
		std::ofstream outFile("particle_rotation_abs.txt");
		//for (const auto &monitor_variable : monitor_variable_vector) outFile << monitor_variable.r_ << " " << monitor_variable.z_ << "\n";
		for (int i = 0; i < monitor_variable_vector.size(); i++)
		{
			Monitor_variable monitor_variable = monitor_variable_vector[i];
			//outFile << monitor_variable.x_ << " " << monitor_variable.y_ << " " << monitor_variable.z_ << "\n";
			//outFile << monitor_variable.r_ << " " << monitor_variable.z_ << "\n";
			//outFile << monitor_variable.phi_x << " " << monitor_variable.phi_y << " " << monitor_variable.phi_z << "\n";
			outFile << monitor_variable.r_ << " " << monitor_variable.z_ << " " << monitor_variable.e_dot_mu << "\n";
		}
	}
}
