#include "PSF_planewave.h"
constexpr double PI = 3.1415159265;



PSF_planewave::PSF_planewave()
{
}

PSF_planewave::PSF_planewave(Experiment * exp)
{
	exp_ = exp;
}

double PSF_planewave::get_probability_absorption(Particle *p)
{
	delta_phi = exp_->angular_step_;
	delta_theta = exp_->angular_step_;
	double phi_max = 2 * PI;
	double theta_max = asin(exp_->opticalSetup_.objective_.NA / exp_->solvent_.refractive_index_);

	//TODO complex polarization.
	gsl_vector_memcpy(laser_polarization_, exp_->opticalSetup_.laser_exc_.polarization_);

	double c_phi, s_phi, c_theta, s_theta, c_theta_sqrt;
	double phi, theta;
	double a, b, c, d, e, f;
	double k = 2 * PI / exp_->opticalSetup_.laser_exc_.wl;
	complex<double> i(0, 1);
	complex<double> shared_term;
	E_x = 0;
	E_y = 0;
	E_z = 0;

	double x, y, z;
	x = p->x_;
	y = p->y_;
	z = p->z_;

	phi = 0;
	while (phi < phi_max)
	{
		c_phi = cos(phi);
		s_phi = sin(phi);

		theta = 0;
		while (theta < theta_max)
		{
			c_theta = cos(theta);
			s_theta = sin(theta);
			c_theta_sqrt = sqrt(c_theta);


			a = 1 + (c_theta - 1) * c_phi * c_phi;
			b = (c_theta - 1) * c_phi * s_phi;
			c = -s_theta * c_phi;
			d = 1 + (c_theta - 1) * s_phi * s_phi;
			e = -s_theta * s_phi;
			f = c_theta;

			gsl_matrix_set(mat_rot, 0, 0, a);
			gsl_matrix_set(mat_rot, 0, 1, b);
			gsl_matrix_set(mat_rot, 0, 2, c);
			gsl_matrix_set(mat_rot, 1, 0, b);
			gsl_matrix_set(mat_rot, 1, 1, d);
			gsl_matrix_set(mat_rot, 1, 2, e);
			gsl_matrix_set(mat_rot, 2, 0, -c);
			gsl_matrix_set(mat_rot, 2, 1, -e);
			gsl_matrix_set(mat_rot, 2, 2, f);

			gsl_blas_dgemv(CblasNoTrans, 1.0, mat_rot, laser_polarization_, 0.0, vector_polarisation);

			shared_term = s_theta * c_theta_sqrt * exp(i * k * (z * c_theta + s_theta * (x * c_phi + y * s_phi))) * delta_phi * delta_theta;

			E_x += gsl_vector_get(vector_polarisation, 0) * shared_term;
			E_y += gsl_vector_get(vector_polarisation, 1) * shared_term;
			E_z += gsl_vector_get(vector_polarisation, 2) * shared_term;

			theta += delta_theta;
		}
		phi += delta_phi;
	}
	/*
	I_x = norm(E_x);
	I_y = norm(E_y);
	I_z = norm(E_z);
	*/

	return 0.0;
}


void PSF_planewave::vizualise_psf(string file, string mode, double dim1_nm, double dim2_nm, double step_nm)
{
	int nb_step_1 = dim1_nm / step_nm;
	int nb_step_2 = dim2_nm / step_nm;
	
	double** intensity_map = new double*[nb_step_1];
	for (int i = 0; i < nb_step_1; i++)
		intensity_map[i] = new double[nb_step_2];
	
	Particle p;
	for (int i = 0; i < nb_step_1; i++)
	{
		p.x_ = i * step_nm;
		for (int j = 0; j < nb_step_2; j++)
		{
			p.y_ = j * step_nm;
			//TODO
			//intensity_map[i][j] = get_laser_intensity(exp_, p);
		}
	}

	if (mode == "xx")
	{

	}
}

void PSF_planewave::explore_psf()
{
}

PSF_planewave::~PSF_planewave()
{
}

