#pragma once
#include "PSF.h"


/**
 * \brief      State of the art PSF evaluation based on the interference of the incident  laser plane wave  around the focus point.
 *
 * \details		Give some article.
 * \sa	PSF.h
 */
class PSF_planewave :
	public PSF
{

public:

	gsl_matrix* mat_rot = gsl_matrix_calloc(3, 3);
	gsl_vector* laser_polarization_ = gsl_vector_calloc(3);
	gsl_vector* vector_polarisation = gsl_vector_calloc(3);
	double delta_theta, delta_phi;
	complex<double> E_x, E_y, E_z;



public:
	PSF_planewave();
	PSF_planewave(Experiment *exp);
	virtual double get_probability_absorption(Particle *p);	
	void vizualise_psf(string file = "psf_map.txt", string mode="xy", double dim1_nm=1500, double dim2_nm = 5000, double step_nm = 10);
	virtual void explore_psf();
	~PSF_planewave();
};

