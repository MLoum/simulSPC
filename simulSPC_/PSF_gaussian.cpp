#include "PSF_gaussian.h"
#include <gsl/gsl_math.h>



E_field PSF_gaussian::get_E_field(Particle & p)
{
	return E_field();
}

PSF_gaussian::PSF_gaussian()
{
}

PSF_gaussian::PSF_gaussian(OpticalSetup OpticalSetup_)
{
	int nb_step_r = (int) (exp_->solvent_.box_size_axial_ / exp_->space_step_);
	int nb_step_z = (int)(exp_->solvent_.box_size_radial_ / exp_->space_step_);


	//allocation (https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new)
	psf_ = new double*[nb_step_r];
	for (int i = 0; i < nb_step_r; ++i)
		psf_[i] = new double[nb_step_z];

	double airy_diameter_ = 1.22 * exp_->opticalSetup_.laser_exc_.wl / exp_->opticalSetup_.objective_.NA;
	double w_xy_ = (airy_diameter_ / 2)*0.66;

	//FIXME ! Should be read from an ini paramater file.
	double w_z_ = 10 * w_xy_;


	// Une Gaussienne	lorentzienne centrée en 	(x = TailleTableau/2,y =TailleTableau/2, z=TailleTableau/2), soit r=TailleTableau/2 z=TailleTableau/2
	// il s'agit ac la symetrie de révolution de I(r,z) -> coordonnée cylindrique

	// Le profil est normalisé de tel sorte qu'en r = nb_step_r/2 et z = nb_step_z/2, le profil vale 1/(Pi w_0²) soit la surface du waist

	double radial_waist_square_nm = w_xy_ * w_xy_;
	double axial_waist_square_nm = w_z_ * w_z_;
	double Pi_waist_square_nm = w_xy_ * w_xy_* M_PI;

	// Gaussian-Gaussian
	for (int z = 0; z < nb_step_z; z++)
	{
		double z_c = (z - nb_step_z/2) * exp_->space_step_;
		//double inter = pow(z_c * exp_->opticalSetup_.laser_exc_.wl/Pi_waist_square_nm, 2);

		for (int r = 0; r < nb_step_r; r++)
		{
			double r_c = (r - nb_step_r/2) * exp_->space_step_;
			psf_[r][z] = exp(-2 * r_c*r_c / radial_waist_square_nm)*exp(-2 * z_c*z_c / axial_waist_square_nm);
		}
	}
	//TODO gaussian beam with rayleigh length for the axial direction
}



PSF_gaussian::~PSF_gaussian()
{
	/*
	int nb_step_r = exp_->solvent_.box_size_axial_ / (OpticalSetup_->laser_exc_.wl / 10);
	if (psf_ != NULL)
	{
		for (int i = 0; i < nb_step_r; ++i) {
			delete[] psf_[i];
		}
	*/
}
