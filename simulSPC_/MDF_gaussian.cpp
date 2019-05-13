#include "MDF_gaussian.h"
#define PLANCK_CST  6.626070040E-34
#define SPEED_OF_LIGHT 299792458


double MDF_gaussian::get_detection_efficiency(Particle *p)
{
	double r_absolue = sqrt((p->x_ * p->x_) + (p->y_ * p->y_));
	double z_absolue = p->z_;

	int r = (int)r_absolue / exp_->space_step_;
	int z = (int)z_absolue / exp_->space_step_;

	double photon_probability = mdf_[r][z] * p->abs_cross_section_ * p->quantum_yield_;

	return photon_probability;

}

MDF_gaussian::MDF_gaussian(Experiment *experiment, double r)
{
	exp_ = experiment;
	r_ = r;

	int nb_step_r = (int)(sqrt(2) * exp_->solvent_.box_size_axial_ / exp_->space_step_);
	int nb_step_z = (int)(exp_->solvent_.box_size_radial_ / exp_->space_step_);

	//allocation (https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new)
	mdf_ = new double*[nb_step_r];
	for (int i = 0; i < nb_step_r; ++i)
		mdf_[i] = new double[nb_step_z];

	double airy_diameter_ = 1.22 * exp_->opticalSetup_.laser_exc_.wl / exp_->opticalSetup_.objective_.NA;
	double w_xy_ = (airy_diameter_ / 2)*0.66;

	double w_z_ = r_ * w_xy_;

	// Une Gaussienn-gaussienne centrée en 	(x = TailleTableau/2,y =TailleTableau/2, z=TailleTableau/2), soit r=TailleTableau/2 z=TailleTableau/2
	// il s'agit ac la symetrie de révolution de I(r,z) -> coordonnée cylindrique

	// Le profil est normalisé de tel sorte qu'en r = nb_step_r/2 et z = nb_step_z/2, le profil vale 1/(Pi w_0²) soit la surface du waist

	double radial_waist_square_nm = w_xy_ * w_xy_;
	double axial_waist_square_nm = w_z_ * w_z_;
	double Pi_waist_square_nm = w_xy_ * w_xy_* M_PI;
	
	// 1E-3 because intensity is in mW
	double norm_coeff = exp_->opticalSetup_.laser_exc_.intensity*1E-3 / Pi_waist_square_nm/(PLANCK_CST*SPEED_OF_LIGHT/(exp_->opticalSetup_.laser_exc_.wl*1E-9));
	
	// Gaussian-Gaussian
	for (int z = 0; z < nb_step_z; z++)
	{
		double z_c = (z - nb_step_z / 2) * exp_->space_step_;		
		double z_c_squared = z_c * z_c;

		for (int r = 0; r < nb_step_r; r++)
		{
			double r_c = (r - nb_step_r / 2) * exp_->space_step_;
			mdf_[r][z] = norm_coeff *  exp(-2 * r_c*r_c / radial_waist_square_nm - 2*z_c_squared/axial_waist_square_nm);
		}
	}

}

void MDF_gaussian::export_mdf()
{
}


/*
MDF_gaussian::MDF_gaussian()
{
}
*/

MDF_gaussian::~MDF_gaussian()
{
/*
int nb_step_r = exp_->solvent_.box_size_axial_ / exp_->space_step_;
if (mdf_ != NULL)
{
	for (int i = 0; i < nb_step_r; ++i) {
		delete[] mdf_[i];
	}
*/
}
