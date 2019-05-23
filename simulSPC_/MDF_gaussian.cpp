#include "MDF_gaussian.h"
#define PLANCK_CST  6.626070040E-34
#define SPEED_OF_LIGHT 299792458


double MDF_gaussian::get_detection_efficiency(Particle *p)
{
	double r_center_float = sqrt((p->x_ * p->x_) + (p->y_ * p->y_));
	double z_center_float = p->z_;

	int r_center_int = (int)(r_center_float / exp_->space_step_);
	int z_center_int = (int)(z_center_float / exp_->space_step_);
	double photon_probability(0);

	if (p->is_point_)
	{
		if (p->is_isotropic_)
		{
			return mdf_[r_center_int][z_center_int] * p->brightness_;
		}
		else
		{
			//FIXME
			//return mdf_[r_center_int][z_center_int] * p->brightness_;

			// absorption dipole
			//TODO complex polarization
			double dot_product;
			gsl_blas_ddot(exp_->opticalSetup_.laser_exc_.polarization_, p->get_rotated_abs_dipole_vector(), &dot_product);
			return mdf_[r_center_int][z_center_int] * p->brightness_ *dot_product*dot_product;

			//gsl_blas_zdotu
			//gsl_blas_dznrm2
		}
	}


	int x_centre, y_centre, z_centre;
	int x_min, y_min, z_min;
	int x_max, y_max, z_max;
	int nb_step_x, nb_step_y, nb_step_z;
	double R;
	int r;
	double MDF_totale(0), nb_case(0);

	if (p->r_hydro_ < exp_->space_step_ / 2)
	{
		photon_probability = mdf_[r_center_int][z_center_int] * p->brightness_;
	}
	else
	{
		x_centre = (int)(p->x_ / exp_->space_step_);
		y_centre = (int)(p->y_ / exp_->space_step_);
		z_centre = (int)(p->z_ / exp_->space_step_);
		x_min = (int)((p->x_ - p->r_hydro_) / exp_->space_step_);
		y_min = (int)((p->y_ - p->r_hydro_) / exp_->space_step_);
		z_min = (int)((p->z_ - p->r_hydro_) / exp_->space_step_);
		x_max = (int)((p->x_ + p->r_hydro_) / exp_->space_step_);
		y_max = (int)((p->y_ + p->r_hydro_) / exp_->space_step_);
		z_max = (int)((p->z_ + p->r_hydro_) / exp_->space_step_);
		nb_step_x = (int)(exp_->solvent_.box_size_radial_ / exp_->space_step_);
		nb_step_y = (int)(exp_->solvent_.box_size_radial_ / exp_->space_step_);
		nb_step_z = (int)(exp_->solvent_.box_size_axial_ / exp_->space_step_);

		for (int x = x_min; x <= x_max; x++)
		{
			for (int y = y_min; y <= y_max; y++)
			{
				for (int z = z_min; z <= z_max; z++)
				{
					R = sqrt(pow((x_centre - x), 2) + pow((y_centre - y), 2) + pow((z_centre - z), 2));
					{
						if (R <= p->r_hydro_)
						{
							if (x >= 0 && y >= 0 && z >= 0 && x <= nb_step_x && y <= nb_step_y && z <= nb_step_z)
							{
								r = (int)sqrt(pow(x, 2) + pow(y, 2));
								MDF_totale += mdf_[r][z];
							}
							nb_case += 1;
						}
					}
				}
			}
		}
		photon_probability = MDF_totale / nb_case * p->abs_cross_section_ * p->quantum_yield_;
	}

	return photon_probability;

}

MDF_gaussian::MDF_gaussian(Experiment *experiment, double r)
{
	exp_ = experiment;
	r_ = r;

	//FIX ME, I have to add a +1, because there is an overflow of the mdf when x=axial_size and y=axial_size
	int nb_step_r = (int)(sqrt(2) * exp_->solvent_.box_size_radial_ / exp_->space_step_) + 1;
	int nb_step_z = (int)(exp_->solvent_.box_size_axial_ / exp_->space_step_);

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
	double norm_coeff = exp_->opticalSetup_.laser_exc_.intensity*1E-3 / Pi_waist_square_nm / (PLANCK_CST*SPEED_OF_LIGHT / (exp_->opticalSetup_.laser_exc_.wl*1E-9));


	// Gaussian-Gaussian


	for (int z = 0; z < nb_step_z; z++)
	{
		double z_c = (z - nb_step_z / 2) * exp_->space_step_;
		double z_c_squared = z_c * z_c;

		for (int r = 0; r < nb_step_r; r++)
		{
			double r_c = (r - nb_step_r / 2) * exp_->space_step_;
			mdf_[r][z] = norm_coeff * exp(-2 * r_c*r_c / radial_waist_square_nm - 2 * z_c_squared / axial_waist_square_nm);
		}
	}
}

void MDF_gaussian::export_mdf()
{
	//FIXME filename etc...
	std::ofstream outFile("mdf_gaussian.txt");
	
	int nb_step_r = (int)(sqrt(2) * exp_->solvent_.box_size_radial_ / exp_->space_step_) + 1;
	int nb_step_z = (int)(exp_->solvent_.box_size_axial_ / exp_->space_step_);

	for (int i = 0; i < nb_step_r; i++)
	{
		for (int j = 0; j < nb_step_z; j++)
		{
			outFile << mdf_[i][j] << " ";
		}
		outFile << "\n";
	}
}


/*
MDF_gaussian::MDF_gaussian()
{
}
*/

MDF_gaussian::~MDF_gaussian()
{
	int nb_step_r = (int)(sqrt(2) * exp_->solvent_.box_size_radial_ / exp_->space_step_) + 1;
	if (mdf_ != NULL)
	{
		for (int i = 0; i < nb_step_r; ++i) {
			delete[] mdf_[i];
		}
		delete[] mdf_;
	}

}
