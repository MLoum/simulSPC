//
// Created by MatthieuL on 07/11/2018.
//

#include "MDF.h"
#include <math.h>
# define M_PIl          3.141592653589793238462643383279502884L /* pi */


MDF::MDF(Experiment &exp)
{
    // Allocate memory
    *exp_ = exp;
    nb_hpoint_ = int(exp_->solvent_.box_size_radial_ / exp_->space_step_);
    nb_vpoint_ = int(exp_->solvent_.box_size_axial_ / exp_->space_step_);
    mdf_ = new double*[nb_hpoint_];
    for(int i = 0; i < nb_hpoint_; ++i)
        mdf_[i] = new double[nb_vpoint_];
}

MDF::~MDF() {
    for(int i = 0; i < nb_hpoint_; ++i)
        delete(mdf_[i]);
    delete(mdf_);

}



void MDF::create_PSF() {
    airy_diameter_ = 1.22 * exp_->opticalSetup_.laser_exc_.wl/exp_->opticalSetup_.objective_.NA;
    w_xy_ = (airy_diameter_/2)*0.66;
    //FIXME !
    w_z_ = 10*w_xy_;


    int nb_hpoint_over_2 = 	nb_hpoint_ / 2;
    int nb_vpoint_over_2 = 	nb_vpoint_ / 2;

    // Une Gaussienne	lorentzienne	 centrée en 	(x = TailleTableau/2,y =TailleTableau/2, z=TailleTableau/2), soit r=TailleTableau/2 z=TailleTableau/2
    // il s'agit ac la symetrie de révolution de I(r,z) 	coordonnée cylindrique

    // Le profil est normalisé de tel sorte qu'en r = TailleTabRadialeSur2 et z = TailleTabSur2, le profil vale 1/(Pi w_0²) soir la surface du waist

    double radial_waist_square_nm = w_xy_*w_xy_;
    double axial_waist_square_nm  = w_z_*w_z_;
    double Pi_waist_square_nm = w_xy_*w_xy_* M_PIl;

// Gaussian-Gaussian
    for(int z = 0 ; z < nb_vpoint_; z++)
    {
        double z_c = (z - nb_vpoint_over_2) * exp_->space_step_;
        //double inter = pow(z_c * exp_->opticalSetup_.laser_exc_.wl/Pi_waist_square_nm, 2);

        for(int r = 0 ; r < nb_hpoint_; r++)
        {
            double r_c = (r - nb_hpoint_over_2) * exp_->space_step_;
            mdf_[r][z] =  exp(-2*r_c*r_c/radial_waist_square_nm)*exp(-2*z_c*z_c/axial_waist_square_nm);
        }
    }
//TODO gaussian beam with rayleigh length for the axial direction
}

void MDF::create_CEF() {
/*
    LoEmission = 510;
a = 		OpticalSetup->TailleTrouFiltrageEnMicron*1000 /  OpticalSetup->GrandissementObjectif;
aCarre= a * a ;
w0 = LoEmission/(PI * tan(asin(OpticalSetup->OuvertureNumeriqueObj/OpticalSetup->IndiceMilieu))); //en nm
woCarre =  w0 * w0 ;
TailleTabSur2 = Param->TailleTableau/2;
    for(int z = 0 ; z < nb_vpoint_; z++)
    {
        double z_c = (z - nb_vpoint_over_2) * exp_->space_step_;
        //double inter = pow(z_c * exp_->opticalSetup_.laser_exc_.wl/Pi_waist_square_nm, 2);


	z_en_nm = (z - TailleTabSur2) * Param->PasEntreDeuxCaseEnNm;
	R_de_z = w0 * sqrt(1 + pow(z_en_nm *  LoEmission/(PI*woCarre),2));
	R_de_z_carre = R_de_z*R_de_z;
        for(int r = 0 ; r < nb_hpoint_; r++)
		{
		r_en_nm = abs_db((r - TailleTabRadialeSur2) * Param->PasEntreDeuxCaseEnNm);
			if(r_en_nm >= R_de_z + a)
			mdf[r][z] = 0;


			else if(r_en_nm > abs_db(R_de_z - a) && r_en_nm < R_de_z + a)
			{
			Theta1 = acos((-R_de_z_carre + r_en_nm*r_en_nm + aCarre)/(2*a*r_en_nm));
			Theta2 = acos((-aCarre + r_en_nm*r_en_nm + R_de_z_carre)/(2*R_de_z*r_en_nm));
			Delta =  sqrt((a + r_en_nm + R_de_z) * (-a + r_en_nm + R_de_z) * (a - r_en_nm + R_de_z) * (a + r_en_nm - R_de_z)/2);
			mdf[r][z] *= pow(max_db(a,w0),2)  * ( Theta1*aCarre + Theta2*R_de_z_carre -  Delta) / (PI * aCarre * R_de_z_carre);
			}

			else if( r_en_nm < abs_db(R_de_z - a))
			mdf[r][z] *= pow(max_db(a,w0),2) / pow(max_db(a,R_de_z),2);

		}
	}
*/
}