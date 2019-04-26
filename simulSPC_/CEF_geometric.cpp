#include "CEF_geometric.h"



double CEF_geometric::getCollectionProbability(Particle & p)
{
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
	return 0;
}

CEF_geometric::CEF_geometric()
{
}


CEF_geometric::~CEF_geometric()
{
}
