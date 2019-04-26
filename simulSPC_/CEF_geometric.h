#pragma once
/**
 * \file          CEF_geometric.h
 * \author    Moltech-Anjou UMR 6200 Samson
 * \version   1.0
 * \date       2019 04 26
 * \brief      Implementation of Quian and Elson geometric Collection Efficiency Function
 *
 * \details    See the article Qian, H., & Elson, E. L. (1991). Analysis of confocal laser-microscope optics for 3-D fluorescence correlation spectroscopy. Applied Optics, 30(10), 1185–1195. https://doi.org/https://doi.org/10.1364/AO.30.001185
 *
 */

#include "CEF.h"

 /**
  * \date       2019 04 26
  * \brief      Implementation of Quian and Elson geometric Collection Efficiency Function
  *
  * \details    This CEF doesn't take into account the vectorial nature of the electrical field. See the article Qian, H., & Elson, E. L. (1991). Analysis of confocal laser-microscope optics for 3-D fluorescence correlation spectroscopy. Applied Optics, 30(10), 1185–1195. https://doi.org/https://doi.org/10.1364/AO.30.001185
  *
  */
class CEF_geometric :
	public CEF
{
private:
	double **mdf_;
	int nb_hpoint_, nb_vpoint_;
	double airy_diameter_; // Airy disk diameter -> 1.22 * lamdba/NA
	double w_xy_, w_z_; // radial and axial waist
public:
	/**
	 * \brief       Implements Quian and Elson CEF
	 * \details    see article ...
	 *                  
	 * \param    Particle &p         Particle with which to calculate the CEF
	 * \return   A double with the probability ([0,1]) of detecting a photon at the position and with the orientation of the p particle.
	 */
	double getCollectionProbability(Particle &p);
	CEF_geometric();
	~CEF_geometric();
};

