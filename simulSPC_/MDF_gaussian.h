#ifndef SIMULSPC_MDF_GAUSSIAN_H
#define SIMULSPC_MDF_GAUSSIAN_H
#include "MDF.h"

class MDF;

class MDF_gaussian :
	public MDF
{
private:
	double w_xy_; //!< in nanometer
	double w_z_; //!< in nanometer
	double r_; //!< ratio of Axial waist (w_z_) over radial waist (w_xy_). This is an empirical value.
public:
	double get_detection_efficiency(Particle *p);
	MDF_gaussian(Experiment *experiment, double r=2);
	virtual void export_mdf();
	//MDF_gaussian();
	~MDF_gaussian();
};

#endif //SIMULSPC_MDF_GAUSSIAN_H

