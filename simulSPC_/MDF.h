#ifndef SIMULSPC_MDF_H
#define SIMULSPC_MDF_H

#include "PSF.h"
#include "CEF.h"
#include "Particle.h"
#include "Experiment.h"
#include <gsl/gsl_math.h>

class PSF;
class CEF;
class Particle;
class Experiment;

class MDF
{
protected :
	PSF *psf_;
	CEF *cef_;
	Experiment *exp_;
	double **mdf_;


public:
	virtual double get_detection_efficiency(Particle *p) = 0;
	virtual void export_mdf() = 0;
	MDF(Experiment *exp);
	MDF();
	~MDF();
};

#endif //SIMULSPC_MDF_H

