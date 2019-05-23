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
	virtual double get_detection_efficiency(Particle *p);
	virtual void export_mdf();
	virtual void explore_mdf();
	MDF(Experiment *exp);
	MDF();
	~MDF();
	void set_psf(PSF *psf) { psf_ = psf;}
	void set_cef(CEF *cef) { cef_ = cef; }
};

#endif //SIMULSPC_MDF_H

