#include "MDF.h"



double MDF::get_detection_efficiency(Particle *p)
{
	double psf = psf_->get_probability_absorption(p);
	double cef = cef_->get_collection_probability(p);
	return cef * psf;
}

void MDF::export_mdf()
{
}

void MDF::explore_mdf()
{
}

MDF::MDF(Experiment *exp)
{
	exp_ = exp;
}

MDF::MDF()
{
}


MDF::~MDF()
{
}
