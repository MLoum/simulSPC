#include "CEF_planewave.h"



CEF_planewave::CEF_planewave()
{
}

CEF_planewave::CEF_planewave(Experiment * exp)
{
	exp_ = exp;
}


CEF_planewave::~CEF_planewave()
{
}

double CEF_planewave::get_collection_probability(Particle * p)
{
	//TODO
	return 0.0;
}
