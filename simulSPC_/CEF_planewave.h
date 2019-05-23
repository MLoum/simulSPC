#pragma once
#include "CEF.h"


class CEF_planewave :
	public CEF
{
public:
	CEF_planewave();
	CEF_planewave(Experiment *exp);
	~CEF_planewave();
	virtual double get_collection_probability(Particle *p);
};

