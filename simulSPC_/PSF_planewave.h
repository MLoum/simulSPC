#pragma once
#include "PSF.h"


/**
 * \brief      State of the art PSF evaluation based on the interference of the incident  laser plane wave  around the focus point.
 *
 * \details		Give some article.
 * \sa	PSF.h
 */
class PSF_planewave :
	public PSF
{
public:
	PSF_planewave();
	~PSF_planewave();
};

