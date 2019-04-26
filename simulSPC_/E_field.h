#pragma once
#include <gsl/gsl_complex.h>

typedef struct
{
	gsl_complex Ex;
	gsl_complex Ey;
	gsl_complex Ez;
}E_field;