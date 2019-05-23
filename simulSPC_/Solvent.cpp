//
// Created by MatthieuL on 07/11/2018.
//

#include "Solvent.h"
#include "Particle.h"
#include "Experiment.h"
constexpr double PI = 3.1415159265;

Solvent::Solvent(Experiment *exp)
{
	exp_ = exp;
}

void Solvent::create_velocity_map()
{
	int nb_step_y = (int)(box_size_radial_ / exp_->space_step_) + 1;
	int nb_step_z = (int)(box_size_axial_ / exp_->space_step_) + 1;

	double h, l, y_m, z_m, G, Bn, sumflow;

	h = box_size_radial_ * 1E-9;
	l = box_size_axial_ * 1E-9;

	sumflow = pow(h, 3) * l / (12 * viscosity_);
	for (int n = 1; n <= 20; n++)
	{
		Bn = (2 * n - 1) * PI / h;
		sumflow -= ((16 * pow(h, 4)) / (pow(PI, 5) * viscosity_ * pow((2 * n - 1), 5))) * ((cosh(Bn * l) - 1) / sinh(Bn * l));
	}

	G = meanFlow_ * 1E-9 / (60 * sumflow);

	velocity_map_ = new double*[nb_step_y];
	for (int i = 0; i < nb_step_y; i++)
		velocity_map_[i] = new double[nb_step_z];

	for (int y = 0; y < nb_step_y; y++)
	{
		y_m = (((float)y) / ((float)nb_step_y)) * h;
		for (int z = 0; z < nb_step_z; z++)
		{
			z_m = (((float)z) / ((float)nb_step_z)) * l;
			//TODO rectangular Poiseuille flow. See for example https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation#Poiseuille_flow_through_some_non-circular_cross-sections[10]
			velocity_map_[y][z] = (G * y_m * (h - y_m)) / (2 * viscosity_);
			for (int n = 1; n <= 20; n++)
			{
				Bn = (2 * n - 1) * PI / h;
				velocity_map_[y][z] -= (4 * G*h*h / (viscosity_ * pow(PI, 3) * pow((2 * n - 1), 3))) * ((sinh(Bn*z_m) + sinh(Bn*(l - z_m))) / (sinh(Bn*l))) * sin(Bn*y_m);
			}
		}
	}
}

double Solvent::get_convection(Particle *p)
{
	return velocity_map_[(int)(p->y_ / exp_->space_step_)][(int)(p->z_ / exp_->space_step_)] * exp_->time_step_;
}

Solvent::~Solvent() {
	int nb_step_y = (int)(box_size_radial_ / exp_->space_step_);
	if (velocity_map_ != NULL)
	{
		for (int i = 0; i < nb_step_y; ++i) {
			delete[] velocity_map_[i];
		}
		delete[] velocity_map_;
	}

}