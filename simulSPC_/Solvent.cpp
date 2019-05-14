//
// Created by MatthieuL on 07/11/2018.
//

#include "Solvent.h"
#include "Particle.h"
#include "Experiment.h"

Solvent::Solvent(Experiment *exp)
{
	exp_ = exp;


}

void Solvent::create_velocity_map()
{
	int nb_step_y = (int)(box_size_radial_ / exp_->space_step_) + 1;
	int nb_step_z = (int)(box_size_axial_ / exp_->space_step_) + 1;

	velocity_map_ = new double*[nb_step_y];
	for (int i = 0; i < nb_step_y; i++)
		velocity_map_[i] = new double[nb_step_z];

	for (int y = 0; y < nb_step_y; y++)
	{
		for (int z = 0; z < nb_step_z; z++)
			//TODO rectangular Poiseuille flow. See for example https://en.wikipedia.org/wiki/Hagen%E2%80%93Poiseuille_equation#Poiseuille_flow_through_some_non-circular_cross-sections[10]
			velocity_map_[y][z] = 0;
	}
}

double Solvent::get_convection(Particle *p)
{
	return velocity_map_[(int)(p->y_ / exp_->space_step_)][(int)(p->z_ / exp_->space_step_)];
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