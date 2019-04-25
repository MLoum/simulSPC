//
// Created by MatthieuL on 07/11/2018.
//

#include "Particle.h"
#include <stdio.h>
#include <random>

Particle::Particle(Experiment &exp) : rand_walk_gaussian_(0,1) {
    // initial position is random
    //http://www.cplusplus.com/reference/random/normal_distribution/

    exp_ = &exp;

    std::uniform_real_distribution<double> ini_pos_xy(0, exp_->solvent_.box_size_radial_);
    std::uniform_real_distribution<double> ini_pos_z(0, exp_->solvent_.box_size_axial_);


    //TODO from File
    Dx_ = 5E-3;
    Dy_ = 5E-3;
    Dz_ = 5E-3;

    x_ = ini_pos_xy(generator_);
    y_ = ini_pos_xy(generator_);
    z_ = ini_pos_z(generator_);
}

void Particle::initParticleParam(Experiment *exp)
{
    exp_ = exp;
    //rand_walk_gaussian_(0,1);

    std::uniform_real_distribution<double> ini_pos_xy(0, exp_->solvent_.box_size_radial_);
    std::uniform_real_distribution<double> ini_pos_z(0, exp_->solvent_.box_size_axial_);


    //TODO from File
    Dx_ = 5E-3;
    Dy_ = 5E-3;
    Dz_ = 5E-3;

    x_ = ini_pos_xy(generator_);
    y_ = ini_pos_xy(generator_);
    z_ = ini_pos_z(generator_);
}

void Particle::move() {

    //TODO Dx * exp_->time_step_ change unit to spare a multiplciation
    x_ += Dx_ * rand_walk_gaussian_(generator_) * exp_->time_step_;
    y_ += Dy_ * rand_walk_gaussian_(generator_) * exp_->time_step_;
    z_ += Dz_ * rand_walk_gaussian_(generator_) * exp_->time_step_;

    printf("%f %f %f\n", x_, y_, z_);

    // Periodic boundary condition
    if(x_ < 0)
    {x_ += exp_->solvent_.box_size_radial_;}
    else if(x_ >= exp_->solvent_.box_size_radial_)
    {x_ -= exp_->solvent_.box_size_radial_;}

    if(y_ < 0)
    {y_ += exp_->solvent_.box_size_radial_;}
    else if(y_ >= exp_->solvent_.box_size_radial_)
    {y_ -= exp_->solvent_.box_size_radial_;}

    if(z_ < 0)
    {z_ += exp_->solvent_.box_size_axial_;}
    else if(z_ >= exp_->solvent_.box_size_axial_)
    {z_ -= exp_->solvent_.box_size_axial_;}

}

Particle::~Particle() {

}

void Particle::light_matter_interraction() {

}

void Particle::inter_particle_interraction(){

}


Particle::Particle()
{

}
