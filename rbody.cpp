
/********************************************
*  Rigid Body Class                         *
*  Version : 0.1                            *
*  Author: Emant                            *
*  Year : 2014                              *
*********************************************/

#include "rbody.h"
#include <math.h>

Rbody::Rbody( FLOAT mass,
              vect diagonal_inertia_tensor,
              vstate init_state)
{
    imass = 1/mass;
    for(int i = 0; i<3; ++i)
        iinertia[i] = 1/diagonal_inertia_tensor[i];
    for(int i = 0; i<13; ++i)
        SV[i] = init_state[i];
    FLOAT cs = cos(SV[6]*PI/360);
    FLOAT sn = sin(SV[6]*PI/360);
    SV[6] = cs;
    for(int i = 0; i<3; ++i)
        SV[7+i] *= sn;
    rpmatrix[3]  = 0;
    rpmatrix[7]  = 0;
    rpmatrix[11] = 0;
    rpmatrix[15] = 1;
    update_dstate(DSV[0]);
    started = 0;
    k=0;l=1;m=0;
}

//Rbody::~Rbody(){}

void Rbody::get_omega(vect& v)
{
    for(int i = 0; i<3; ++i)
        v[i] = 2*homega[i];
}

void Rbody::get_angular_momentum(vect& L)
{
    for(int i = 0; i<3; ++i)
        L[i] = SV[i+10];
}

void Rbody::update_dstate(vstate &DSV_)
{
    // Update Velocity
     DSV_[0] = imass*SV[3]; DSV_[1] = imass*SV[4]; DSV_[2] = imass*SV[5];
    // Normalize Quaternion
        FLOAT den = 0;
    for(int i = 6; i<10 ;++i)
        den += SV[i]*SV[i];
    den = 1.0/sqrt(den);
    for(int i = 6; i<10 ;++i)
        SV[i] = SV[i]*den;
    // Update Rot Matrix;
    rpmatrix[0] = 1 - 2*(SV[8]*SV[8]+SV[9]*SV[9]);
    rpmatrix[1] = 2*(SV[7]*SV[8] + SV[6]*SV[9]);
    rpmatrix[2] = 2*(SV[7]*SV[9] - SV[6]*SV[8]);
    rpmatrix[4] = 2*(SV[7]*SV[8] - SV[6]*SV[9]);
    rpmatrix[5] = 1 - 2*(SV[7]*SV[7]+SV[9]*SV[9]);
    rpmatrix[6] = 2*(SV[8]*SV[9] + SV[6]*SV[7]);
    rpmatrix[8] = 2*(SV[7]*SV[9] + SV[6]*SV[8]);
    rpmatrix[9] = 2*(SV[8]*SV[9] - SV[6]*SV[7]);
    rpmatrix[10] = 1 - 2*(SV[7]*SV[7]+SV[8]*SV[8]);
    rpmatrix[12] = SV[0]; rpmatrix[13] = SV[1]; rpmatrix[14] = SV[2];
    // Calculate World Inverse Inertia Tensor. Is a symmatric tensor, so only 6 term are needed
    wrld_iin[0] = rpmatrix[0]*rpmatrix[0]*iinertia[0] + rpmatrix[4]*rpmatrix[4]*iinertia[1] + rpmatrix[8]*rpmatrix[8]*iinertia[2];
    wrld_iin[1] = rpmatrix[1]*rpmatrix[1]*iinertia[0] + rpmatrix[5]*rpmatrix[5]*iinertia[1] + rpmatrix[9]*rpmatrix[9]*iinertia[2];
    wrld_iin[2] = rpmatrix[2]*rpmatrix[2]*iinertia[0] + rpmatrix[6]*rpmatrix[6]*iinertia[1] + rpmatrix[10]*rpmatrix[10]*iinertia[2];
    wrld_iin[3] = rpmatrix[0]*rpmatrix[1]*iinertia[0] + rpmatrix[4]*rpmatrix[5]*iinertia[1] + rpmatrix[8]*rpmatrix[9]*iinertia[2];
    wrld_iin[4] = rpmatrix[0]*rpmatrix[2]*iinertia[0] + rpmatrix[4]*rpmatrix[6]*iinertia[1] + rpmatrix[8]*rpmatrix[10]*iinertia[2];
    wrld_iin[5] = rpmatrix[1]*rpmatrix[2]*iinertia[0] + rpmatrix[5]*rpmatrix[6]*iinertia[1] + rpmatrix[9]*rpmatrix[10]*iinertia[2];
    //omega = [I^-1_world]*[L], [I^-1_world] = [R][I^-1][R^-1], homega = 0.5*omega;
    homega[0] = 0.5*(wrld_iin[0]*SV[10] + wrld_iin[3]*SV[11] + wrld_iin[4]*SV[12]);
    homega[1] = 0.5*(wrld_iin[3]*SV[10] + wrld_iin[1]*SV[11] + wrld_iin[5]*SV[12]);
    homega[2] = 0.5*(wrld_iin[4]*SV[10] + wrld_iin[5]*SV[11] + wrld_iin[2]*SV[12]);
    // Calculate dR/dt = homega*R ; homega = [0, 0.5*omega]
    // Quaternion Mult: q0*q1 = [s0 s1 - v0*v1, s0 v1 + s1 v0 + v0 x v1]. !!! s0 scalar term is 0 in homega !!!
    DSV_[6] =  -(homega[0]*SV[7] + homega[1]*SV[8] + homega[2]*SV[9]);
    DSV_[7] =  homega[0]*SV[6] + homega[1]*SV[9] - homega[2]*SV[8];
    DSV_[8] =  homega[1]*SV[6] + homega[2]*SV[7] - homega[0]*SV[9];
    DSV_[9] =  homega[2]*SV[6] + homega[0]*SV[8] - homega[1]*SV[7];
}

matrix& Rbody::integrate(FLOAT dt)
{

    if(started == 2)                    //ADAMS-BASHFORTH 3-steps during normal working
    {
        for(int i  = 0; i<13; ++i)
        {
            SV[i] += dt/12*(23*DSV[k][i]-16*DSV[l][i]+5*DSV[m][i]);
        }
        update_dstate(DSV[m]);
        tmp=k; k=m; m=l; l=tmp;
    }
    else if (started == 1)              //ADAMS-BASHFORTH 2-steps : second integration step
    {
         for(int i  = 0; i<13; ++i)
        {
            SV[i] += dt*0.5*(3*DSV[1][i]-DSV[1][i]);
        }
        update_dstate(DSV[2]);
        started = 2; k=2;
    }
    else                                //Explicit EULER : First integration step.
    {
        for(int i  = 0; i<13; ++i)
        {
            SV[i] += DSV[0][i]*dt;
        }
        update_dstate(DSV[1]);
        started = 1; k=1;
    }
    return rpmatrix;
}

void Rbody::clear_force()
{
    for(int i = 3; i<6; ++i)
        DSV[k][i] = 0;
}

void Rbody::clear_torque()
{
    for(int i = 10; i<13; ++i)
        DSV[k][i] = 0;
}

void Rbody::add_force_world_cs(vect const& force,
                          vect const& pos )
{
    for(int i = 3; i<6; ++i)
    {
        DSV[k][i] += force[i];
    }
    DSV[k][10] += (pos[1]-SV[1])*force[2] - (pos[2]-SV[2])*force[1];
    DSV[k][11] += (pos[2]-SV[2])*force[0] - (pos[0]-SV[0])*force[2];
    DSV[k][12] += (pos[0]-SV[0])*force[1] - (pos[1]-SV[1])*force[0];
}

void Rbody::add_force_body_cs(vect const& force,
                          vect const& pos )
{
    FLOAT t0,t1,t2;
    t0 = pos[1]*force[2]-pos[2]*force[1];
    t1 = pos[2]*force[0]-pos[0]*force[2];
    t2 = pos[0]*force[1]-pos[1]*force[0];
    DSV[k][3] += rpmatrix[0]*force[0] + rpmatrix[4]*force[1] + rpmatrix[8]*force[2];
    DSV[k][4] += rpmatrix[1]*force[0] + rpmatrix[5]*force[1] + rpmatrix[9]*force[2];
    DSV[k][5] += rpmatrix[2]*force[0] + rpmatrix[6]*force[1] + rpmatrix[10]*force[2];
    DSV[k][10] += rpmatrix[0]*t0 + rpmatrix[4]*t1 + rpmatrix[8]*t2;
    DSV[k][11] += rpmatrix[1]*t0 + rpmatrix[5]*t1 + rpmatrix[9]*t2;
    DSV[k][12] += rpmatrix[2]*t0 + rpmatrix[6]*t1 + rpmatrix[10]*t2;
}

void Rbody::add_g_force_world_cs(vect const& force)
{
    for(int i = 3; i<6; ++i)
    {
        DSV[k][i] += force[i];
    }
}

void Rbody::add_g_force_body_cs(vect const& force)
{
    DSV[k][3] += rpmatrix[0]*force[0] + rpmatrix[4]*force[1] + rpmatrix[8]*force[2];
    DSV[k][4] += rpmatrix[1]*force[0] + rpmatrix[5]*force[1] + rpmatrix[9]*force[2];
    DSV[k][5] += rpmatrix[2]*force[0] + rpmatrix[6]*force[1] + rpmatrix[10]*force[2];
}

void Rbody::add_torque_world_cs(vect const& torque)
{
    for(int i = 10; i<13; ++i)
    {
        DSV[k][i] += torque[i];
    }
}

void Rbody::add_torque_body_cs(vect const& torque)
{
    DSV[k][10] += rpmatrix[0]*torque[0] + rpmatrix[4]*torque[1] + rpmatrix[8]*torque[2];
    DSV[k][11] += rpmatrix[1]*torque[0] + rpmatrix[5]*torque[1] + rpmatrix[9]*torque[2];
    DSV[k][12] += rpmatrix[2]*torque[0] + rpmatrix[6]*torque[1] + rpmatrix[10]*torque[2];
}
