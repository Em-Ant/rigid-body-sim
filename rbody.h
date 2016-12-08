
/********************************************
*  Header for Rigid Body Class              *
*  Version : 0.1                            *
*  Author: Emant                            *
*  Year : 2014                              *
*********************************************/

#ifndef RBODY_H
#define RBODY_H

#include "defs.h"

class Rbody final
{
    public:
        Rbody(FLOAT mass,
              vect diagonal_inertia_tensor,
              vstate init_state);

        //virtual ~Rbody();

        // Interface
        matrix& integrate(FLOAT dt);
        void get_omega(vect& v);
        void get_angular_momentum(vect& L);
        void clear_force();
        void clear_torque();
        void add_force_world_cs(vect const& force, vect const& pos);
        void add_force_body_cs(vect const& force, vect const& pos);
        void add_g_force_world_cs(vect const& force);
        void add_g_force_body_cs(vect const& force);
        void add_torque_world_cs(vect const& torque);
        void add_torque_body_cs(vect const& torque);

    private:
        // BODY DATA
        FLOAT   imass;          //Inverse Mass m^-1
        vect    iinertia;       //Inverse Inertia Tensor I^-1 in Body Central Principal Frame
        // STATE VECTOR
        vstate SV;              // SV = {vect X, vect P, quat R, vect L};
                                //  X = Position of Mass Center, P = Momentum,
                                //  R = Rotation Quaternion, L = Angular Momentum
        // DERIVATE STATE VECTOR
        vstate DSV[3];          // DSV = {vect V, vect F, quat dR, vect T};
                                // Velocity V = P/m, Input Forces F: dP/dt = F,
                                // dR/dt = 1/2*omega*R (quat omega = [0, vect omega]),
                                // Input Torque T : dL/dt = T
                                // The array stores values from two preceding time steps.
                                // Needed by integration routine.
        // FLAGS & INDEXES
        short started;          // Needed by Integration Routine to calc  working condition
        short k,l,m,tmp;

        // Aux Data
        vect    homega;         // Half Angular velocity : omega = I^-1 * L,
                                //   homega = 0.5*omega. Needed to calculate dR/dt
        vect6   wrld_iin;       // Inverse Inertia Tensor in World Frame.
                                //   Needed for updating L
        matrix  rpmatrix;       // OpenGL 4x4 Rotation + Position Matrix

        // Functions
        void update_dstate(vstate &DSV_);
};

#endif // RBODY_H
