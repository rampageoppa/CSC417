#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::VectorXd &dV, const Eigen::VectorXd &q, double stiffness) {
    dV.resize(1);
    
    // (d/dq)1/2kq^2 => kq
    dV = (q * stiffness);
}