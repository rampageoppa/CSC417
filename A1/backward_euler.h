#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force, STIFFNESS &stiffness) {
    Eigen::VectorXd temp_f;
    Eigen::MatrixXd temp_k;
    force(temp_f, q, qdot);
    stiffness(temp_k, q, qdot);
    qdot = (qdot + dt * temp_f / mass) / (1 - dt * dt * temp_k(0, 0) / mass);
    q = q + dt * qdot;
}