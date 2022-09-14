//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd k_1;
    Eigen::VectorXd k_2;
    Eigen::VectorXd k_3;
    Eigen::VectorXd k_4;

    Eigen::VectorXd force_k1;
    force(force_k1, q, qdot);
    k_1 = force_k1 / mass;

    Eigen::VectorXd force_k2;
    Eigen::VectorXd temp_1 = dt * k_1/2;
    force(force_k2, q, temp_1);
    k_2 = force_k2 / mass;

    Eigen::VectorXd force_k3;
    Eigen::VectorXd temp_2 = dt * k_2/2;
    force(force_k3, q, temp_2);
    k_3 = force_k3 / mass;

    Eigen::VectorXd force_k4;
    Eigen::VectorXd temp_3 = dt * k_3;
    force(force_k4, q, temp_3);
    k_4 = force_k4 / mass;

    qdot = qdot + dt * (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6;
    q = q + dt * qdot;
}