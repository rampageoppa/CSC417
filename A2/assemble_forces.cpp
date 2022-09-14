#include <assemble_forces.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) {
    f.resize(q.size());
    f.setZero();

    // https://eigen.tuxfamily.org/dox/group__TutorialBlockOperations.html
    // E: the spring connectivity matrix
    // need to consider gravity?

    for(int i = 0; i < E.rows(); i++){
        Eigen::Vector3d q0 = q.segment<3>(3.0 * E(i, 0));
        Eigen::Vector3d q1 = q.segment<3>(3.0 * E(i, 1));
        Eigen::Vector6d temp_force;

        dV_spring_particle_particle_dq(temp_force, q0, q1, l0[i], k);

        f.segment<3>(3 * E(i, 0)) -= temp_force.head<3>();
        f.segment<3>(3 * E(i, 1)) -= temp_force.tail<3>();
    }
};