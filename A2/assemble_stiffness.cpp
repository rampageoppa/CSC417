#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) {

    std::vector< Eigen::Triplet<double> > temp_triplet;
    Eigen::Matrix66d H;

    for(int s = 0; s < E.rows(); s++){
        Eigen::Vector3d q0 = q.segment<3>(3.0 * E(s, 0));
        Eigen::Vector3d q1 = q.segment<3>(3.0 * E(s, 1));

        d2V_spring_particle_particle_dq2(H, q0, q1, l0[s], k);
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                temp_triplet.push_back(Eigen::Triplet<double>(3.0 * E(s,0) + i, 3 * E(s, 0) + j, H(i, j)));
                temp_triplet.push_back(Eigen::Triplet<double>(3.0 * E(s,0) + i, 3 * E(s, 1) + j, H(i, j + 3)));
                temp_triplet.push_back(Eigen::Triplet<double>(3.0 * E(s,1) + i, 3 * E(s, 0) + j, H(i + 3, j)));
                temp_triplet.push_back(Eigen::Triplet<double>(3.0 * E(s,1) + i, 3 * E(s, 1) + j, H(i + 3, j + 3)));
            }
        }
    }
    K.resize(q.size(), q.size());
    K.setZero();
    
    K.setFromTriplets(temp_triplet.begin(), temp_triplet.end());
};