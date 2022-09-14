#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

    M.resize(q.size(), q.size());
    M.setIdentity();
    for(int i = 0; i < q.size(); i++){
        M.insert(i,i) *= mass;
    }
}
