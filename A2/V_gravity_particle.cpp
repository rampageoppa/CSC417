#include <V_gravity_particle.h>

void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    
    // V = mgh
    // 加负号？
    V = - mass * q.transpose() * g;
}