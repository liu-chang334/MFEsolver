#include <iostream> 

#include <elem/C3D8.h>
#include <elem/MaterialPoint3D.h>
#include <cons/IsoElastic.h>




int main() {
    C3D8 e(1);
    e.mat = new IsoElastic(1000, 0.3, 1000, 3);
    Eigen::VectorXd dstrain = Eigen::VectorXd::Zero(6);
    dstrain << 0.1, 0.1, 0.1, 0, 0, 0;
    for (int i = 0; i < e.mp.size(); i++) {
        std::cout << e.mp[i].strain.transpose() << std::endl;
        std::cout << e.mp[i].stress.transpose() << std::endl;
        e.mat->updatemp(dstrain, e.mp[i]);
    }
    for (int i = 0; i < e.mp.size(); i++) {
        std::cout << e.mp[i].strain.transpose() << std::endl;
        std::cout << e.mp[i].stress.transpose() << std::endl;
    }

    return 0;
}
