#include <elem/C3D8.h>


C3D8::C3D8(int eid) : ElemBase(eid, "C3D8") {
    mp.resize(8);
    for (int i = 0; i < mp.size(); i++) mp[i] = MaterialPoint3D(i);
}