// #pragma once
#include <vector>


#include <elem/ElemBase.h>
#include <elem/MaterialPoint3D.h>
#include <cons/MaterialBase.h>

class C3D8 : public ElemBase{
public:
    MaterialBase *mat;
    std::vector<MaterialPoint3D> mp;

public:
    C3D8(int eid);
    ~C3D8() = default;



};