//for the M and N single site. Because the matrix is sparse, so construct a class especially.
#ifndef SINGLE_SUB_H
#define SINGLE_SUB_H
#include "SOP.h"
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;


class SingleSub
 {

 private:
        SOP _System;
        SOP _SysA;
        SOP _SysAdag;

 public:
        const SOP& System()const{return _System;};
        const SOP& SysA()const{return _SysA;};
        const SOP& SysAdag()const{return _SysAdag;};
         SingleSub(){};
         ~SingleSub(){};
         SingleSub(const Parameter& para)
        {
                SOP a(para, Annihilation), adag(para, Creation), iden(para, Iden);
                SOP sigmaz(para, SigmaZ), sigmap(para, SigmaP), sigmam(para, SigmaM);
                SOP sigmai(para, SigmaI);

                _SysA.Kron(a, sigmai);
                _SysAdag.Kron(adag, sigmai);


                SOP SysSigmaZ(iden, sigmaz), SysSigmap(iden, sigmap), SysSigmam(iden, sigmam);

                _System=_SysAdag*_SysA*para.omega0()+SysSigmap*SysSigmam*para.omegaq()
                +para.gr()*(_SysAdag*SysSigmam+SysSigmap*_SysA)
                +para.gcr()*(_SysAdag*SysSigmap+SysSigmam*_SysA);


}
         
 }; 


#endif // SINGLE_SUB_H
