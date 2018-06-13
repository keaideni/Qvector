#ifndef SUB_H
#define SUB_H

#include "OP.h"

using  namespace std; 
using namespace Eigen;
string itos(const int& i);
class Sub
{
private:
        int _Orbital;
        OP _System;
        OP _SysA;
        OP _SysAdag;
        OP _SysEye;
        OP _SysA1;
        OP _SysAdag1;

public:
        static int nmax;
        const int& Orbital()const{return _Orbital;};
        const OP& System()const{return _System;};
        const OP& SysA()const{return _SysA;};
        const OP& SysAdag()const{return _SysAdag;};
        const OP& SysA1()const{return _SysA1;};
        const OP& SysAdag1()const{return _SysAdag1;};
        const OP& SysEye()const{return _SysEye;};


        Sub(){};
        ~Sub(){};
        Sub(const Sub& a):
        _Orbital(a._Orbital),
        _System(a._System),
        _SysA(a._SysA),
        _SysAdag(a._SysAdag),
        _SysA1(a._SysA1),
        _SysAdag1(a._SysAdag1)
        {}

        Sub(const Parameter& para, const int& orbital);
        //Sub(const Parameter& para, const Sub& SubL, const Sub& SubR, const int& orbital);

        //const Sub& operator=(const Sub& a);

        //void Trunc(const MatrixXd& U);
        //void Save()const;
        //void Read(const int& orbital);
        void Show()const;

        //void ChangeOrbital(const int& orb){_Orbital=orb;};

};



#endif // SUB_H
