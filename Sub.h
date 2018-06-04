#ifndef SUB_H
#define SUB_H
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include "Parameter.h"
#include<sstream>
#include<string>


using  namespace std; 
using namespace Eigen;
string itos(const int& i);
class Sub
{
private:
        int _Orbital;
        MatrixXd _System;
        MatrixXd _SysA;
        MatrixXd _SysAdag;
        MatrixXd _SysEye;
        MatrixXd _SysA1;
        MatrixXd _SysAdag1;

        void MatrixSave(const MatrixXd& A, ofstream& outfile)const;
        void MatrixRead(MatrixXd& A, ifstream& infile);
        void Kron(MatrixXd& ab,const MatrixXd& a, const MatrixXd& b);
public:
        static int nmax;
        const int& Orbital()const{return _Orbital;};
        const MatrixXd& System()const{return _System;};
        const MatrixXd& SysA()const{return _SysA;};
        const MatrixXd& SysAdag()const{return _SysAdag;};
        const MatrixXd& SysA1()const{return _SysA1;};
        const MatrixXd& SysAdag1()const{return _SysAdag1;};
        const MatrixXd& SysEye()const{return _SysEye;};


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
        Sub(const Parameter& para, const Sub& SubL, const Sub& SubR, const int& orbital);

        const Sub& operator=(const Sub& a);

        void Trunc(const MatrixXd& U);
        void Save()const;
        void Read(const int& orbital);
        void Show()const;

        void ChangeOrbital(const int& orb){_Orbital=orb;};

};



#endif // SUB_H
