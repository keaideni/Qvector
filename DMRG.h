#ifndef DMRG_H
#define DMRG_H
#include "SuperEnergy.h"



class DMRG
{
        QWave InitWave;
        OP InitOPWave;
        Sub Sys;
        Sub Env;
        Sub m;
        Sub n;
        
        ofstream SaveAll;
        double _FEnergy;
        double _Entropy;
	double _Excited;
        const Parity& pari;
public:
        const double& FEnergy()const{return _FEnergy;};
        const double& Entropy()const{return _Entropy;};
        const double& Excited()const{return _Excited;};
        


        //DMRG(){};
        ~DMRG(){};

        DMRG(Parameter& para, const Parity& _pari);

        void Initialize(const Parameter& para, const int& dir, const int& Gdir, int& OS, int & OE);



        //=================periodic condition==================
        void BuildUp(Parameter& para, int& OS, int& OE);
        
        //================sweep======================
        void Sweep(Parameter& para, int& OS, int& OE);
        void CalcuEnergy(Parameter& para, int& OS, int& OE, const int& dir, const int& Gdir);
        //===========one site sweep===================
        //void OneSiteSweep(Parameter& para, int& OS, int& OE);
        
        
        
};









#endif // DMRG_H
