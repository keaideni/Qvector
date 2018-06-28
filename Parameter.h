#ifndef PARAMETER_H
#define PARAMETER_H
#include<fstream>
#include<string>
#include<iostream>

enum OpType{
	Creation, Annihilation, Iden, SigmaZ, SigmaP, SigmaM, SigmaI
};
class Parameter
{
private:
        double _omega0, _omegaq,_gr, _gcr, _Jr, _Jcr;
        int _D;
        int _LatticeSize;
        int _nmax;
public:
        const double& omega0()const{return _omega0;};
        const double& omegaq()const{return _omegaq;};
        const double& gr()const{return _gr;};
        const double& gcr()const{return _gcr;};
        const double& Jr()const{return _Jr;};
        const double& Jcr()const{return _Jcr;};
        const int& LatticeSize()const{return _LatticeSize;};
        const int& D()const{return _D;};
        const int& nmax()const{return _nmax;};

        void ChangeD(const int& dd){_D=dd;};
        //Parameter();
        ~Parameter(){};
        double Energy;

        Parameter()
        {
                std::ifstream infile("Parameter");
                if(!infile)
                {
                        std::cerr<<"The file Parameter doesn't exit!"<<std::endl;
                        exit(true);
                }

                std::string temp;
                infile>>temp>>_nmax>>temp>>_D>>temp>>_LatticeSize>>temp>>_omega0>>
                        temp>>_omegaq>>temp>>_gr>>temp>>_gcr
                >>temp>>_Jr>>temp>>_Jcr;
        };
        
};

std::string itos(const int& i);
#endif // PARAMETER_H
