#ifndef Q_WAVE_H
#define Q_WAVE_H
#include<vector>
#include "Sub.h"
#include "SingleSub.h"
#include <ctime>
struct dimension
{
        int SDim, EDim, MDim, NDim;
};


class QWave
{
private:


        //========================reshape===============================

        const int order(const string& a1, const string& a2, const string& a3, const string& a4)const
        {
                return s2n(a1)*pow(2, 3)+s2n(a2)*pow(2, 2)+s2n(a3)*pow(2,1)+s2n(a4);
        }

        const int s2n(const string& s)const
        {
                if(s=="positive")
                        return 0;
                else
                        return 1;
        }

        const vector<string> n2s(const int& n)const
        {
                vector<string> haha;
                if(n/8==1)haha.push_back("negative");
                else haha.push_back("positive");

                if((n%8)/4==1)haha.push_back("negative");
                else haha.push_back("positive");
                
                if(((n%8)%4)/2==1)haha.push_back("negative");
                else haha.push_back("positive");
                
                if(((n%8)%4)%2==1)haha.push_back("negative");
                else haha.push_back("positive");

                return haha;
        }


        //=====================================================
        ////Wave=Sys*Wave and only used for no parity changed system operator;
        void SysOPWave(const OP& Sys, const Parity& pari);
        //_Wave=Sys*wave;
        void SysOPWave(const OP& Sys, const QWave& wave, const Parity& pari);
        //_wave+=wave*Env.
        void EnvOPWave(const OP& Env, const QWave& wave, const Parity& pari);
        //_wave+=M*wave.
        void MOPWave(const SOP& M, const QWave& wave, const Parity& pari);
        void NOPWave(const SOP& N, const QWave& wave, const Parity& pari);
        //_Wave+=g*Sys*wave*M;
        void SysMOPWave(const OP& Sys, const SOP& M, const QWave& Wave, const double& g, const Parity& pari);      
        //_Wave+=g*M*wave*Env;
        void EnvMOPWave(const OP& Env, const SOP& M, const QWave& Wave, const double& g, const Parity& pari);       
        //_Wave+=g*wave*N*Env;
        void EnvNOPWave(const OP& Sys, const SOP& M, const QWave& Wave, const double& g, const Parity& pari);       
        //_Wave+=g*Sys*wave*N;
        void SysNOPWave(const OP& Sys, const SOP& M, const QWave& Wave, const double& g, const Parity& pari);        

        void add(const QWave& wavei, const Parity& pari);
//=====================================================
        const vector<unordered_map<string, int>> DimOES()const;//dimension of parity for each site.
        

        const string QAdd(const string& a, const string& b)const;

        const string OppS(const string& s)const;//return the opposite of s.

public:
        const vector<vector<vector<MatrixXd>>>& Wave()const
        {return _Wave;};
        const vector<dimension>& Dim()const{return _Dim;};
        //QWave(){};
        ~QWave(){};
        QWave(const QWave&a):
        _Dim(a._Dim),
        _Wave(a._Wave)
        {};
        QWave(const Sub& sys, const SingleSub& m, const SingleSub& n, const Sub& env);
        
        //_Wave=H*wave; befor this must copy wave to _Wave
        void Hamiltanian(const Sub& Sys, const SingleSub& M, const SingleSub& n, const Sub& Env, const QWave& wave, const Parameter& para, const Parity& pari);
        QWave(const Sub& Sys, const SingleSub& M, const SingleSub& n, const Sub& Env, const QWave& wave, const Parameter& para, const Parity& pari);


        //======== Trunc=================================
        void Wave2SMEN(OP& A, const Parity& pari)const;
        void Wave2NSME(OP& A, const Parity& pari)const;
        void NSME2Wave(const OP& A, const Parity& pari);
        void SMEN2Wave(const OP& A, const Parity& pari);
        //void NSME(MatrixXd& A)const;


        


        void Wave2f(vector<double>& f, const Parity& pari)const;
        void f2Wave(const vector<double>& f, const Parity& pari);

        int waveDim(const Parity& pari)const;
//        void f2Wave(const VectorXd& f);

        const QWave& operator=(const QWave& a)
        {
                _Dim=a._Dim;
                _Wave=a.Wave();
                return *this;
        }
        /*void Norm();
        void Show()const;*/
        /*void Random()
        {
                for(int im=0; im<Dm; ++im)
                {
                        for(int in=0; in<Dn; ++in)
                        {
                                _Wave.at(im).at(in)
                                =MatrixXd::Random(DSys, DEnv);
                        }
                }
        }*/


};



#endif // Q_WAVE_H
