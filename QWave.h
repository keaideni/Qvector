#ifndef Q_WAVE_H
#define Q_WAVE_H
#include<vector>
#include "Sub.h"
#include "SingleSub.h"

struct dimension
{
        int SDim, EDim, MDim, NDim;
};

enum Parity
{
        Positive, Negative
};

class QWave
{
private:

        vector<vector<vector<MatrixXd>>> _Wave;//i<order<m<n,mat>>>
        vector<dimension> _Dim;
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
        const vector<unordered_map<string, int>> DimOES()const //dimension of parity for each site.
        {
                vector<unordered_map<string, int>> tempdim(4);//to store the dimension for sys(0)... Env(3) of different parity.
                for(int i=0; i<16; ++i)
                {
                        vector<string> temps(n2s(i));
                        for(int k=0; i<4; ++k)
                        {
                                auto it=tempdim.at(k).find(temps.at(k));
                                if(it==tempdim.at(k).end())
                                {
                                        int tempd;
                                        switch (k)
                                        {
                                                case 0:
                                                {
                                                        tempd=_Dim.at(i).SDim;
                                                        break;
                                                }
                                                case 1:
                                                {
                                                        tempd=_Dim.at(i).MDim;
                                                        break;
                                                }
                                                case 2:
                                                {
                                                        tempd=_Dim.at(i).NDim;
                                                        break;
                                                }
                                                case 3:
                                                {
                                                        tempd=_Dim.at(i).EDim;
                                                        break;
                                                }
                                        }
                                        tempdim.at(k).insert(pair<string, int>(temps.at(k), tempd));
                                }
                        }
                }
                return tempdim;
        }


public:

        const vector<vector<vector<MatrixXd>>>& Wave()const
        {return _Wave;};
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
        //void NSME(MatrixXd& A)const;


        


        void Wave2f(vector<double>& f, const Parity& pari)const;
        void f2Wave(const vector<double>& f, const Parity& pari);
//        void f2Wave(const VectorXd& f);
/*
        const QWave& operator=(const QWave& a)
        {
                _Wave=a.Wave();
                return *this;
        }
        void Norm();
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
