#ifndef Q_WAVE_H
#define Q_WAVE_H
#include<vector>
#include "Sub.h"
#include "SingleSub.h"


class QWave
{
private:

        vector<vector<MatrixXd>> _Wave;//<m<n,mat>>
        const int DSys, DEnv, Dm, Dn;
        //========================reshape===============================




public:

        const vector<vector<MatrixXd>>& Wave()const
        {return _Wave;};
        //QWave(){};
        ~QWave(){};
        QWave(const QWave&a):
        DSys(a.DSys),
        DEnv(a.DEnv),
        Dm(a.Dm),
        Dn(a.Dn),
        _Wave(a._Wave)
        {};
        QWave(const int& sys, const int& m, const int& n, const int& env):
        DSys(sys),
        DEnv(env),
        Dm(m),
        Dn(n)
        {
                MatrixXd temp(MatrixXd::Zero(DSys, DEnv));
                vector<MatrixXd> tempv;
                for(int i=0; i<Dn; ++i)
                {
                        tempv.push_back(temp);
                }
                for(int i=0; i<Dm; ++i)
                {
                        _Wave.push_back(tempv);
                }
        };

        void SMEN2Wave(const MatrixXd& A);
        void SMEN(MatrixXd& A)const;
        void NSME(MatrixXd& A)const;
        void TruncL(MatrixXd& truncU, const int& D)const;

//=====================================================
        void SysOPWave(const MatrixXd& Sys);//Wave=Sys*Wave;
        //_Wave+=Sys*wave;
        void SysOPWave(const MatrixXd& Sys, const QWave& Wave);        
        //_wave+=wave*Env.
        void EnvOPWave(const MatrixXd& Env, const QWave& wave);
        //_wave+=M*wave.
        void MOPWave(const SpMat& M, const QWave& wave);
        //M*_Wave;
        const QWave MOPWave(const SpMat& M, const double&)const;
        //_wave+=wave*N
        void NOPWave(const SpMat& N, const QWave& wave);
        //_Wave*N;
        const QWave NOPWave(const SpMat& N, const double&)const;
        void add(const QWave& wave);
//=====================================================

        


        void Wave2f(vector<double>& f)const;
        void f2Wave(const vector<double>& f);
        void f2Wave(const VectorXd& f);

        const QWave& operator=(const QWave& a)
        {
                _Wave=a.Wave();
                return *this;
        }
        void Norm();
        void Show()const;
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
