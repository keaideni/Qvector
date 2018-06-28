#ifndef SUPER_H
#define SUPER_H
#include "QWave.h"

class Super
{
private:
        double sigma_;
        //double Jr, Jcr;
        
        
        
        const Sub& Sys;
        const Sub& Env;
        const SingleSub& m;
        const SingleSub& n;
        const Parameter& para;
        const Parity& pari;
        QWave _Wave;
        

        
public:

        clock_t tsys, tsm, tsn, tem, ten;

        int rows() { return Dim; };
        int cols() { return Dim; };
        void set_shift(double sigma) { sigma_ = sigma; }
        void perform_op(double *x_in, double *y_out)
        {
                f1tof2(x_in, y_out);
        };
        const QWave& Wave()const{return _Wave;};
        int Dim;
        Super(const Parameter& _para, const Sub& _Sys, const SingleSub& _m, const SingleSub& _n, const Sub& _Env, const Parity& _pari):
        para(_para),
        m(_m),
        n(_n),
        Sys(_Sys),
        Env(_Env),
        _Wave(_Sys, _m, _n, _Env),
        pari(_pari)
        {
               time_t start=clock();
                tsys=tsm=tsn=tem=ten=start-start;
               Dim=_Wave.waveDim(pari); 
        };

        void f1tof2(double* f, double* g)
        {
                vector<double> ff, gg;
                for(int i=0; i<Dim; ++i)
                {
                        ff.push_back(f[i]);
                }
                f1tof2(ff,gg);
                for(int i=0; i<Dim; ++i)
                {
                        g[i]=gg.at(i);
                }
        };

        
        void f1tof2(const vector<double>& f, vector<double>& g)
        {
                _Wave.f2Wave(f, pari);
                QWave temp(Sys, m, n, Env, _Wave, para, pari);
                _Wave=temp;
                _Wave.Wave2f(g, pari);
                tsys+=temp.tsys;
                tsm+=temp.tsm;
                tsn+=temp.tsn;
                tem+=temp.tem;
                ten+=temp.ten;
                //cout<<"haha"<<endl;
        };

};


#endif // SUPER_H
