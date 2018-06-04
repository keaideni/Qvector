#include "Super.h"
void SaveTruncM(const MatrixXd&, const int&);
void Super::f1tof2(const vector<double>& f, vector<double>& g)
{
        _Wave.f2Wave(f);
        OneIteration();
        _Wave.Wave2f(g);
}

void Super::OneIteration()
{
        


        QWave temp(_Wave);//_Wave.Show();int nn; cin>>nn;
        _Wave.SysOPWave(Sys.System());

        /*MatrixXd temp1;
        _Wave.SMEN(temp1);
        cout<<temp1<<endl;
        int nn;cin>>nn;*/
        
        _Wave.EnvOPWave(Env.System(), temp);
        

        
        _Wave.MOPWave(m.System(), temp);
       

        
        _Wave.NOPWave(n.System(), temp);

//==========Sys-m=====================

        _Wave.SysOPWave(Sys.SysA(), temp.MOPWave(m.SysAdag(), -1*Jr));

        _Wave.SysOPWave(Sys.SysAdag(), temp.MOPWave(m.SysA(), -1*Jr));

        _Wave.SysOPWave(Sys.SysAdag(), temp.MOPWave(m.SysAdag(), -1*Jcr));

        _Wave.SysOPWave(Sys.SysA(), temp.MOPWave(m.SysA(), -1*Jcr));
        //cout<<"2"<<endl;

//=================m-Env=====================

        _Wave.EnvOPWave(Env.SysA1(), temp.MOPWave(m.SysAdag(), -1*Jr));

        _Wave.EnvOPWave(Env.SysAdag1(), temp.MOPWave(m.SysA(), -1*Jr));

        _Wave.EnvOPWave(Env.SysAdag1(), temp.MOPWave(m.SysAdag(), -1*Jcr));

        _Wave.EnvOPWave(Env.SysA1(), temp.MOPWave(m.SysA(), -1*Jcr));
        //cout<<"3"<<endl;

//====================Env-n==============================

        _Wave.EnvOPWave(Env.SysA(), temp.NOPWave(n.SysAdag(), -1*Jr));

        _Wave.EnvOPWave(Env.SysAdag(), temp.NOPWave(n.SysA(), -1*Jr));

        _Wave.EnvOPWave(Env.SysAdag(), temp.NOPWave(n.SysAdag(), -1*Jcr));

        _Wave.EnvOPWave(Env.SysA(), temp.NOPWave(n.SysA(), -1*Jcr));
        //cout<<"4"<<endl;

//=============for the periodic bound condition=====================
        //===========n-Sys=====================

        _Wave.SysOPWave(Sys.SysA1(), temp.NOPWave(n.SysAdag(), -1*Jr));

        _Wave.SysOPWave(Sys.SysAdag1(), temp.NOPWave(n.SysA(), -1*Jr));

        _Wave.SysOPWave(Sys.SysAdag1(), temp.NOPWave(n.SysAdag(), -1*Jcr));

        _Wave.SysOPWave(Sys.SysA1(), temp.NOPWave(n.SysA(), -1*Jcr));


//==================================

        //MatrixXd wave;
        //_Wave.SMEN(wave);
        //SaveTruncM(wave, 1000);
        //exit(true);
//==================================

}


void Super::f1tof2(double* f, double* g)
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
}