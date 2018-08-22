#include "DMRG.h"
//#include "Trunc.h"

#include <iomanip>
#include <ctime>

string ptos(const Parity& _pari);
string ptos(const Parity& _pari)
{
        if(_pari==Positive)
        return "p";
        else return "n";
}


DMRG::DMRG(Parameter& para, const Parity& _pari):
Sys(para, 1),
Env(para, para.LatticeSize()),
m(para, 2),
n(para, para.LatticeSize()-1),
SaveAll("./result/SaveAll"+ptos(_pari)),
_FEnergy(100000),
_Entropy(0),
pari(_pari)
{
        SaveAll.precision(15);
        Sys.Save(); Env.Save();

        int OS(1), OE(para.LatticeSize());

        cout<<"===================The growth process:==================="<<endl;
        SaveAll<<"===================The growth process:==================="<<endl;
        Parameter paraup;
        paraup.ChangeD(40);
        BuildUp(paraup, OS, OE);
        cout<<"===================The sweep process:==================="<<endl;
        SaveAll<<"===================The sweep process:==================="<<endl;


        //para.ChangeD(200);
        OS-=1;OE+=1;//This one for the IniWave works.
        //OS=4; OE=7;
        Sys.Read(OS); Env.Read(OE);
        Sweep(para, OS, OE);

        //cout<<"===================The sweep process finished!================="<<endl;
        SaveAll<<"===================The sweep process finished!================="<<endl;





        //OneSiteSweep(para, OS, OE);
        SaveAll.close();




};





void DMRG::BuildUp(Parameter& para, int& OS, int& OE)
{
        while(OE-OS>1)
        {

                Sys.Read(OS);Env.Read(OE);//if(OS==2)Sys.Show();
                //Env.read(OE);
                //Sys.SysEye().show();
                //Sub SysNew(para, Sys, m, OS+1);
                
                
                time_t start, end;
                time(&start);
                SingleSub haha(para); 
                Super Sup(para, Sys, haha, haha, Env, pari);
                SuperEnergy Supp(para, Sup, pari);
                time(&end);
                //cout<<"The process of getting eigenstate takes "<<(end-start)<<"s."<<endl;
                //cout<<"the system operator takes "<<(double)(Sup.tsys)/CLOCKS_PER_SEC<<"s."<<endl;        
                //cout<<"the sm operator takes "<<(double)(Sup.tsm)/CLOCKS_PER_SEC<<"s."<<endl;        
                //cout<<"the sn operator takes "<<(double)(Sup.tsn)/CLOCKS_PER_SEC<<"s."<<endl;        
                //cout<<"the em operator takes "<<(double)(Sup.tem)/CLOCKS_PER_SEC<<"s."<<endl;        
                //cout<<"the en operator takes "<<(double)(Sup.ten)/CLOCKS_PER_SEC<<"s."<<endl;        
                

                cout.precision(15);
                cout<<"OS=\t"<<setw(10)<<OS<<"; OE=\t"<<setw(10)<<OE<<"; The energy=\t"
                <<setw(18)<<para.Energy<<endl;
                SaveAll<<"OS=\t"<<setw(10)<<OS<<"; OE=\t"<<setw(10)<<OE<<"; The energy=\t"
                <<setw(18)<<para.Energy<<endl;
                //int nn; cin>>nn;

                OP OPwave;
                
                time(&start);
                Supp.wave.Wave2SMEN(OPwave, pari);
                
                OP MatrixU;
                double err;
                MatrixU.DenTruncL(OPwave, para.D(), err);
                OP MatrixV;
                MatrixV.DenTruncR(OPwave, para.D(), err);


                Sub SysNew(para, Sys, m, OS+1);
                SysNew.Trunc(MatrixU);
                SysNew.Save();
                MatrixU.TruncSave(SysNew.Orbital());

                Sub EnvNew(para, Env, n, OE-1);
                EnvNew.Trunc(MatrixV);

                //EnvNew.SysEye().show();
                //MatrixV.show();

                EnvNew.Save();
                MatrixV.TruncSave(EnvNew.Orbital());
                time(&end);
                //cout<<"The process of truncate takes "<<(end-start)<<"s."<<endl;
                m.ChangeOrbital(m.Orbital()+1);
                n.ChangeOrbital(n.Orbital()-1);
                OS+=1;
                OE-=1;

                if(OE-OS==1)
                {
                        InitWave=Supp.wave;
                }

        }


}



void DMRG::Sweep(Parameter& para, int& OS, int& OE)
{
        int dir(1);//for the Sub lattice growth. 1 means the System grows and -1 means the Environment.
        int Gdir(-1);//for the grow direction, -1 means left, 1 means right.

        //double menergy(0);
        double err(1);
        int SweepNo(1);
        bool stop(false);
        

        while(!stop)
        {
                m.ChangeOrbital(OS+dir); n.ChangeOrbital(OE+dir);
                if(OE==para.LatticeSize()|OS==1)
                {
                        dir*=-1;Gdir*=-1;
                }
                
                CalcuEnergy(para, OS, OE, dir, Gdir);
                


                //cout<<dir<<endl;
                Initialize(para, dir, Gdir, OS, OE);
                //cout<<dir<<endl;

                if((OS==(para.LatticeSize()/2)&dir==1)|(OE==(para.LatticeSize()/2+1)&dir==-1))
                {
                        err=abs(para.Energy-_FEnergy);SaveAll<<err<<endl;
                        _FEnergy=para.Energy;
                        Gdir*=-1;
                        stop=(err<1e-6);

                        if(!stop)
                        cout<<"==========the "<<SweepNo<<"th sweeps=============="<<endl;
                        SaveAll<<"==========the "<<SweepNo++<<"th sweeps=============="<<endl;
                        if(SweepNo==2)Gdir=-1;
                        
                }
                




        }

}




void DMRG::CalcuEnergy(Parameter& para, int& OS, int& OE, const int& dir, const int& Gdir)
{
        
        

        time_t start, end;
        
        time(&start); 
        SingleSub haha(para); 
        Super Sup(para, Sys, haha, haha, Env, pari);
        time(&end);
        //cout<<"The process of constructing Super takes "<<(end-start)<<"s."<<endl;
        //SaveAll<<"The process of constructing Super takes "<<(end-start)<<"s."<<endl;
        
        time(&start);
        SuperEnergy Supp(para, Sup, pari, InitWave);
       
        //==============to save the final wave=======================
        if((OS==(para.LatticeSize()/2-1)&dir==1)|(OE==(para.LatticeSize()/2+2)&dir==-1))
        {
                OP finalwave;
                Supp.wave.Wave2SMEN(finalwave, pari);
                //if(pari==Positive)
                finalwave.TruncSave(10000);//The file "/Trunc/10000" is the final wave.
                //else
                //finalwave.TruncSave(10001);
        }
         //SuperEnergy Supp(para, Sup, pari, InitWave);

	//_Excited=Supp.excited();
        time(&end);
        //cout<<"The process of getting eigenstate takes "<<(end-start)<<"s."<<endl;
        


        cout.precision(15);
        cout<<"OS=\t"<<setw(10)<<OS<<"; OE=\t"<<setw(10)<<OE<<"; The energy=\t"
        <<setw(18)<<para.Energy<<endl;
        SaveAll<<"OS=\t"<<setw(10)<<OS<<"; OE=\t"<<setw(10)<<OE<<"; The energy=\t"
        <<setw(18)<<para.Energy<<endl;

        OP wave;

        if(Gdir==1)
        {
                if(dir==1)
                {
                        
                        Sub SysNew(para, Sys, m, OS+dir);
                        
                        Supp.wave.Wave2SMEN(wave, pari);
                        OP matrixT;
                        matrixT.DenTruncL(wave, para.D(), _Entropy);
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        matrixT.TruncSave(SysNew.Orbital());
                        InitOPWave.LWavetime(matrixT, wave);
                        //cout<<"haha"<<endl;
                }else
                {
                        
                        Sub SysNew(para, Env, n, OE+dir);
                        
                        Supp.wave.Wave2SMEN(wave, pari);
                        OP matrixT;
                        matrixT.DenTruncR(wave, para.D(), _Entropy);
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        matrixT.TruncSave(SysNew.Orbital());
                        InitOPWave.RWavetime(wave, matrixT);
                }
        }else
        {
                
                if(dir==1)
                {
                        
                        Sub SysNew(para, n, Sys, OS+dir);
                        
                        Supp.wave.Wave2NSME(wave, pari);
                        OP matrixT;
                        matrixT.DenTruncL(wave, para.D(), _Entropy);
                        //cout<<matrixT.rows()<<"x"<<matrixT.cols()<<endl;
                        //int nn; cin>>nn;
                        SysNew.Trunc(matrixT);
                        //SysNew.SysEye().show();
                        SysNew.Save();
                        matrixT.TruncSave(SysNew.Orbital());
                        //matrixT.show();
                        //wave.show();
                        InitOPWave.LWavetime(matrixT, wave);
                }else
                {
                        
                        Sub SysNew(para, m, Env, OE+dir);
                        
                        Supp.wave.Wave2NSME(wave, pari);
                        OP matrixT;
                        matrixT.DenTruncR(wave, para.D(), _Entropy);
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        matrixT.TruncSave(SysNew.Orbital());
                        InitOPWave.RWavetime(wave, matrixT);
                }
        }

        //exit(true);

}






void DMRG::Initialize(const Parameter& para, const int& dir, const int& Gdir, int& OS, int& OE)
{
        SingleSub haha(para);
        OP tempOPWave;
        
        if(Gdir==1)
        {
                if(dir==1)
                {

                        OS+=dir;
                        //Sys.SysEye().show();
                        Sys.Read(OS);
                        QWave temp(Sys, haha, haha, Env);
                        //InitOPWave.show();
                        //Env.SysEye().show();
                        //Sys.SysEye().show();
                        //exit(true);
                        temp.d1g1(tempOPWave, InitOPWave, pari);
                        OP tempOPE;
                        tempOPE.TruncRead(OE);
                        InitOPWave.RWavetime2(tempOPWave, tempOPE);
                        OE+=dir;
                        Env.Read(OE);

                        QWave NewInitWave(Sys, haha, haha, Env);
                        NewInitWave.NSME2Wave(InitOPWave, pari);
                        //cout<<"haha"<<endl;
                        InitWave=NewInitWave;
                        //exit(true);
                }else
                {
                        OE+=dir;
                        //Sys.SysEye().show();
                        Env.Read(OE);
                        QWave temp(Sys, haha, haha, Env);
                        //InitOPWave.show();
                        //Env.SysEye().show();
                        //Sys.SysEye().show();
                        //exit(true);
                        temp.dmg1(tempOPWave, InitOPWave, pari);
                        OP tempOPS;
                        tempOPS.TruncRead(OS);
                        InitOPWave.LWavetime2(tempOPS, tempOPWave);

                        OS+=dir;
                        Sys.Read(OS);
                        //exit(true);
                        

                        QWave NewInitWave(Sys, haha, haha, Env);
                        NewInitWave.NSME2Wave(InitOPWave, pari);
                        InitWave=NewInitWave;

                }
                
        }else
        {
                if(dir==1)
                {
                        OS+=dir;
                        //Sys.SysEye().show();
                        Sys.Read(OS);
                        QWave temp(Sys, haha, haha, Env);
                        //InitOPWave.show();
                        //Env.SysEye().show();
                        //Sys.SysEye().show();
                        //exit(true);
                        temp.d1gm(tempOPWave, InitOPWave, pari);
                        OP tempOPE;
                        //cout<<"haha"<<OE<<endl;
                        tempOPE.TruncRead(OE);
                        //tempOPWave.show();
                        //tempOPE.show();
                        InitOPWave.RWavetime2(tempOPWave, tempOPE);
                        
                        OE+=dir;
                        Env.Read(OE);

                        QWave NewInitWave(Sys, haha, haha, Env);
                        NewInitWave.SMEN2Wave(InitOPWave, pari);
                        InitWave=NewInitWave;
                        //exit(true);
                }else
                {
                        OE+=dir;
                        //Sys.SysEye().show();
                        Env.Read(OE);
                        QWave temp(Sys, haha, haha, Env);
                        //InitOPWave.show();
                        //Env.SysEye().show();
                        //Sys.SysEye().show();
                        //exit(true);
                        temp.dmgm(tempOPWave, InitOPWave, pari);
                        
                        OP tempOPS;
                        tempOPS.TruncRead(OS);
                        InitOPWave.LWavetime2(tempOPS, tempOPWave);

                        OS+=dir;
                        Sys.Read(OS);
                        //exit(true);

                        QWave NewInitWave(Sys, haha, haha, Env);
                        NewInitWave.SMEN2Wave(InitOPWave, pari);
                        InitWave=NewInitWave;

                }
        }
}


