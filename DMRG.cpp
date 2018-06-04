#include "DMRG.h"
#include "Trunc.h"

#include <iomanip>
#include <ctime>




void SaveTruncM(const MatrixXd& A, const int& logo);
void SaveTruncM(const MatrixXd& A, const int& logo)
{
        string filename("./trunc/");
        filename+=itos(logo);

        ofstream outfile(filename);
        outfile<<A.rows()<<endl<<A.cols()<<endl;
        outfile.precision(20);
        outfile<<A<<endl;
        
        outfile.close();
}


void ReadTruncM(MatrixXd& A, const int& logo);
void ReadTruncM(MatrixXd& A, const int& logo)
{
        string filename("./trunc/");
        filename+=itos(logo);

        ifstream infile(filename);
        if(!infile)
        {
                cerr<<"Error! The file "<<filename<<" doesn't exist!!"<<endl;
                exit(true);
        }
        int rows, cols;
        infile>>rows>>cols;
        A=MatrixXd::Zero(rows, cols);

        for(int i=0; i<rows; ++i)
        {
                for(int j=0; j<cols; ++j)
                {
                        infile>>A(i,j);
                }
        }
        infile.close();

        
}





DMRG::DMRG(Parameter& para):
Sys(para, 1),
Env(para, para.LatticeSize()),
m(para, 2),
n(para, para.LatticeSize()-1),
SaveAll("./result/SaveAll"),
_FEnergy(100000),
_Entropy(0)
{
        SaveAll.precision(15);
        Sys.Save(); Env.Save();

        int OS(1), OE(para.LatticeSize());

        //cout<<"===================The growth process:==================="<<endl;
        SaveAll<<"===================The growth process:==================="<<endl;
        Parameter paraup;
        paraup.ChangeD(20);
        BuildUp(paraup, OS, OE);
        //cout<<"===================The sweep process:==================="<<endl;
        SaveAll<<"===================The sweep process:==================="<<endl;


        //para.ChangeD(200);
        OS-=1;OE+=1;//This one for the IniWave works.
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

                //Sub SysNew(para, Sys, m, OS+1);
                
                

                Super Sup(para, Sys, Env);
                SuperEnergy Supp(para, Sup);


                cout.precision(15);
                //cout<<"OS="<<setw(10)<<OS<<"; OE="<<setw(10)<<OE<<"; The energy="
                //<<setw(18)<<para.Energy<<endl;
                SaveAll<<"OS="<<setw(10)<<OS<<"; OE="<<setw(10)<<OE<<"; The energy="
                <<setw(18)<<para.Energy<<endl;
                //int nn; cin>>nn;

                MatrixXd wave;

                Supp.wave.SMEN(wave);
                
                
                MatrixXd MatrixU(DenTruncL(wave, para.D()));
                MatrixXd MatrixV(DenTruncR(wave, para.D()));


                Sub SysNew(para, Sys, m, OS+1);
                SysNew.Trunc(MatrixU);
                SysNew.Save();
                SaveTruncM(MatrixU, SysNew.Orbital());

                Sub EnvNew(para, Env, n, OE-1);
                EnvNew.Trunc(MatrixV);
                EnvNew.Save();
                SaveTruncM(MatrixV, Env.Orbital());

                m.ChangeOrbital(m.Orbital()+1);
                n.ChangeOrbital(n.Orbital()-1);
                OS+=1;
                OE-=1;

                if(OE-OS==1)
                {
                        Supp.wave.SMEN(IniWave);
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
                Sys.Read(OS);Env.Read(OE);m.ChangeOrbital(OS+dir); n.ChangeOrbital(OE+dir);
                
                CalcuEnergy(para, OS, OE, dir, Gdir);
                


                Initialize(dir, Gdir, OS, OE);

                OS+=dir;
                OE+=dir;


                if((OS==(para.LatticeSize()/2)&dir==1)|(OE==(para.LatticeSize()/2+1)&dir==-1))
                {
                        err=abs(para.Energy-_FEnergy);SaveAll<<err<<endl;
                        _FEnergy=para.Energy;
                        Gdir*=-1;
                        stop=(err<1e-6);

                        if(!stop)
                        SaveAll<<"==========the "<<SweepNo++<<"th sweeps=============="<<endl;
                        //cout<<"==========the "<<SweepNo++<<"th sweeps=============="<<endl;
                        if(SweepNo==2)Gdir=-1;
                        
                }
                if(OE==para.LatticeSize()|OS==1)
                {
                        dir*=-1;Gdir*=-1;
                }
                




        }

}




void DMRG::CalcuEnergy(Parameter& para, int& OS, int& OE, const int& dir, const int& Gdir)
{
        
        

        time_t start, end;
        
        
        
        Super Sup(para, Sys, Env);
        
        //cout<<"The process of constructing Super takes "<<(end-start)<<"s."<<endl;
        //SaveAll<<"The process of constructing Super takes "<<(end-start)<<"s."<<endl;
        

        //time(&start);
        SuperEnergy Supp(para, Sup, IniWave);

	//_Excited=Supp.excited();
        //time(&end);
        //cout<<"The process of getting eigenstate takes "<<(end-start)<<"s."<<endl;
        //==============to save the final wave=======================
        if((OS==(para.LatticeSize()/2-1)&dir==1)|(OE==(para.LatticeSize()/2+2)&dir==-1))
        {
                MatrixXd finalwave;
                Supp.wave.SMEN(finalwave);
                SaveTruncM(finalwave, 10000);//The file "/Trunc/10000" is the final wave.
		SuperEnergy Suppp(para, Sup, IniWave, _Excited);
        }



        //cout.precision(15);
        //cout<<"OS="<<setw(10)<<OS<<"; OE="<<setw(10)<<OE<<"; The energy="
        //<<setw(18)<<para.Energy<<endl;
        SaveAll<<"OS="<<setw(10)<<OS<<"; OE="<<setw(10)<<OE<<"; The energy="
        <<setw(18)<<para.Energy<<endl;

        MatrixXd wave;

        if(Gdir==1)
        {
                if(dir==1)
                {
                        
                        Sub SysNew(para, Sys, m, OS+dir);
                        
                        Supp.wave.SMEN(wave);
                        MatrixXd matrixT=(DenTruncL(wave, para.D(), _Entropy));
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        SaveTruncM(matrixT, SysNew.Orbital());
                        IniWave=matrixT.adjoint()*wave;
                }else
                {
                        
                        Sub SysNew(para, Env, n, OE+dir);
                        
                        Supp.wave.SMEN(wave);
                        MatrixXd matrixT=(DenTruncR(wave, para.D(), _Entropy));
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        SaveTruncM(matrixT, SysNew.Orbital());
                        IniWave=wave*matrixT;
                }
        }else
        {
                
                if(dir==1)
                {
                        
                        Sub SysNew(para, n, Sys, OS+dir);
                        
                        Supp.wave.NSME(wave);
                        MatrixXd matrixT=(DenTruncL(wave, para.D(), _Entropy));
                        //cout<<matrixT.rows()<<"x"<<matrixT.cols()<<endl;
                        //int nn; cin>>nn;
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        SaveTruncM(matrixT, SysNew.Orbital());
                        IniWave=matrixT.adjoint()*wave;
                }else
                {
                        
                        Sub SysNew(para, m, Env, OE+dir);
                        
                        Supp.wave.NSME(wave);
                        MatrixXd matrixT=(DenTruncR(wave, para.D(), _Entropy));
                        SysNew.Trunc(matrixT);
                        SysNew.Save();
                        SaveTruncM(matrixT, SysNew.Orbital());
                        IniWave=wave*matrixT;
                }
        }

}






void DMRG::Initialize(const int& dir, const int& Gdir, const int& OS, const int& OE)
{

        int Dm(m.SysEye().cols()), Dn(n.SysEye().cols());


        if(dir==1)
        {

                int De(Env.SysEye().cols()), Ds, Dee;
                {
                        Sub temp;
                        temp.Read(OE+1);
                        Dee=temp.SysEye().cols();
                        temp.Read(OS+1);
                        Ds=temp.SysEye().cols();
                }

                if(Gdir==1)
                {


                        MatrixXd temp(MatrixXd::Zero(Dn*Ds, De));
                        for(int is=0; is<IniWave.rows(); ++is)
                        {
                                for(int in=0; in<Dn;++in)
                                {
                                        for(int ie=0; ie<De; ++ie)
                                        {
                                                temp(in*Ds+is, ie)=IniWave(is, ie*Dn+in);
                                        }
                                }
                        }
                        MatrixXd tempV;
                        ReadTruncM(tempV, OE);
                        temp=temp*tempV.adjoint();


                        IniWave=MatrixXd::Zero(Ds*Dm, Dee*Dn);
                        for(int is=0; is<Ds; ++is)
                        {
                                for(int ie=0; ie<Dee; ++ie)
                                {
                                        for(int im=0; im<Dm; ++im)
                                        {
                                                for(int in=0; in<Dn; ++in)
                                                IniWave(is*Dm+im, ie*Dn+in)
                                                =temp(in*Ds+is, im*Dee+ie);
                                        }
                                }
                        }


                }else
                {
                        MatrixXd temp(MatrixXd::Zero(Ds*Dm, De));
                        for(int is=0; is<Ds; ++is)
                        {
                                for(int im=0; im<Dm;++im)
                                {
                                        for(int ie=0; ie<De; ++ie)
                                        {
                                                temp(is*Dm+im, ie)=IniWave(is, im*De+ie);
                                        }
                                }
                        }
                        MatrixXd tempV;
                        ReadTruncM(tempV, OE);
                        temp=temp*tempV.adjoint();


                        IniWave=MatrixXd::Zero(Dn*Ds, Dm*Dee);
                        for(int is=0; is<Ds; ++is)
                        {
                                for(int ie=0; ie<Dee; ++ie)
                                {
                                        for(int im=0; im<Dm; ++im)
                                        {
                                                for(int in=0; in<Dn; ++in)
                                                IniWave(is*Dm+im, ie*Dn+in)
                                                =temp(is*Dm+im, ie*Dn+in);
                                        }
                                }
                        }
                }
        }else
        {


                int De, Ds(Sys.SysEye().cols()), Dss;
                {
                        Sub temp;
                        temp.Read(OS-1);
                        Dss=temp.SysEye().cols();

                        temp.Read(OE-1);
                        De=temp.SysEye().cols();
                }
                if(Gdir==1)
                {


                        MatrixXd temp(MatrixXd::Zero(Ds, Dm*De));
                        for(int is=0; is<Ds; ++is)
                        {
                                for(int im=0; im<Dn;++im)
                                {
                                        for(int ie=0; ie<De; ++ie)
                                        {
                                                temp(is, im*De+ie)=IniWave(is*Dm+im, ie);
                                        }
                                }
                        }
                        MatrixXd tempU;
                        ReadTruncM(tempU, OS);
                        temp=tempU*temp;


                        IniWave=MatrixXd::Zero(Dss*Dm, De*Dn);
                        for(int is=0; is<Dss; ++is)
                        {
                                for(int ie=0; ie<De; ++ie)
                                {
                                        for(int im=0; im<Dm; ++im)
                                        {
                                                for(int in=0; in<Dn; ++in)
                                                IniWave(is*Dm+im, ie*Dn+in)
                                                =temp(in*Dss+is, im*De+ie);
                                        }
                                }
                        }


                }else
                {
                        MatrixXd temp(MatrixXd::Zero(Ds, De*Dn));
                        for(int is=0; is<Ds; ++is)
                        {
                                for(int in=0; in<Dn;++in)
                                {
                                        for(int ie=0; ie<De; ++ie)
                                        {
                                                temp(is, ie*Dn+in)=IniWave(in*Ds+is, ie);
                                        }
                                }
                        }
                        MatrixXd tempU;
                        ReadTruncM(tempU, OS);
                        temp=tempU*temp;


                        IniWave=MatrixXd::Zero(Dn*Dss, Dm*De);
                        for(int is=0; is<Dss; ++is)
                        {
                                for(int ie=0; ie<De; ++ie)
                                {
                                        for(int im=0; im<Dm; ++im)
                                        {
                                                for(int in=0; in<Dn; ++in)
                                                IniWave(is*Dm+im, ie*Dn+in)
                                                =temp(is*Dm+im, ie*Dn+in);
                                        }
                                }
                        }
                }
        }


}


