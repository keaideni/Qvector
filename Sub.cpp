#include "Sub.h"

string itos(const int& i)
{
        stringstream s;
        s<<i;
        return s.str();
}



Sub::Sub(const Parameter& para, const int& orbital):
_Orbital(orbital), 
{
       OP a(para, Creation), adag(para, Annihilation), iden(para, Iden);
       OP sigmaz(para, SigmaZ), sigmap(para, SigmaP), sigmam(para, SigmaM);
       OP sigmai(para, SigmaI);

       _SysA.Kron(a, sigmai);
       _SysAdag.Kron(adag, sigmai);
       _SysEye.Kron(iden, sigmai);
       _SysA1.Kron(a, sigmai);
       _SysAdag1.Kron(adag, sigmai);

}


/*Sub::Sub(const Parameter& para, const Sub& SubL, const Sub& SubR, const int& orbital):
_Orbital(orbital)
{
        Kron(_System, SubL.System(), SubR.SysEye());MatrixXd temp;
        Kron(temp, SubL.SysEye(), SubR.System());_System+=temp;
        Kron(temp, SubL.SysAdag(), SubR.SysA1());temp*=para.Jr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysA(), SubR.SysAdag1());temp*=para.Jr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysAdag(), SubR.SysAdag1());temp*=para.Jcr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysA(), SubR.SysA1());temp*=para.Jcr();temp*=-1;_System+=temp;

        Kron(_SysA, SubL.SysEye(), SubR.SysA());
        Kron(_SysAdag, SubL.SysEye(), SubR.SysAdag());
        Kron(_SysA1, SubL.SysA1(), SubR.SysEye());
        Kron(_SysAdag1, SubL.SysAdag1(), SubR.SysEye());

        Kron(_SysEye, SubL.SysEye(), SubR.SysEye());
}



const Sub& Sub::operator=(const Sub& a)
{
        _Orbital=a._Orbital;
        _System=a._System;
        _SysA=a._SysA;
        _SysAdag=a._SysAdag;
        _SysEye=a._SysEye;
        _SysA1=a._SysA1;
        _SysAdag1=a._SysAdag1;
}





void Sub::Trunc(const MatrixXd& truncU)
{
        
        //truncU.col(0)<<endl;
        _System=truncU.adjoint()*_System*truncU;
        
        _SysA=truncU.adjoint()*_SysA*truncU;

        _SysAdag=truncU.adjoint()*_SysAdag*truncU;
        _SysEye=truncU.adjoint()*_SysEye*truncU;
        _SysA1=truncU.adjoint()*_SysA1*truncU;
        _SysAdag1=truncU.adjoint()*_SysAdag1*truncU;
}


void Sub::Save()const
{
        string filename("./block/");
        filename+=itos(_Orbital);

        ofstream outfile(filename);
        outfile<<_Orbital<<endl;
        MatrixSave(_System, outfile);
        MatrixSave(_SysA, outfile);
        MatrixSave(_SysAdag, outfile);
        MatrixSave(_SysA1, outfile);
        MatrixSave(_SysAdag1, outfile);
        MatrixSave(_SysEye, outfile);
        outfile.close();
}



void Sub::Read(const int& orbital)
{
        _Orbital=orbital;
        string filename("./block/");
        filename+=itos(orbital);
        ifstream infile(filename);
        if(!infile)
        {
                cerr<<"error: unable to open the file: "<<filename<<endl;
                exit(true);
        }
        infile>>_Orbital;
        MatrixRead(_System, infile);
        MatrixRead(_SysA, infile);
        MatrixRead(_SysAdag, infile);
        MatrixRead(_SysA1, infile);
        MatrixRead(_SysAdag1, infile);
        MatrixRead(_SysEye, infile);

        infile.close();

}



void Sub::Show()const
{
        cout<<"The site of Sub block is "<<_Orbital<<endl;
        cout<<"The System:"<<endl;
        cout<<_System<<endl;//.rows()<<"X"<<_System.cols()<<endl;
        cout<<"The SysA:"<<endl;
        cout<<_SysA.rows()<<"X"<<_SysA.cols()<<endl;
        cout<<"The SysAdag:"<<endl;
        cout<<_SysAdag.rows()<<"X"<<_SysAdag.cols()<<endl;
        cout<<"The SysA1:"<<endl;
        cout<<_SysA1.rows()<<"X"<<_SysA1.cols()<<endl;
        cout<<"The SysAdag1:"<<endl;
        cout<<_SysAdag1.rows()<<"X"<<_SysAdag1.cols()<<endl;
        cout<<"The SysEye:"<<endl;
        cout<<_SysEye.rows()<<"X"<<_SysEye.cols()<<endl;
}*/
