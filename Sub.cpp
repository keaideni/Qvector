#include "Sub.h"

string itos(const int& i)
{
        stringstream s;
        s<<i;
        return s.str();
}

void Sub::MatrixSave(const MatrixXd& A, ofstream& outfile)const
{
        outfile<<A.rows()<<endl
        <<A.cols()<<endl;
        outfile.precision(20);
        outfile<<A<<endl;
}

void Sub::MatrixRead(MatrixXd& A, ifstream& infile)
{
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
}


void Sub::Kron(MatrixXd& ab,const MatrixXd& a, const MatrixXd& b)
{
        ab=MatrixXd::Zero(a.rows()*b.rows(),a.cols()*b.cols());

        int sizer(b.rows()),sizec(b.cols());
        for(int i=0; i<a.rows(); ++i)
        {
                for(int j=0; j<a.cols(); ++j)
                {
                        int startr(i*b.rows()),startc(j*b.cols());

                        ab.block(startr,startc,sizer, sizec)=a(i,j)*b;
                }
        }
}





Sub::Sub(const Parameter& para, const int& orbital):
_Orbital(orbital),
_SysEye(MatrixXd::Identity((nmax+1)*2,(nmax+1)*2))
{
        MatrixXd adag=MatrixXd::Zero(nmax+1,nmax+1),
        a(MatrixXd::Zero(nmax+1,nmax+1)),
        eye(MatrixXd::Identity(nmax+1, nmax+1));
        Matrix2d sigmaplu, sigmamin, sigmaeye;
        sigmaplu<<0,0,1,0;sigmamin<<0,1,0,0;sigmaeye<<1,0,0,1;//cout<<sigmaeye<<endl;
        for(int i=0; i<nmax; ++i)
        {
                a(i,i+1)=sqrt(i+1);
                adag(i+1,i)=sqrt(i+1);
        }

        Kron(_SysA, a, sigmaeye);//cout<<_SysA<<endl;
        Kron(_SysAdag, adag, sigmaeye);
        _SysA1=_SysA;_SysAdag1=_SysAdag;

        Kron(_System, adag*a, sigmaeye);
        //test===========================
        //_System*=-1;
        //===============================
        MatrixXd temp;
        Kron(temp, eye, sigmaplu*sigmamin); _System+=temp;
        Kron(temp, adag, sigmamin);temp*=para.gr(); _System+=temp;
        Kron(temp, a, sigmaplu);temp*=para.gr(); _System+=temp;
        Kron(temp, adag, sigmaplu); temp*=para.gcr();_System+=temp;
        Kron(temp, a, sigmamin); temp*=para.gcr(); _System+=temp;
}


Sub::Sub(const Parameter& para, const Sub& SubL, const Sub& SubR, const int& orbital):
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
}