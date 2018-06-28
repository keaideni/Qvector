#include "Sub.h"




Sub::Sub(const Parameter& para, const int& orbital):
_Orbital(orbital) 
{
       OP a(para, Annihilation), adag(para, Creation), iden(para, Iden);
       OP sigmaz(para, SigmaZ), sigmap(para, SigmaP), sigmam(para, SigmaM);
       OP sigmai(para, SigmaI);

       _SysA.Kron(a, sigmai);//_SysA.show();exit(true);
       _SysAdag.Kron(adag, sigmai);
       _SysEye.Kron(iden, sigmai);
       _SysA1.Kron(a, sigmai);
       _SysAdag1.Kron(adag, sigmai);


       OP SysSigmaZ(iden, sigmaz), SysSigmap(iden, sigmap), SysSigmam(iden, sigmam);

       _System=_SysAdag*_SysA*para.omega0()+SysSigmap*SysSigmam*para.omegaq()
               +para.gr()*(_SysAdag*SysSigmam+SysSigmap*_SysA)
               +para.gcr()*(_SysAdag*SysSigmap+SysSigmam*_SysA);

}


Sub::Sub(const Parameter& para, const Sub& SubL, const Sub& SubR, const int& orbital):
_Orbital(orbital),
_SysA(SubL._SysEye, SubR._SysA),
_SysAdag(SubL._SysEye, SubR._SysAdag),
_SysEye(SubL._SysEye, SubR._SysEye),
_SysA1(SubL._SysA1, SubR._SysEye),
_SysAdag1(SubL._SysAdag1, SubR._SysEye),
_System(SubL._System, SubR._SysEye)
{
        OP temp;
        temp.Kron(SubL._SysEye, SubR._System);
        

        _System.add(temp);

        temp.Kron(SubL._SysAdag, SubR._SysA1);
        _System.add(-1*temp*para.Jr());

        temp.Kron(SubL._SysA, SubR._SysAdag1);
        _System.add(-1*temp*para.Jr());

        temp.Kron(SubL._SysAdag, SubR._SysAdag1);
        _System.add(-1*temp*para.Jcr());

        temp.Kron(SubL._SysA, SubR._SysA1);
        _System.add(-1*temp*para.Jcr());

        

        //temp.Kron(SubL._SysAdag1, SubR._SysA);
        //_System.add(temp*para.Jr());

        //temp.Kron(SubL._SysA1, SubR._SysAdag);
        //_System.add(temp*para.Jr());

        //temp.Kron(SubL._SysAdag1, SubR._SysAdag);
        //_System.add(temp*para.Jcr());

        //temp.Kron(SubL._SysA1, SubR._SysA);
        //_System.add(temp*para.Jcr());

        //_System=temp.Kron(SubL._System, SubR._SysEye);//+temp.Kron(SubL._SysEye, SubR._System)+para.gr()*temp.Kron(SubL._SysAdag, SubR._SysA1)+para.gr()*temp.Kron(SubL._SysA, SubR._SysAdag1)+para.gcr()*temp.Kron(SubL._SysAdag, SubR._SysAdag1)+para.gcr()*temp.Kron(SubL._SysA, SubR._SysA1)+para.gr()*temp.Kron(SubL._SysAdag1, SubR._SysA)+para.gr()*temp.Kron(SubL._SysA1, SubR._SysAdag)+para.gcr()*temp.Kron(SubL._SysAdag1, SubR._SysAdag)+para.gcr()*temp.Kron(SubL._SysA1, SubR._SysA);

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





void Sub::Trunc(const OP& truncU)
{
        
        
        _System.TruncU(truncU);
        
        _SysA.TruncU(truncU);

        _SysAdag.TruncU(truncU);
        _SysEye.TruncU(truncU);
        _SysA1.TruncU(truncU);
        _SysAdag1.TruncU(truncU);
}


void Sub::Save()const
{
        string filename("./block/");
        filename+=itos(_Orbital);

        ofstream outfile(filename);
        outfile<<_Orbital<<endl;
        _System.save(outfile);
        _SysA.save(outfile);
        _SysAdag.save(outfile);
        _SysA1.save(outfile);
        _SysAdag1.save(outfile);
        _SysEye.save(outfile);
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
        _System.read(infile);
        _SysA.read(infile);
        _SysAdag.read(infile);
        _SysA1.read(infile);
        _SysAdag1.read(infile);
        _SysEye.read(infile);

        infile.close();

}



void Sub::Show()const
{
        cout<<"The site of Sub block is "<<_Orbital<<endl;
        //cout<<"The System:"<<endl;
        cout<<"The System:"<<endl; 
        _System.show();
        cout<<"The SysA:"<<endl;
        _SysA.show();
        cout<<"The SysAdag:"<<endl;
        _SysAdag.show();
        cout<<"The SysA1:"<<endl;
        _SysA1.show();
        cout<<"The SysAdag1:"<<endl;
        _SysAdag1.show();
        cout<<"The SysEye:"<<endl;
        _SysEye.show();
}
