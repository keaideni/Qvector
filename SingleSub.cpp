#include "SingleSub.h"
#include <cmath>

void SingleSub::Kron(SpMat& ab, const SpMat& a, const SpMat& b)
{
        ab.setZero();
        ab.resize(a.rows()*b.rows(), a.cols()*b.cols());

        for (int k=0; k<a.outerSize(); ++k)
                for(SpMat::InnerIterator it(a,k); it; ++it)
                {
                        for(int l=0; l<b.outerSize(); ++l)
                        {
                                for(SpMat::InnerIterator itt(b,l); itt; ++itt)
                                {
                                        ab.insert(it.row()*b.rows()+itt.row(), 
                                                it.col()*b.cols()+itt.col())
                                        =it.value()*itt.value();
                                }
                        }
                }

}

SingleSub::SingleSub(const Parameter& para):
_System(para.nmax()*2, para.nmax()*2),
_SysA(para.nmax()*2, para.nmax()*2),
_SysAdag(para.nmax()*2, para.nmax()*2)
{
        SpMat tempA(para.nmax()+1, para.nmax()+1), 
        tempAdag(para.nmax()+1, para.nmax()+1),
        tempEye(para.nmax()+1, para.nmax()+1);
        for(int i=0; i<para.nmax(); ++i)
        {
                tempA.insert(i, i+1)=sqrt(i+1);
                tempAdag.insert(i+1, i)=sqrt(i+1);
                tempEye.insert(i, i)=1;
        }//cout<<tempA<<endl;
        tempEye.insert(para.nmax(), para.nmax())=1;
        SpMat Sigmamin(2,2), Sigmaplu(2,2), Sigmaeye(2,2);
        Sigmaeye.insert(0,0)=1; Sigmaeye.insert(1,1)=1;
        Sigmamin.insert(0,1)=1; Sigmaplu.insert(1,0)=1;

        Kron(_SysA, tempA, Sigmaeye);Kron(_SysAdag, tempAdag, Sigmaeye);

        Kron(_System, tempAdag*tempA, Sigmaeye);
        //test=================================
        //_System*=-1;
        //=====================================

        SpMat temp;
        Kron(temp, tempEye, Sigmaplu*Sigmamin); _System+=temp;
        Kron(temp, tempAdag, Sigmamin); temp*=para.gr(); _System+=temp;
        Kron(temp, tempA, Sigmaplu); temp*=para.gr(); _System+=temp;

        Kron(temp, tempAdag, Sigmaplu); temp*=para.gcr(); _System+=temp;
        Kron(temp, tempA, Sigmamin); temp*=para.gcr(); _System+=temp;

}
