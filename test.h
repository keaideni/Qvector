#ifndef TEST_H
#define TEST_H
#include "QWave.h"
#include "SingleSub.h"
#include "DMRG.h"
#include <ctime>
const MatrixXd TruncL(const MatrixXd& _Wave, const int& D);
const MatrixXd TruncR(const MatrixXd& _Wave, const int& D);

void SaveTruncM(const MatrixXd& A, const int& logo);
void ReadTruncM(MatrixXd& A, const int& logo);
void test(Parameter& para);

void test(Parameter& para)
{
        
        /*SingleSub m(para);
        Sub Sys(para, 1);
        //MatrixXd System;
        MatrixXd IniWave;

        //ReadTruncM(System, 1002);
        ReadTruncM(IniWave, 1003);
        
        

        Super Sup(para, Sys, Sys);

        vector<double> f;
        for(int i=0; i<100; ++i)
        {
            for(int j=0; j<100; ++j)f.push_back(IniWave(i,j));
        }
        Sup._Wave.f2Wave(f);
        
        Sup.OneIteration();
        Sup.Wave().SMEN(IniWave);
        
        SaveTruncM(IniWave, 1000);*/
        
        MatrixXd A, B;
        ReadTruncM(A, 1000);ReadTruncM(B, 1001);
        double sum(0);
        for(int i=0; i<B.rows(); ++i)
        {
            for(int j=0; j<B.cols(); ++j)
            {
                cout<<A(i,j)<<"\t"<<B(i,j)<<"\t"<<
                abs(A(i,j)-B(i,j))<<endl;
                sum+=abs(A(i,j)-B(i,j));
            }
        }cout<<sum<<endl;//<<A.transpose()*B<<endl;;

        
        

}

#endif // TEST_H