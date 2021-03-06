/*************************************************************************
	> File Name: SuperEnergy.cpp
	> Author: 
	> Mail: 
	> Created Time: 2018 Мог 28, Пү 22:34:04
 ************************************************************************/

#include<iostream>
#include "SuperEnergy.h"
using namespace std;
SuperEnergy::SuperEnergy(Parameter&para,Super& sup, const Parity& pari):
wave(sup.Wave())
{
                
        int a(6);
        if(sup.Dim < 6)a=4;
        SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, a);
        time_t start, end;
        time(&start);
        eigs.init();
        eigs.compute(100000);
        time(&end);
        //cout<<"The calculation takes "<<end-start<<"s."<<endl;
        if (eigs.info() == SUCCESSFUL)
        {

                vector<double> f;double sum(0);
                time(&start);

                VectorXd tempvec(eigs.eigenvectors(1));//Here is the must steps, if use eigs directly, the process will take lots of time;
                for(int i=0; i<tempvec.size(); ++i)
                {
                        f.push_back(tempvec(i));
                }
                time(&end);

		wave.f2Wave(f, pari);
                              
                                        
                          
                para.Energy = eigs.eigenvalues()(0);
		
                //std::cout << eigs.num_iterations() << std::endl;
        }
        //cout<<"The values copy takes "<<end-start<<"s."<<endl;

                
}
SuperEnergy::SuperEnergy(Parameter&para,Super& sup, const Parity& pari, const QWave& initwave):
wave(sup.Wave())
{
        
        std::vector<double> f;
        initwave.Wave2f(f, pari);
        double *pt = new double [sup.Dim];
        for(int i = 0; i < sup.Dim; ++i)pt[i] = f.at(i);
                
        int a(6);
        if(sup.Dim < 6)a=4;
        SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, a);
        time_t start, end;
        time(&start);
        eigs.init(pt);
        eigs.compute(100000);
        time(&end);
        //cout<<"The calculation takes "<<end-start<<"s."<<endl;
        if (eigs.info() == SUCCESSFUL)
        {
		vector<double> f;double sum(0);
                time(&start);

                VectorXd tempvec(eigs.eigenvectors(1));//Here is the must steps, if use eigs directly, the process will take lots of time;
                for(int i=0; i<tempvec.size(); ++i)
                {
                        f.push_back(tempvec(i));
                }
                time(&end);

		wave.f2Wave(f, pari);
                              
                                        
                          
                para.Energy = eigs.eigenvalues()(0);
        }
        //cout<<"The values copy takes "<<end-start<<"s."<<endl;

                
};

