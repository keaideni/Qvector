//The Initial wave must be coherent with the QWave order!!!!.

#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>
#include <iostream>
#include "Super.h"

using namespace Spectra;

#ifndef SUPERENERGY_H
#define SUPERENERGY_H
class SuperEnergy
{
private:
	//double _excited;
public:
        QWave wave; 
	//const double excited()const{
		//return _excited;
	//};


        //SuperEnergy(){};
         
        SuperEnergy(Parameter& para, Super& sup, const Parity& pari);
        SuperEnergy(Parameter& para, Super& sup, const Parity& pari, const QWave& InitWave);

        
        
	
/*	SuperEnergy(Parameter&para,Super& sup, const MatrixXd& initwave, double& excited):
        wave(sup.Wave())
        {
                
                std::vector<double> f;
                //wave=initwave;
                int Dm(wave.Wave().size());auto it=wave.Wave().begin();
                int Dn(it->size()); auto itt=it->begin();
                int DSys(itt->rows()), DEnv(itt->cols());
                for(int is=0; is<DSys; ++is)
                {
                        for(int im=0; im<Dm; ++im)
                        {
                                for(int ie=0; ie<DEnv; ++ie)
                                {
                                        for(int in=0; in<Dn; ++in)
                                        f.push_back(initwave(is*Dm+im, ie*Dn+in));
                                }
                        }
                }
                double *pt = new double [sup.Dim];
                for(int i = 0; i < sup.Dim; ++i)pt[i] = f.at(i);
                
                int a(8);
                if(sup.Dim < 6)a=4;
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 2, a);
                eigs.init(pt);
                eigs.compute(10000);
                if (eigs.info() == SUCCESSFUL)
                {
			//if(eigs.eigenvalues()(0)<eigs.eigenvalues()(1))
			//{

				//wave.f2Wave(eigs.eigenvectors(1).col(0));
			//}else
			//{
			//	wave.f2Wave(eigs.eigenvectors(2).col(1));
			//}

                        //para.Energy = eigs.eigenvalues()(0);
                        //para.Energy = eigs.eigenvalues()(0);
                        excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
                }

                        std::cout << eigs.num_iterations() << std::endl;

                
        };
        
	
	SuperEnergy(Parameter&para,Super& sup, const QWave& initwave):
        wave(sup.Wave())
        {
                
                std::vector<double> f;
                wave=initwave;
                wave.Wave2f(f);
                double *pt = new double [sup.Dim];
                for(int i = 0; i < sup.Dim; ++i)pt[i] = f.at(i);
                
                int a(6);
                if(sup.Dim < 6)a=4;
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, a);
                eigs.init(pt);
                eigs.compute(10000);
                if (eigs.info() == SUCCESSFUL)
                {
			//if(eigs.eigenvalues()(0)<eigs.eigenvalues()(1))
			//{

				wave.f2Wave(eigs.eigenvectors(1).col(0));
			//}else
			//{
				//wave.f2Wave(eigs.eigenvectors(2).col(1));
			//}
                        //para.Energy = eigs.eigenvalues()(0);
                        //std::cout << eigs.num_iterations() << std::endl;
                        para.Energy = eigs.eigenvalues()(0);
//                        _excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
                }

                
        };*/
        
};

#endif
SuperEnergy::SuperEnergy(Parameter&para,Super& sup, const Parity& pari):
        wave(sup.Wave())
        {
                
                int a(6);
                if(sup.Dim < 6)a=4;
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, a);
                time_t start, end;
                time(&start);
                eigs.init();
                eigs.compute(10000);
                time(&end);
                cout<<"The calculation takes "<<end-start<<"s."<<endl;
                if (eigs.info() == SUCCESSFUL)
                {
			//if(eigs.eigenvalues()(0)<eigs.eigenvalues()(1))
			//{

                        vector<double> f;double sum(0);
                        time(&start);

                        VectorXd tempvec(eigs.eigenvectors(1));//Here is the must steps, if use eigs directly, the process will take lots of time;
                        for(int i=0; i<tempvec.size(); ++i)
                        {
                                f.push_back(tempvec(i));
                        }
                        time(&end);

				wave.f2Wave(f, pari);
                              
                                        
                          
			//}else
			//{
				//wave.f2Wave(eigs.eigenvectors(2).col(1));
			//}
                        para.Energy = eigs.eigenvalues()(0);
                        //_excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
		
                        //std::cout << eigs.num_iterations() << std::endl;
                }
                cout<<"The values copy takes "<<end-start<<"s."<<endl;

                
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
                eigs.init(pt);
                time_t start, end;
                time(&start);
                eigs.init();
                eigs.compute(10000);
                time(&end);
                cout<<"The calculation takes "<<end-start<<"s."<<endl;
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
                              
                                        
                          
			//}else
			//{
				//wave.f2Wave(eigs.eigenvectors(2).col(1));
			//}
                        para.Energy = eigs.eigenvalues()(0);
                        //_excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
                }
                cout<<"The values copy takes "<<end-start<<"s."<<endl;

                
        };
        
