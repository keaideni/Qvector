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
        SuperEnergy(Parameter&para,Super& sup):
        wave(sup.Wave())
        {
                
                int a(6);
                if(sup.Dim < 6)a=4;
                SymEigsSolver<double, SMALLEST_ALGE, Super> eigs(&sup, 1, a);
                eigs.init();
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
                        para.Energy = eigs.eigenvalues()(0);
                        //_excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
		
                        //std::cout << eigs.num_iterations() << std::endl;
                }

                
        };
        

        
        SuperEnergy(Parameter&para,Super& sup, const MatrixXd& initwave):
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
			//	wave.f2Wave(eigs.eigenvectors(2).col(1));
			//}

                        //para.Energy = eigs.eigenvalues()(0);
                        //std::cout << eigs.num_iterations() << std::endl;
                        para.Energy = eigs.eigenvalues()(0);
                        //_excited= eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
                }

                
        };
        
	
	SuperEnergy(Parameter&para,Super& sup, const MatrixXd& initwave, double& excited):
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

                
        };
        
};

#endif
