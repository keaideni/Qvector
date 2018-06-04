//!!!!!!!!the three functions must be coherent!!!!!!!!
//wave2f, f2wave, f2wave!!!!!!!!!



#include "QWave.h"
struct Eigstruct
{
        double lamda;
        VectorXd state;
};



bool comp(const Eigstruct& a, const Eigstruct& b);
//============================Reshape===================================
void QWave::SMEN2Wave(const MatrixXd& A)
{

        for(int ie=0; ie<DEnv; ++ie)
        {
                for(int in=0; in<Dn; ++in)
                {
                        for(int is=0; is<DSys; ++is)
                        {
                                for(int im=0; im<Dm; ++im)
                                {
                                        _Wave.at(im).at(in)(is,ie)=A(is*Dm+im, ie*Dn+in);
                                }
                        }
                }
        }
}

void QWave::SMEN(MatrixXd& A)const
{
        A=MatrixXd::Zero(DSys*Dm, DEnv*Dn);
        for(int ie=0; ie<DEnv; ++ie)
        {
                for(int in=0; in<Dn; ++in)
                {
                        for(int is=0; is<DSys; ++is)
                        {
                                for(int im=0; im<Dm; ++im)
                                {
                                        A(is*Dm+im, ie*Dn+in)=_Wave.at(im).at(in)(is,ie);
                                }
                        }
                }
        }
}


void QWave::NSME(MatrixXd& A)const
{
        A=MatrixXd::Zero(Dn*DSys, Dm*DEnv);
        for(int im=0; im<Dm; ++im)
        {
                for(int ie=0; ie<DEnv; ++ie)
                {
                        for(int in=0; in<Dn; ++in)
                        {
                                for(int is=0; is<DSys; ++is)
                                {
                                        A(in*DSys+is, im*DEnv+ie)=_Wave.at(im).at(in)(is,ie);
                                }
                        }
                }
        }
}


void QWave::TruncL(MatrixXd& truncU, const int& D)const
{
        MatrixXd Wave;
        SMEN(Wave);

        vector<Eigstruct> denmat;

        if(Wave.cols()*Wave.rows()>16)
        {
                BDCSVD<MatrixXd> svd(Wave, ComputeFullU);
                for(int i=0; i<svd.singularValues().rows(); ++i)
                {
                        Eigstruct base;
                        base.lamda=svd.singularValues()(i);
                        base.state=svd.matrixU().col(i);

                        denmat.push_back(base);
                }
        }else
        {
                JacobiSVD<MatrixXd> svd(Wave, ComputeFullU);
                for(int i=0; i<svd.singularValues().rows(); ++i)
                {
                        Eigstruct base;
                        base.lamda=svd.singularValues()(i);
                        base.state=svd.matrixU().col(i);

                        denmat.push_back(base);
                }
        }

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.rows());
        int ncol(D<denmat.size()?D:denmat.size());
        truncU=(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncU.col(i)=denmat.at(i).state;
        }



}


//===========================================================================
void QWave::SysOPWave(const MatrixXd& Sys)
{
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        MatrixXd temp(_Wave.at(im).at(in));
                        _Wave.at(im).at(in)=Sys*_Wave.at(im).at(in);
                }
        }
}

void QWave::SysOPWave(const MatrixXd& Sys, const QWave& wave)
{
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        _Wave.at(im).at(in)+=Sys*wave.Wave().at(im).at(in);
                }
        }
}

void QWave::EnvOPWave(const MatrixXd& Env, const QWave& wave)
{
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        _Wave.at(im).at(in)+=wave.Wave().at(im).at(in)*Env.transpose();
                }
        }
}

void QWave::MOPWave(const SpMat& M, const QWave& wave)
{
        for(int in=0; in<Dn; ++in)
        {
                for (int k=0; k<M.outerSize(); ++k)
                for (SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
                {
                
                       _Wave.at(it.row()).at(in)+=it.value()*wave.Wave().at(it.col()).at(in);
                }
                
                
        }
}

const QWave QWave::MOPWave(const SpMat& M, const double& j)const
{
        QWave temp(DSys, Dm, Dn, DEnv);
        for(int in=0; in<Dn; ++in)
        {
                for (int k=0; k<M.outerSize(); ++k)
                for (SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
                {
                
                       temp._Wave.at(it.row()).at(in)+=it.value()*_Wave.at(it.col()).at(in)*j;
                }
                
                
        }
        return temp;

}


void QWave::NOPWave(const SpMat& N, const QWave& wave)
{       
        for(int im=0; im<Dm; ++im)
        {
                for (int k=0; k<N.outerSize(); ++k)
                for (SparseMatrix<double>::InnerIterator it(N,k); it; ++it)
                {
                
                       _Wave.at(im).at(it.row())+=it.value()*wave.Wave().at(im).at(it.col());
                }
                
                
        }
}

const QWave QWave::NOPWave(const SpMat& N, const double& j)const
{
        QWave temp(DSys, Dm, Dn, DEnv);

        for (int k=0; k<N.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(N,k); it; ++it)
        {
                for(int im=0; im<Dm; ++im)
                {
                       temp._Wave.at(im).at(it.row())=it.value()*Wave().at(im).at(it.col())*j;
                }
                
                
        }
        return temp;

}


void QWave::add(const QWave& Wave)
{
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        _Wave.at(im).at(in)+=Wave._Wave.at(im).at(in);
                }
        }
}

//===========================================================================

void QWave::Wave2f(vector<double>& f)const
{
        for(int is=0; is<DSys; ++is)
        {
                for(int im=0; im<Dm; ++im)
                {
                        for(int ie=0; ie<DEnv; ++ie)
                        {
                                for(int in=0; in<Dn; ++in)
                                {
                                        f.push_back((_Wave.at(im).at(in))(is, ie));
                                }
                        }
                }
        }
}


void QWave::f2Wave(const vector<double>& f)
{
        int i(0);//for(int i=0; i<f.size(); ++i)cout<<f.at(i)<<endl;
        for(int is=0; is<DSys; ++is)
        {
                //vector<MatrixXd> tempv;
                for(int im=0; im<Dm; ++im)
                {
                        
                        for(int ie=0; ie<DEnv; ++ie)
                        {
                                for(int in=0; in<Dn; ++in)
                                {
                                        _Wave.at(im).at(in)(is, ie)=f.at(i++);
                                }
                        }
                        
                }
                

        }
        //Show();int nn; cin>>nn;
}


void QWave::f2Wave(const VectorXd& f)
{
        int i(0);
        for(int is=0; is<DSys; ++is)
        {
                
                for(int im=0; im<Dm; ++im)
                {
                        
                        for(int ie=0; ie<DEnv; ++ie)
                        {
                                for(int in=0; in<Dn; ++in)
                                {
                                        _Wave.at(im).at(in)(is, ie)=f(i++);
                                }
                        }
                       
                }
               

        }
}


void QWave::Show()const
{
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        cout<<"<"<<im<<", "<<in<<">:"<<endl<<_Wave.at(im).at(in)<<endl;
                }
        }
}



void QWave::Norm()
{
        double sum(0);
        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        for(int is=0; is<DSys; ++is)
                        {
                                for(int ie=0; ie<DEnv; ++ie)
                                    sum+=pow(_Wave.at(im).at(in)(is, ie),2);
                        }
                }
        }

        for(int im=0; im<Dm; ++im)
        {
                for(int in=0; in<Dn; ++in)
                {
                        for(int is=0; is<DSys; ++is)
                        {
                                for(int ie=0; ie<DEnv; ++ie)
                                    _Wave.at(im).at(in)(is, ie)/=sqrt(sum);
                        }
                }
        }
}