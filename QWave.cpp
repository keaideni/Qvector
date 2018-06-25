//!!!!!!!!the three functions must be coherent!!!!!!!!
//wave2f, f2wave, f2wave!!!!!!!!!



#include "QWave.h"
struct Eigstruct
{
        double lamda;
        VectorXd state;
};



bool comp(const Eigstruct& a, const Eigstruct& b);



QWave::QWave(const Sub& sys, const SingleSub& m, const SingleSub& n, const Sub& env):
_Wave(16),
_Dim(16)
{
        for(auto its=sys.SysEye().PDim().begin(); its!=sys.SysEye().PDim().end();++its)
        {
                for(auto itm=m.System().PDim().begin(); itm!=m.System().PDim().end();++itm)
                {
                        for(auto itn=n.System().PDim().begin(); itn!=n.System().PDim().end();++itn)
                        {
                                for(auto ite=env.SysEye().PDim().begin(); ite!=env.SysEye().PDim().end();++ite)
                                {
                                        int ord(order(its->first, itm->first, itn->first, ite->first));
                                        _Dim.at(ord).SDim=its->second;
                                        _Dim.at(ord).MDim=itm->second;
                                        _Dim.at(ord).NDim=itn->second;
                                        _Dim.at(ord).EDim=ite->second;
                                        MatrixXd tempM(its->second, ite->second);
                                        vector<vector<MatrixXd>> temp2;
                                        for(int mm=0; mm<itm->second; ++mm)
                                        {
                                                vector<MatrixXd> temp1;
                                                for(int nn=0; nn<itn->second; ++nn)
                                                {
                                                        temp1.push_back(tempM);
                                                }
                                                temp2.push_back(temp1);
                                        }

                                        _Wave.at(ord)=temp2;
                                }

                        }
                }
        }
}


//====================Trunc===========================================
void QWave::SMENTruncL(OP& A, const int& D, const Parity& pari)const
{
        
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        vector<map<string, int>> tempdim(4);//to store the dimension for sys(0)... Env(3) of different parity.
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                for(int k=0; i<4; ++k)
                {
                        auto it=tempdim.at(k).find(temps.at(k));
                        if(it==tempdim.at(k).end())
                        {
                                int tempd;
                                switch (k)
                                {
                                        case 0:
                                        {
                                                tempd=_Dim.at(i).SDim;
                                                break;
                                        }
                                        case 1:
                                        {
                                                tempd=_Dim.at(i).MDim;
                                                break;
                                        }
                                        case 2:
                                        {
                                                tempd=_Dim.at(i).NDim;
                                                break;
                                        }
                                        case 3:
                                        {
                                                tempd=_Dim.at(i).EDim;
                                                break;
                                        }
                                }
                                
                                tempdim.at(k).insert(pair<string, int>(temps.at(k), tempd));
                        }
                }
        }
}








//============================Reshape===================================
/*void QWave::SMEN2Wave(const MatrixXd& A)
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
*/

//===========================================================================
void QWave::SysOPWave(const OP& Sys, const Parity& pari)
{
        auto it=Sys.PRL().begin();
        if(it->first!=it->second)
        {
                cerr<<"The SysOPWave can't be used for the Sysop have changed the parity!!"<<endl;
                exit(true);
        }
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<_Dim.at(i).MDim; ++im)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        //cout<<temps.size()<<endl;        
                                        _Wave.at(i).at(im).at(in)=Sys.PMat().at(temps.at(0))*_Wave.at(i).at(im).at(in);
                                }
                        }
                }
        }
}

void QWave::SysOPWave(const OP& Sys, const QWave& wave, const Parity& pari)
{

        auto it=Sys.PRL().begin();
        if(it->first!=it->second)
        {
                cerr<<"The SysOPWave can't be used for the Sysop have changed the parity!!"<<endl;
                exit(true);
        }
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<_Dim.at(i).MDim; ++im)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        //cout<<temps.size()<<endl;        
                                        _Wave.at(i).at(im).at(in)=Sys.PMat().at(temps.at(0))*wave._Wave.at(i).at(im).at(in);
                                }
                        }
                }
        }
}
void QWave::EnvOPWave(const OP& Env, const QWave& wave, const Parity& pari)
{
        auto it=Env.PRL().begin();
        if(it->first!=it->second)
        {
                cerr<<"The EnvOPWave can't be used for the Sysop have changed the parity!!"<<endl;
                exit(true);
        }
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<_Dim.at(i).MDim; ++im)
                        {
                                for(int in=0; in<_Dim.at(i).MDim; ++in)
                                {
                                        _Wave.at(i).at(im).at(in)+=wave.Wave().at(i).at(im).at(in)*Env.PMat().at(temps.at(3)).transpose();
                                }
                        }
                }
        }
}

void QWave::MOPWave(const SOP& M, const QWave& wave, const Parity& pari)
{
        auto it=M.PRL().begin();
        if(it->first!=it->second)
        {
                cerr<<"The MOPWave can't be used for the Sysop have changed the parity!!"<<endl;
                exit(true);
        }
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int in=0; in<wave._Dim.at(i).NDim; ++in)
                        {
                                for (int k=0; k<M.PMat().at(temps.at(1)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(M.PMat().at(temps.at(1)),k); it; ++it)
                                {
                
                                        _Wave.at(i).at(it.row()).at(in)+=it.value()*wave.Wave().at(i).at(it.col()).at(in);
                                }
                
                
                        }

                }
        }
}

void QWave::NOPWave(const SOP& N, const QWave& wave, const Parity& pari)
{       
        auto it=N.PRL().begin();
        if(it->first!=it->second)
        {
                cerr<<"The NOPWave can't be used for the Sysop have changed the parity!!"<<endl;
                exit(true);
        }
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<wave._Dim.at(i).MDim; ++im)
                        {
                                for (int k=0; k<N.PMat().at(temps.at(2)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(N.PMat().at(temps.at(2)),k); it; ++it)
                                {
                
                                        _Wave.at(i).at(im).at(it.row())+=it.value()*wave.Wave().at(i).at(im).at(it.col());
                                }
                
                
                        }

                }
        }
}

void QWave::SysMOPWave(const OP& Sys, const SOP& M, const QWave& wave, const double& g, const Parity& pari)
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int in=0; in<wave._Dim.at(i).NDim; ++in)
                        {
                                for (int k=0; k<M.PMat().at(temps.at(1)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(M.PMat().at(temps.at(1)),k); it; ++it)
                                {
                                        int j(order(Sys.PRL().at(temps.at(0)),M.PRL().at(temps.at(1)) , temps.at(2), temps.at(3)));
                                        _Wave.at(j).at(it.row()).at(in)+=g*it.value()*Sys.PMat().at(temps.at(0))*wave.Wave().at(i).at(it.col()).at(in);       
        
                                }
                        }
                }
        }
}

void QWave::EnvMOPWave(const OP& Env, const SOP& M, const QWave& wave, const double& g, const Parity& pari)
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int in=0; in<wave._Dim.at(i).NDim; ++in)
                        {
                                for (int k=0; k<M.PMat().at(temps.at(1)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(M.PMat().at(temps.at(1)),k); it; ++it)
                                {
                                        int j(order(temps.at(0),M.PRL().at(temps.at(1)) , temps.at(2), Env.PRL().at(temps.at(3))));
                                        _Wave.at(j).at(it.row()).at(in)+=g*it.value()*wave.Wave().at(i).at(it.col()).at(in)*Env.PMat().at(temps.at(3)).transpose();       
        
                                }
                        }
                }
        }
}



void QWave::EnvNOPWave(const OP& Env, const SOP& N, const QWave& wave, const double& g, const Parity& pari)
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<wave._Dim.at(i).MDim; ++im)
                        {
                                for (int k=0; k<N.PMat().at(temps.at(2)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(N.PMat().at(temps.at(2)),k); it; ++it)
                                {
                                        int j(order(temps.at(0),temps.at(1) , N.PRL().at(temps.at(2)), Env.PRL().at(temps.at(3))));
                                        _Wave.at(j).at(im).at(it.row())+=g*it.value()*wave.Wave().at(i).at(im).at(it.col())*Env.PMat().at(temps.at(3)).transpose();       
        
                                }
                        }
                }
        }
}



void QWave::SysNOPWave(const OP& Sys, const SOP& N, const QWave& wave, const double& g, const Parity& pari)
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int im=0; im<wave._Dim.at(i).MDim; ++im)
                        {
                                for (int k=0; k<N.PMat().at(temps.at(2)).outerSize(); ++k)
                                for (SparseMatrix<double>::InnerIterator it(N.PMat().at(temps.at(2)),k); it; ++it)
                                {
                                        int j(order(Sys.PRL().at(temps.at(0)),temps.at(1) , N.PRL().at(temps.at(2)), temps.at(3)));
                                        _Wave.at(j).at(im).at(it.row())+=g*it.value()*Sys.PMat().at(temps.at(0))*wave.Wave().at(i).at(im).at(it.col());       
        
                                }
                        }
                }
        }
}


void QWave::add(const QWave& Wave, const Parity& pari)
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {


                        bool ts(_Dim.at(i).SDim!=Wave._Dim.at(i).SDim);
                        bool te(_Dim.at(i).EDim!=Wave._Dim.at(i).EDim);
                        bool tm(_Dim.at(i).MDim!=Wave._Dim.at(i).MDim);
                        bool tn(_Dim.at(i).NDim!=Wave._Dim.at(i).NDim);

                        if(ts|te|tm|tn)
                        {
                                cerr<<"The wave add error!!"<<endl;
                                exit(true);
                        }


                        for(int im=0; im<_Dim.at(i).MDim; ++im)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        _Wave.at(i).at(im).at(in)+=Wave._Wave.at(i).at(im).at(in);
                                }
                        }

                }

        }
}

void QWave::Hamiltanian(const Sub& Sys, const SingleSub& M, const SingleSub& N, const Sub& Env, const QWave& wave, const Parameter& para, const Parity& pari)
{
        SysOPWave(Sys.System(), pari);
        MOPWave(M.System(), wave, pari);
        NOPWave(N.System(), wave, pari);
        EnvOPWave(Env.System(), wave, pari);

        SysMOPWave(Sys.SysA(), M.SysAdag(), wave, para.gr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysA(), wave, para.gr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysAdag(), wave, para.gcr(), pari);
        SysMOPWave(Sys.SysA(), M.SysA(), wave, para.gcr(), pari);


        EnvMOPWave(Env.SysA1(), M.SysAdag(), wave, para.gr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysA(), wave, para.gr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysAdag(), wave, para.gcr(), pari);
        EnvMOPWave(Env.SysA1(), M.SysA(), wave, para.gcr(), pari);


        EnvNOPWave(Env.SysA(), N.SysAdag(), wave, para.gr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysA(), wave, para.gr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysAdag(), wave, para.gcr(), pari);
        EnvNOPWave(Env.SysA(), N.SysA(), wave, para.gcr(), pari);

        SysNOPWave(Sys.SysA1(), N.SysAdag(), wave, para.gr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysA(), wave, para.gr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysAdag(), wave, para.gcr(), pari);
        SysNOPWave(Sys.SysA1(), N.SysA(), wave, para.gcr(), pari);


}

QWave::QWave(const Sub& Sys, const SingleSub& M, const SingleSub& N, const Sub& Env, const QWave& wave, const Parameter& para, const Parity& pari):
QWave(Sys, M, N, Env)
{
        SysOPWave(Sys.System(), wave, pari);
        MOPWave(M.System(), wave, pari);
        NOPWave(N.System(), wave, pari);
        EnvOPWave(Env.System(), wave, pari);

        SysMOPWave(Sys.SysA(), M.SysAdag(), wave, para.gr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysA(), wave, para.gr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysAdag(), wave, para.gcr(), pari);
        SysMOPWave(Sys.SysA(), M.SysA(), wave, para.gcr(), pari);


        EnvMOPWave(Env.SysA1(), M.SysAdag(), wave, para.gr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysA(), wave, para.gr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysAdag(), wave, para.gcr(), pari);
        EnvMOPWave(Env.SysA1(), M.SysA(), wave, para.gcr(), pari);


        EnvNOPWave(Env.SysA(), N.SysAdag(), wave, para.gr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysA(), wave, para.gr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysAdag(), wave, para.gcr(), pari);
        EnvNOPWave(Env.SysA(), N.SysA(), wave, para.gcr(), pari);

        SysNOPWave(Sys.SysA1(), N.SysAdag(), wave, para.gr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysA(), wave, para.gr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysAdag(), wave, para.gcr(), pari);
        SysNOPWave(Sys.SysA1(), N.SysA(), wave, para.gcr(), pari);


}
//===========================================================================

void QWave::Wave2f(vector<double>& f, const Parity& pari)const
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                        {
                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                {
                                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                                        {
                                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                                {
                                                        f.push_back((_Wave.at(i).at(im).at(in))(is, ie));
                                                }
                                        }
                                }
                        }
                }        
                
        }

}


void QWave::f2Wave(const vector<double>& f, const Parity& pari)
{
        int fi(0);//for(int i=0; i<f.size(); ++i)cout<<f.at(i)<<endl;


        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)/2);
                if(nnpari==npari)
                {
                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                        {
                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                {
                                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                                        {
                                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                                {
                                                        _Wave.at(i).at(im).at(in)(is, ie)=f.at(fi++);
                                                }
                                        }
                                }
                        }
                }        
                
        }

}

/*
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
}*/
