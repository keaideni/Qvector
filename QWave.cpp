//!!!!!!!!the three functions must be coherent!!!!!!!!
//wave2f, f2wave, f2wave!!!!!!!!!



#include "QWave.h"






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


//====================Transform===========================================
void QWave::Wave2NSME(OP& A, const Parity& pari)const
{
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(1), temps.at(3)));
                        string PL(QAdd(temps.at(0), temps.at(2)));
                        auto it=A._PRL.find(PR);
                        if(it==A._PRL.end())
                               A._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=A._PMat.find(PR);
                        if(itt==A._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0))*dimoes.at(2).at(temps.at(2)));
                                dimL+=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(2).at(OppS(temps.at(2)));
                                int dimR(dimoes.at(3).at(temps.at(3))*dimoes.at(1).at(temps.at(1)));
                                dimR+=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(1).at(OppS(temps.at(1)));
                                MatrixXd tempmat(MatrixXd::Zero(dimL, dimR));
                                A._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int startL, startR;
                                                        if(temps.at(1)=="positive")
                                                                startR=0;
                                                        else startR=dimoes.at(1).at(OppS(temps.at(1)))*dimoes.at(3).at((OppS(temps.at(3))));
                                                        if(temps.at(2)=="positive")
                                                                startL=0;
                                                        else startL=dimoes.at(2).at(OppS(temps.at(2)))*dimoes.at(0).at(OppS(temps.at(0)));
                                                        //cout<<startL<<endl;
                                                        //cout<<startR<<endl;
                                                        //cout<<startL+in*_Dim.at(i).SDim+is<<endl<< startR+im*_Dim.at(i).EDim+ie<<endl;
                                                        A._PMat.at(PR)(startL+in*_Dim.at(i).SDim+is, startR+im*_Dim.at(i).EDim+ie)=_Wave.at(i).at(im).at(in)(is, ie);
                                                }
                                        }
                                }
                        }
                }
        }

}



void QWave::Wave2SMEN(OP& A, const Parity& pari)const
{
        
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(2), temps.at(3)));
                        string PL(QAdd(temps.at(0), temps.at(1)));
                        auto it=A._PRL.find(PR);
                        if(it==A._PRL.end())
                               A._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=A._PMat.find(PR);
                        if(itt==A._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0))*dimoes.at(1).at(temps.at(1)));
                                dimL+=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                int dimR(dimoes.at(3).at(temps.at(3))*dimoes.at(2).at(temps.at(2)));
                                dimR+=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                MatrixXd tempmat(dimL, dimR);
                                A._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int startL, startR;
                                                        if(temps.at(3)=="positive")
                                                                startR=0;
                                                        else startR=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                                        if(temps.at(0)=="positive")
                                                                startL=0;
                                                        else startL=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                                        //cout<<startL+is*_Dim.at(i).MDim+im<<endl<< startR+ie*_Dim.at(i).NDim+in<<endl;
                                                        A._PMat.at(PR)(startL+is*_Dim.at(i).MDim+im, startR+ie*_Dim.at(i).NDim+in)=_Wave.at(i).at(im).at(in)(is, ie);
                                                }
                                        }
                                }
                        }
                }
        } 
}

void QWave::SMEN2Wave(const OP& A, const Parity& pari)
{
        
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(2), temps.at(3)));
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int startL, startR;
                                                        if(temps.at(3)=="positive")
                                                                startR=0;
                                                        else startR=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                                        if(temps.at(0)=="positive")
                                                                startL=0;
                                                        else startL=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                                        //cout<<startL+is*_Dim.at(i).MDim+im<<endl<< startR+ie*_Dim.at(i).NDim+in<<endl;
                                                        _Wave.at(i).at(im).at(in)(is, ie)=A._PMat.at(PR)(startL+is*_Dim.at(i).MDim+im, startR+ie*_Dim.at(i).NDim+in);
                                                }
                                        }
                                }
                        }
                }
        } 
}


void QWave::NSME2Wave(const OP& A, const Parity& pari)
{
       
        QWave test(*this);
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        //double sum(0);
        //double sum1(0);
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(1), temps.at(3)));
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int startL, startR;
                                                        if(temps.at(1)=="positive")
                                                                startR=0;
                                                        else startR=dimoes.at(1).at(OppS(temps.at(1)))*dimoes.at(3).at((OppS(temps.at(3))));
                                                        if(temps.at(2)=="positive")
                                                                startL=0;
                                                        else startL=dimoes.at(2).at(OppS(temps.at(2)))*dimoes.at(0).at(OppS(temps.at(0)));
                                                        //cout<<startL+in*_Dim.at(i).SDim+is<<endl<< startR+im*_Dim.at(i).EDim+ie<<endl;
                                                        //sum+=abs(_Wave.at(i).at(im).at(in)(is, ie)-A._PMat.at(PR)(startL+in*_Dim.at(i).SDim+is, startR+im*_Dim.at(i).EDim+ie));
                                                        //cout<<sum<<endl;
                                                        //OP test(A);
                                                        _Wave.at(i).at(im).at(in)(is,ie)=A.PMat().at(PR)(startL+in*_Dim.at(i).SDim+is, startR+im*_Dim.at(i).EDim+ie);
                                                        //sum1+=abs(test._Wave.at(i).at(im).at(in)(is,ie)-_Wave.at(i).at(im).at(in)(is,ie));
                                                      //cout<<"haha"<<sum1<<endl;
                                                }
                                        }
                                }
                        }
                }
        } 
}

const string QWave::OppS(const string& s)const
{
        if(s=="positive")
        {
                return "negative";
        }else
        {
                return "positive";
        }
}

const string QWave::QAdd(const string& a, const string& b)const
{
                        double tellp1(2*(s2n(a)-1./2)), tellp2(2*(s2n(b)-1./2));
                        int tellP((int)(tellp1*tellp2));
                        if(tellP==1)
                        {
                                return "positive";
                        }else
                        {
                                return "negative";
                        }
}
const vector<unordered_map<string, int>> QWave::DimOES()const
{
        vector<unordered_map<string, int>> tempdim(4);//to store the dimension for sys(0)... Env(3) of different parity.
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                for(int k=0; k<4; ++k)
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
        return tempdim;
}

//========================To reshape the waveop by the QWave========================
//=======================This part works for DMRG::Initialize=======================


void QWave::d1g1(OP& tar, const OP& sta, const Parity& pari)const
{
        tar._PMat.clear();
        tar._PDim.clear();
        tar._PRL.clear();

        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(temps.at(3));
                        string PL(QAdd(temps.at(0), temps.at(2)));
                        string staPR(QAdd(temps.at(3), temps.at(2)));
                        auto it=tar._PRL.find(PR);
                        if(it==tar._PRL.end())
                               tar._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=tar._PMat.find(PR);
                        if(itt==tar._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0))*dimoes.at(2).at(temps.at(2)));
                                dimL+=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(2).at(OppS(temps.at(2)));
                                int dimR(dimoes.at(3).at(temps.at(3)));
                                MatrixXd tempmat(dimL, dimR);
                                tar._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                                //cout<<PR<<"=>"<<dimL<<"x"<<dimR<<endl;
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                //for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                //{
                                                        int tarstartL, stastartR;
                                                        if(temps.at(3)=="positive")
                                                                stastartR=0;
                                                        else stastartR=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                                        if(temps.at(2)=="positive")
                                                                tarstartL=0;
                                                        else tarstartL=dimoes.at(2).at(OppS(temps.at(2)))*dimoes.at(0).at(OppS(temps.at(0)));
                                                        //cout<<tarstartL+in*_Dim.at(i).SDim+is<<"\t"<< ie<<"\t"<<is<<"\t"<<stastartR+ie*_Dim.at(i).NDim+in<<endl;
                                                        tar._PMat.at(PR)(tarstartL+in*_Dim.at(i).SDim+is,ie)=sta._PMat.at(staPR)(is, stastartR+ie*_Dim.at(i).NDim+in);
                                                //}
                                        }
                                }
                        }
                }
        } 
}


void QWave::dmg1(OP& tar, const OP& sta, const Parity& pari)const
{
        tar._PMat.clear();
        tar._PDim.clear();
        tar._PRL.clear();

        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(3), temps.at(1)));
                        string PL(temps.at(0));
                        string staPR(temps.at(3));
                        auto it=tar._PRL.find(PR);
                        if(it==tar._PRL.end())
                               tar._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=tar._PMat.find(PR);
                        if(itt==tar._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0)));
                                int dimR(dimoes.at(3).at(temps.at(3))*dimoes.at(1).at(temps.at(1)));
                                dimR+=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(1).at(OppS(temps.at(1)));
                                MatrixXd tempmat(dimL, dimR);
                                tar._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                //for(int in=0; in<_Dim.at(i).NDim; ++in)
                                //{
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int tarstartR, stastartL;
                                                        if(temps.at(1)=="positive")
                                                                tarstartR=0;
                                                        else tarstartR=dimoes.at(1).at(OppS(temps.at(1)))*dimoes.at(3).at(OppS(temps.at(3)));
                                                        if(temps.at(0)=="positive")
                                                                stastartL=0;
                                                        else stastartL=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                                        //cout<<startL+is*_Dim.at(i).MDim+im<<endl<< startR+ie*_Dim.at(i).NDim+in<<endl;
                                                        tar._PMat.at(PR)(is,tarstartR+im*_Dim.at(i).EDim+ie)=sta._PMat.at(staPR)(stastartL+is*_Dim.at(i).MDim+im, ie);
                                                }
                                        }
                                //}
                        }
                }
        } 
}




void QWave::d1gm(OP& tar, const OP& sta, const Parity& pari)const
{
        tar._PMat.clear();
        tar._PDim.clear();
        tar._PRL.clear();

        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                //cout<<i<<endl;
                vector<string> temps(n2s(i));
                int nnpari((i/8+(i%8)/4+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(temps.at(3));
                        string PL(QAdd(temps.at(0), temps.at(1)));
                        string staPR(QAdd(temps.at(1), temps.at(3)));
                        auto it=tar._PRL.find(PR);
                        if(it==tar._PRL.end())
                               tar._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=tar._PMat.find(PR);
                        if(itt==tar._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0))*dimoes.at(1).at(temps.at(1)));
                                dimL+=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                int dimR(dimoes.at(3).at(temps.at(3)));
                                MatrixXd tempmat(dimL, dimR);
                                tar._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                                //cout<<PR<<"=>"<<dimL<<"x"<<dimR<<endl;
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                //for(int in=0; in<_Dim.at(i).NDim; ++in)
                                //{
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                {
                                                        int tarstartL, stastartR;
                                                        if(temps.at(1)=="positive")
                                                                stastartR=0;
                                                        else stastartR=dimoes.at(1).at(OppS(temps.at(1)))*dimoes.at(3).at(OppS(temps.at(3)));
                                                        if(temps.at(0)=="positive")
                                                                tarstartL=0;
                                                        else tarstartL=dimoes.at(0).at(OppS(temps.at(0)))*dimoes.at(1).at(OppS(temps.at(1)));
                                                        //cout<<startL+is*_Dim.at(i).MDim+im<<endl<< startR+ie*_Dim.at(i).NDim+in<<endl;
                                                        //cout<<tarstartL+is*_Dim.at(i).MDim+im<<"\t"<<ie<<"\t"<<is<<"\t"<<im*_Dim.at(i).EDim+ie<<endl;
                                                        tar._PMat.at(PR)(tarstartL+is*_Dim.at(i).MDim+im,ie)=sta._PMat.at(staPR)(is, stastartR+im*_Dim.at(i).EDim+ie);
                                                }
                                        }
                                //}
                        }
                }
        } 
}


void QWave::dmgm(OP& tar, const OP& sta, const Parity& pari)const
{
        tar._PMat.clear();
        tar._PDim.clear();
        tar._PRL.clear();

        int npari;
        if(pari==Positive)npari=0;
        else npari=1;
        vector<unordered_map<string, int>> dimoes(DimOES());
        for(int i=0; i<16; ++i)
        {
                vector<string> temps(n2s(i));
                int nnpari((i/8+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        string PR(QAdd(temps.at(3), temps.at(2)));
                        string PL(temps.at(0));
                        string staPR(temps.at(3));
                        auto it=tar._PRL.find(PR);
                        if(it==tar._PRL.end())
                               tar._PRL.insert(pair<string, string>(PR, PL));

                        auto itt=tar._PMat.find(PR);
                        if(itt==tar._PMat.end())
                        {
                                int dimL(dimoes.at(0).at(temps.at(0)));
                                int dimR(dimoes.at(3).at(temps.at(3))*dimoes.at(2).at(temps.at(2)));
                                dimR+=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                MatrixXd tempmat(dimL, dimR);
                                tar._PMat.insert(pair<string, MatrixXd>(PR, tempmat));
                        }
                        for(int ie=0; ie<_Dim.at(i).EDim; ++ie)
                        {
                                for(int in=0; in<_Dim.at(i).NDim; ++in)
                                {
                                        for(int is=0; is<_Dim.at(i).SDim; ++is)
                                        {
                                                //for(int im=0; im<_Dim.at(i).MDim; ++im)
                                                //{
                                                        int tarstartR, stastartL;
                                                        if(temps.at(3)=="positive")
                                                                tarstartR=0;
                                                        else tarstartR=dimoes.at(3).at(OppS(temps.at(3)))*dimoes.at(2).at(OppS(temps.at(2)));
                                                        if(temps.at(2)=="positive")
                                                                stastartL=0;
                                                        else stastartL=dimoes.at(2).at(OppS(temps.at(2)))*dimoes.at(0).at(OppS(temps.at(0)));
                                                        //cout<<startL+is*_Dim.at(i).MDim+im<<endl<< startR+ie*_Dim.at(i).NDim+in<<endl;
                                                        tar._PMat.at(PR)(is,tarstartR+ie*_Dim.at(i).NDim+in)=sta._PMat.at(staPR)(stastartL+in*_Dim.at(i).SDim+is, ie);
                                                //}
                                        }
                                }
                        }
                }
        } 
}
//==================================================================================

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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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

        SysMOPWave(Sys.SysA(), M.SysAdag(), wave, -0.5*para.Jr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysA(), wave, -0.5*para.Jr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysAdag(), wave, -0.5*para.Jcr(), pari);
        SysMOPWave(Sys.SysA(), M.SysA(), wave, -0.5*para.Jcr(), pari);


        EnvMOPWave(Env.SysA1(), M.SysAdag(), wave, -0.5*para.Jr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysA(), wave, -0.5*para.Jr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysAdag(), wave, -0.5*para.Jcr(), pari);
        EnvMOPWave(Env.SysA1(), M.SysA(), wave, -0.5*para.Jcr(), pari);


        EnvNOPWave(Env.SysA(), N.SysAdag(), wave, -0.5*para.Jr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysA(), wave, -0.5*para.Jr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysAdag(), wave, -0.5*para.Jcr(), pari);
        EnvNOPWave(Env.SysA(), N.SysA(), wave, -0.5*para.Jcr(), pari);

        SysNOPWave(Sys.SysA1(), N.SysAdag(), wave, -0.5*para.Jr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysA(), wave, -0.5*para.Jr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysAdag(), wave, -0.5*para.Jcr(), pari);
        SysNOPWave(Sys.SysA1(), N.SysA(), wave, -0.5*para.Jcr(), pari);


}

QWave::QWave(const Sub& Sys, const SingleSub& M, const SingleSub& N, const Sub& Env, const QWave& wave, const Parameter& para, const Parity& pari):
QWave(Sys, M, N, Env)
{
        clock_t start, end;
        start=clock();tsys=tsm=tsn=tem=ten=start-start;
        SysOPWave(Sys.System(), wave, pari);
        MOPWave(M.System(), wave, pari);
        NOPWave(N.System(), wave, pari);
        EnvOPWave(Env.System(), wave, pari);
        end=clock();
        tsys+=end-start;

        start=clock();
        SysMOPWave(Sys.SysA(), M.SysAdag(), wave, -0.5*para.Jr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysA(), wave, -0.5*para.Jr(), pari);
        SysMOPWave(Sys.SysAdag(), M.SysAdag(), wave, -0.5*para.Jcr(), pari);
        SysMOPWave(Sys.SysA(), M.SysA(), wave, -0.5*para.Jcr(), pari);
        end=clock();
        tsm+=end-start;

        start=clock();
        EnvMOPWave(Env.SysA1(), M.SysAdag(), wave, -0.5*para.Jr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysA(), wave, -0.5*para.Jr(), pari);
        EnvMOPWave(Env.SysAdag1(), M.SysAdag(), wave, -0.5*para.Jcr(), pari);
        EnvMOPWave(Env.SysA1(), M.SysA(), wave, -0.5*para.Jcr(), pari);
        end=clock();
        tem+=(end-start);

        start=clock();
        EnvNOPWave(Env.SysA(), N.SysAdag(), wave, -0.5*para.Jr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysA(), wave, -0.5*para.Jr(), pari);
        EnvNOPWave(Env.SysAdag(), N.SysAdag(), wave, -0.5*para.Jcr(), pari);
        EnvNOPWave(Env.SysA(), N.SysA(), wave, -0.5*para.Jcr(), pari);
        end=clock();
        ten+=end-start;

        start=clock();
        SysNOPWave(Sys.SysA1(), N.SysAdag(), wave, -0.5*para.Jr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysA(), wave, -0.5*para.Jr(), pari);
        SysNOPWave(Sys.SysAdag1(), N.SysAdag(), wave, -0.5*para.Jcr(), pari);
        SysNOPWave(Sys.SysA1(), N.SysA(), wave, -0.5*para.Jcr(), pari);
        end=clock();
        tsn+=end-start;

}
//===========================================================================

void QWave::Wave2f(vector<double>& f, const Parity& pari)const
{
        f.clear();
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
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



int QWave::waveDim(const Parity& pari)const
{
        int re(0);
        int npari;
        if(pari==Positive)npari=0;
        else npari=1;

        for(int i=0; i<16; ++i)
        {
                int nnpari((i/8+(i%8)/4+((i%8)%4)/2+((i%8)%4)%2)%2);
                if(nnpari==npari)
                {
                        re+=_Dim.at(i).SDim*_Dim.at(i).MDim*_Dim.at(i).NDim*_Dim.at(i).EDim;
                        


                }
        }
        return re;

}
