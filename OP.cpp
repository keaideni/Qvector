/*************************************************************************
    > File Name: OP.cpp
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年06月05日 星期二 11时05分26秒
 ************************************************************************/

#include "OP.h"
string itos(const int& i)
{
        stringstream s;
        s<<i;
        return s.str();
}


OP::OP(const Parameter& para, const OpType& type)
{
	switch (type)
	{
		case Creation:
		{
				
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			MatrixXd tempp(MatrixXd::Zero(_PDim.at("negative"), _PDim.at("positive")));
			MatrixXd tempn(MatrixXd::Zero(_PDim.at("positive"), _PDim.at("negative")));
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				if(i<_PDim.at("negative"))
				tempp(i, i)=sqrt(2*i+1);
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				if(i+1<_PDim.at("positive"))
				tempn(i+1, i)=sqrt(2*(i+1));
			}

			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
			break;
		}
		case Annihilation:
		{

			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			MatrixXd tempp(MatrixXd::Zero(_PDim.at("negative"), _PDim.at("positive")));
			MatrixXd tempn(MatrixXd::Zero(_PDim.at("positive"), _PDim.at("negative")));
			for(int i=0; i<_PDim.at("positive")-1; ++i)
			{
				if(i<_PDim.at("negative"))
				tempp(i, i+1)=sqrt(2*(i+1));
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				if(i<_PDim.at("positive"))
				tempn(i, i)=sqrt(2*i+1);
			}

			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
			break;
		}
		case Iden:
		{

			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			MatrixXd tempp(MatrixXd::Zero(_PDim.at("positive"), _PDim.at("positive")));
			MatrixXd tempn(MatrixXd::Zero(_PDim.at("negative"), _PDim.at("negative")));
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				//if(i<_PDim.at("negative"))
				tempp(i, i)=1;
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				//if(i+1<_PDim.at("positive"))
				tempn(i, i)=1;
			}

			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
			break;
		}
		case SigmaZ:
		{
			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			MatrixXd tempp(MatrixXd::Identity(1,1));
			MatrixXd tempn(MatrixXd::Identity(1,1));
			tempp*=-1;
			
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
			break;



		}
		case SigmaP:
		{
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			MatrixXd tempp(MatrixXd::Identity(1,1));
			//MatrixXd tempn(MatrixXd::Identity(1,1));
			//tempp*=-1;
			
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			//_PMat.insert(pair<string, MatrixXd>("negative", tempn));
			break;


		}

		case SigmaM:
		{
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			MatrixXd tempp(MatrixXd::Identity(1,1));
			//MatrixXd tempn(MatrixXd::Identity(1,1));
			//tempp*=-1;
			
			//_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempp));
			break;


		}

		case SigmaI:
		{
			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			MatrixXd tempp(MatrixXd::Identity(1,1));
			MatrixXd tempn(MatrixXd::Identity(1,1));
			//tempp*=-1;
			
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));

			break;

		}
	}
}



const OP& OP::Kron(const OP& a, const OP& b)
{
        _PRL.clear();
        _PDim.clear();
        _PMat.clear();
	if(a._PRL.at("positive")==b._PRL.at("positive"))
	{
		_PRL.insert(pair<string, string>("positive", "positive"));
		_PRL.insert(pair<string, string>("negative", "negative"));

	}else
	{

		_PRL.insert(pair<string, string>("positive", "negative"));
		_PRL.insert(pair<string, string>("negative", "positive"));
	}

	int pp(a._PDim.at("positive")*b._PDim.at("positive")),
	    pn(a._PDim.at("positive")*b._PDim.at("negative")),
	    np(a._PDim.at("negative")*b._PDim.at("positive")),
	    nn(a._PDim.at("negative")*b._PDim.at("negative"));

	
	
	_PDim.insert(pair<string, int>("positive", pp+nn));
	_PDim.insert(pair<string, int>("negative", pn+np));
	

	MatrixXd tempp(MatrixXd::Zero(_PDim.at(_PRL.at("positive")), _PDim.at("positive")));
	MatrixXd tempn(MatrixXd::Zero(_PDim.at(_PRL.at("negative")), _PDim.at("negative")));
	bool tell1, tell2;
	if(_PRL.at("positive")=="positive")
	{
		if(a._PRL.at("positive")=="positive")
		{
			MatrixXd temp;
			{
				auto ita=a._PMat.find("positive");
				auto itb=b._PMat.find("positive");
				tell1=ita!=(a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				tempp.block(0,0, pp, pp)=temp;}
			}

			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				tempp.block(pp, pp, nn, nn)=temp;}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				tempn.block(0,0, pn, pn)=temp;}
			}
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				tempn.block(pn, pn, np, np)=temp;}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
		}else
		{
			MatrixXd temp;
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				tempp.block(0,pp, pp, nn)=temp;}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				tempp.block(pp, 0, nn, pp)=temp;}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));

			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				tempn.block(0,pn, pn, np)=temp;}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2)
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				tempn.block(pn, 0, np, pn)=temp;
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));

		}
	}else
	{
		if(a._PRL.at("positive")=="positive")
		{
			
			MatrixXd temp;
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				tempp.block(0,0, pn, pp)=temp;}
			}
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				tempp.block(pn, pp, np, nn)=temp;}
			}
			
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				tempn.block(0,0, pp, pn)=temp;}
			}
			{
				
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				tempn.block(pp, pn, nn, np)=temp;}
			}

			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));
		}else
		{
			MatrixXd temp;
			{	
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				tempp.block(0,pp, pn, nn)=temp;}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				tempp.block(pn, 0, np, pp)=temp;}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("positive", tempp));


			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				tempn.block(0,pn, pp, np)=temp;}
			}

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				tempn.block(pp, 0, nn, pn)=temp;}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, MatrixXd>("negative", tempn));

		}

	}
}



const OP& OP::add(const OP& a)
{
	if(_PRL.at("positive")!=a._PRL.at("positive"))
	{
		cerr<<"Error: The PRL of the two operator conflict!"<<endl;
		exit(true);
	}
	if(_PDim.at("positive")!=a._PDim.at("positive")|_PDim.at("negative")!=a._PDim.at("negative"))
	{
		cerr<<"Error: The dimension for each parity conflict!"<<endl;
		exit(true);
	}
	
	if((! Matexist(*this, "positive"))&(Matexist(a, "positive")))
	{
		_PMat.at("positive")=a._PMat.at("positive");
	}else if(Matexist(*this, "positive")&Matexist(a, "positive"))
	{
		_PMat.at("positive")+=a._PMat.at("positive");

	}

	if((! Matexist(*this, "negative"))&(Matexist(a, "negative")))
	{
		_PMat.at("negative")=a._PMat.at("negative");
	}else if(Matexist(*this, "negative")&Matexist(a, "negative"))
	{
		_PMat.at("negative")+=a._PMat.at("negative");

	}
	return *this;
	
}

const double OP::trace()const
{
        for(auto it=_PRL.begin(); it!=_PRL.end(); ++it)
        {
                if(it->first!=it->second)
                {
                        return 0;
                }
        }
        double temp(0);
        for(auto it=_PMat.begin(); it!=_PMat.end(); ++it)
        {
                temp+=it->second.trace();
        }

        return temp;
}


const OP OP::adjoint()const
{
        OP temp;

        temp._PDim=(*this)._PDim;

        for(auto it=_PRL.begin(); it!=_PRL.end(); ++it)
        {
                temp._PRL.insert(pair<string, string>(it->first, it->second));

                auto itt=_PMat.find(it->first);
                if(itt!=_PMat.end())
                        temp._PMat.insert(pair<string, MatrixXd>(it->second, itt->second));
        }

        return temp;
}

const OP OP::operator+(const OP& a)const
{
	OP aa(*this);

	return aa.add(a);

	//aa.show();

	//return aa;

}




const bool OP::Matexist(const OP& a, const string& b)
{
	auto it=a._PMat.find(b);
	return it!=a._PMat.end();
}



const OP& OP::operator=(const OP& a)
{
	_PDim=a._PDim;
	_PRL=a._PRL;
	_PMat=a._PMat;

	return *this;
}



const OP& OP::time(const OP& a, const OP& b)
{
        _PDim.clear();
        _PRL.clear();
        _PMat.clear();
	_PDim=(a._PDim);
	
	if(_PDim.at("positive")!=a._PDim.at("positive")|_PDim.at("negative")!=a._PDim.at("negative"))
	{
		cerr<<"Error: The dimension for each parity conflict!"<<endl;
		exit(true);
	}

	for(auto itb=b._PMat.begin(); itb!=b._PMat.end(); ++itb)
	{
		auto ita=a._PMat.find(b._PRL.at(itb->first));
		if(ita!=a._PMat.end())
			_PMat.insert(pair<string, MatrixXd>(itb->first, ita->second*itb->second));
	}
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
		_PRL.insert(pair<string, string>(itb->first, a._PRL.at(itb->second)));
	}
	return *this;
}


const OP& OP::LWavetime2(const OP& a, const OP& b)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
	        for(auto ita=a._PRL.begin(); ita!=a._PRL.end(); ++ita)
                {
		        if(itb->second==ita->first)
                        {
			        _PMat.insert(pair<string, MatrixXd>(itb->first, a._PMat.at(ita->first)*b._PMat.at(itb->first)));

		                _PRL.insert(pair<string, string>(itb->first, ita->second));

                        }
                }
	}
	return *this;
}



const OP& OP::LWavetime(const OP& a, const OP& b)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
	        for(auto ita=a._PRL.begin(); ita!=a._PRL.end(); ++ita)
                {
		        if(itb->second==ita->second)
                        {
                                _PMat.insert(pair<string, MatrixXd>(itb->first, a._PMat.at(ita->first).transpose()*b._PMat.at(itb->first)));

		                _PRL.insert(pair<string, string>(itb->first, ita->first));
                        }	
                }
	}
	return *this;
}

const OP& OP::RWavetime(const OP& a, const OP& b)
{
        _PDim.clear();
        _PRL.clear();
        _PMat.clear();
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
		auto ita=a._PMat.find(b._PRL.at(itb->first));
		if(ita!=a._PMat.end())
                {
        	        _PMat.insert(pair<string, MatrixXd>(itb->first, a._PMat.at(ita->first)*b._PMat.at(itb->first)));

		        _PRL.insert(pair<string, string>(itb->first, a._PRL.at(ita->first)));                
                }
	
	}
	return *this;
}


const OP& OP::RWavetime2(const OP& a, const OP& b)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
	        for(auto ita=a._PRL.begin(); ita!=a._PRL.end(); ++ita)
                {
		        if(itb->first==ita->first)
                        {
                	        _PMat.insert(pair<string, MatrixXd>(itb->second,  a._PMat.at(ita->first)*b._PMat.at(itb->first).transpose()));

		                _PRL.insert(pair<string, string>(itb->second, ita->second));                
                        }
		
                }
	}
	return *this;
}


const double OP::AverageL(const OP& wave, const OP& a)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
	for(auto itwave=wave._PRL.begin(); itwave!=wave._PRL.end(); ++itwave)
	{
	        for(auto ita=a._PRL.begin(); ita!=a._PRL.end(); ++ita)
                {
		        if(itwave->second==ita->first)
                        {
                                for(auto itwavepri=wave._PRL.begin(); itwavepri!=wave._PRL.end();++itwavepri)
                                {
                                        if(itwavepri->second==ita->second)
                                        {
                	                        _PMat.insert(pair<string, MatrixXd>(itwave->first, wave._PMat.at(itwavepri->first).adjoint()*a._PMat.at(ita->first)*wave._PMat.at(itwave->first)));
		                                _PRL.insert(pair<string, string>(itwave->first, itwavepri->first));                
                                        }
                                }

                        }
		
                }
	}
	return trace();
}

const double OP::AverageR(const OP& wave, const OP& a)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
	for(auto itwave=wave._PRL.begin(); itwave!=wave._PRL.end(); ++itwave)
	{
	        for(auto ita=a._PRL.begin(); ita!=a._PRL.end(); ++ita)
                {
		        if(itwave->first==ita->first)
                        {
                	        _PMat.insert(pair<string, MatrixXd>(itwave->second, wave._PMat.at(ita->second)*a._PMat.at(ita->first)*wave._PMat.at(itwave->first).adjoint()));

		                _PRL.insert(pair<string, string>(itwave->second, wave._PRL.at(ita->second)));                
                        }
		
                }
	}
	return trace();
}




const double OP::Average(const OP& wave, const OP& OS, const OP& O)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
        OP waveadjoint(wave.adjoint());
        OP OE(O.adjoint());

	for(auto itwave=wave._PRL.begin(); itwave!=wave._PRL.end(); ++itwave)
	{

                _PMat.insert(pair<string, MatrixXd>(itwave->first, OE._PMat.at(waveadjoint.PRL().at(OS.PRL().at(itwave->second)))*waveadjoint._PMat.at(OS.PRL().at(itwave->second)).adjoint()*OS._PMat.at(itwave->second)*wave._PMat.at(itwave->first)));

		_PRL.insert(pair<string, string>(itwave->first, OE.PRL().at(waveadjoint.PRL().at(OS.PRL().at(itwave->second)))));
			
	}
	return trace();
}




const OP OP::operator*(const OP& a)const
{
	OP b;
	b.time(*this, a);
	return b;
}
	






const OP& OP::time(const double& d)
{
	for(auto it=_PMat.begin(); it!=_PMat.end();++it)
		it->second*=d;

	return *this;
}
	


const OP OP::operator*(const double& d)const
{
	OP temp(*this);
	return temp.time(d);
}



const OP operator*(const double& d, const OP& a)
{
	return a*d;
}




void OP::save(ofstream& outfile)const
{
        outfile<<_PRL.size()<<endl;
	for(auto it=_PRL.begin(); it!=_PRL.end(); ++it)
	{
		outfile<<it->first<<"\t"<<it->second<<endl;
	}
        outfile<<_PDim.size()<<endl;
	for(auto it=_PDim.begin(); it!=_PDim.end(); ++it)
	{
		outfile<<it->first<<"\t"<<it->second<<endl;
	}

	outfile.precision(20);
         outfile<<_PMat.size()<<endl;
	for(auto it=_PMat.begin(); it!=_PMat.end();++it)
	{
		outfile<<it->first<<endl
			<<it->second<<endl;
	}
}

void OP::TruncSave(const int& orbital)const
{
        string filename("./trunc/");
        filename+=itos(orbital);
        ofstream outfile(filename);
        outfile<<_PRL.size()<<endl;
	for(auto it=_PRL.begin(); it!=_PRL.end(); ++it)
	{
		outfile<<it->first<<"\t"<<it->second<<endl;
	}
        outfile<<_PDim.size()<<endl;
	for(auto it=_PDim.begin(); it!=_PDim.end(); ++it)
	{
		outfile<<it->first<<"\t"<<it->second<<endl;
	}

	outfile.precision(20);
         outfile<<_PMat.size()<<endl;
	for(auto it=_PMat.begin(); it!=_PMat.end();++it)
	{

		outfile<<it->first<<endl
                        <<it->second.rows()<<endl
                        <<it->second.cols()<<endl
			<<it->second<<endl;
	

        }
        outfile.close();
}



void OP::read(ifstream& infile)
{
	_PRL.clear();
	_PDim.clear();
	_PMat.clear();


	string Rname, Lname;
	int Dim, tempsize;
        infile>>tempsize;
        for(int i=0; i<tempsize; ++i)
        {
       
                infile>>Rname>>Lname;
	        _PRL.insert(pair<string, string>(Rname, Lname));
        }

         infile>>tempsize;
         for(int i=0; i<tempsize; ++i)
        {
                infile>>Rname>>Dim;
	        _PDim.insert(pair<string, int>(Rname, Dim));
        }

        infile>>tempsize;
        for(int i=0; i<tempsize; ++i)
        {
                infile>>Rname;
                MatrixXd temp(MatrixXd::Zero( _PDim.at(_PRL.at(Rname)), _PDim.at(Rname)));

                for(int j=0; j<temp.rows(); ++j)
                {
                 for(int k=0; k<temp.cols(); ++k)
                {
                        infile>>temp(j, k);
                }
                }
	        _PMat.insert(pair<string, MatrixXd>(Rname, temp));
        }
}

void OP::TruncRead(const int& orbital)
{
        _PDim.clear();
        _PRL.clear();
        _PMat.clear();
        string filename("./trunc/");
        filename+=itos(orbital);

        ifstream infile(filename);
        if(!infile)
        {
                cerr<<"the file "<<filename<<" doesn't exist!!!!!"<<endl;
                exit(true);
        }
        _PRL.clear();
	_PDim.clear();
	_PMat.clear();


	string Rname, Lname;
	int Dim, tempsize;
        infile>>tempsize;
        for(int i=0; i<tempsize; ++i)
        {
       
                infile>>Rname>>Lname;
	        _PRL.insert(pair<string, string>(Rname, Lname));
        }

         infile>>tempsize;
         for(int i=0; i<tempsize; ++i)
        {
                infile>>Rname>>Dim;
	        _PDim.insert(pair<string, int>(Rname, Dim));
        }

        infile>>tempsize;
        for(int i=0; i<tempsize; ++i)
        {
                infile>>Rname;
                int temprow, tempcol;
                infile>>temprow>>tempcol;
                MatrixXd temp(MatrixXd::Zero(temprow, tempcol));
                
                for(int j=0; j<temp.rows(); ++j)
                {
                 for(int k=0; k<temp.cols(); ++k)
                {
                        infile>>temp(j, k);
                }
                }
	        _PMat.insert(pair<string, MatrixXd>(Rname, temp));
        }

        infile.close();

}

struct Eigstruct
{
        double lambda;
        string pari; 
        VectorXd state;
};


bool comp(const Eigstruct& a, const Eigstruct& b);
bool comp(const Eigstruct& a, const Eigstruct& b)
{
        return a.lambda>b.lambda;
}

void OP::DenTruncL(const OP& A, const int& D, double& err)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
        vector<Eigstruct> Denmat;
        for(auto it=A._PMat.begin(); it!=A._PMat.end(); ++it)
        {
                _PRL.insert(pair<string, string>(A._PRL.at(it->first), A._PRL.at(it->first)));
                MatrixXd temp(it->second*it->second.transpose());
                SelfAdjointEigenSolver<MatrixXd> es(temp);
                for(int i=0; i<es.eigenvalues().size(); ++i)
                {
                        Eigstruct tempe={es.eigenvalues()(i), A._PRL.at(it->first), es.eigenvectors().col(i)};
                        Denmat.push_back(tempe);
                }
        }
        sort(Denmat.begin(), Denmat.end(), comp);

        int min(Denmat.size()<D?Denmat.size():D);

        for(int i=0; i<min; ++i)
        {
                auto it=_PDim.find(Denmat.at(i).pari);
                if(it!=_PDim.end())
                {
                        it->second+=1;
                }else
                {
                        _PDim.insert(pair<string, int>(Denmat.at(i).pari, 1));
                }
        }

        for(auto it=_PDim.begin(); it!=_PDim.end(); ++it)
        {
                int matR(it->second);
                int matL;
                for(int i=0; i<min; ++i)
                {
                        if(Denmat.at(i).pari==it->first)
                        {
                                matL=Denmat.at(i).state.size();
                                break;
                        }
                }
                MatrixXd temp(matL, matR);
                int order(0);
                for(int i=0; i<min; ++i)
                {
                        if(Denmat.at(i).pari==it->first)
                        {
                                temp.col(order++)=Denmat.at(i).state;
                        }
                }
                _PMat.insert(pair<string, MatrixXd>(it->first, temp));
        }
        double sum(0);
        double sum1(0);
        for(int i=0; i<Denmat.size(); ++i)
        {
                sum+=Denmat.at(i).lambda;
                if(i<min)sum1+=Denmat.at(i).lambda;
        }
        err=sum-sum1;
        
}
void OP::DenTruncR(const OP& A, const int& D, double& err)
{
        _PRL.clear();
        _PMat.clear();
        _PDim.clear();
        vector<Eigstruct> Denmat;
        for(auto it=A._PMat.begin(); it!=A._PMat.end(); ++it)
        {
                _PRL.insert(pair<string, string>(A._PRL.at(it->first), A._PRL.at(it->first)));
                MatrixXd temp(it->second.transpose()*it->second);
                SelfAdjointEigenSolver<MatrixXd> es(temp);
                for(int i=0; i<es.eigenvalues().size(); ++i)
                {
                        Eigstruct tempe={es.eigenvalues()(i), it->first, es.eigenvectors().col(i)};
                        Denmat.push_back(tempe);
                }
        }
        sort(Denmat.begin(), Denmat.end(), comp);

        int min(Denmat.size()<D?Denmat.size():D);

        for(int i=0; i<min; ++i)
        {
                auto it=_PDim.find(Denmat.at(i).pari);
                if(it!=_PDim.end())
                {
                        it->second+=1;
                }else
                {
                        _PDim.insert(pair<string, int>(Denmat.at(i).pari, 1));
                }
        }

        for(auto it=_PDim.begin(); it!=_PDim.end(); ++it)
        {
                int matR(it->second);
                int matL;
                for(int i=0; i<min; ++i)
                {
                        if(Denmat.at(i).pari==it->first)
                        {
                                matL=Denmat.at(i).state.size();
                                break;
                        }
                }
                MatrixXd temp(matL, matR);
                int order(0);
                for(int i=0; i<min; ++i)
                {
                        if(Denmat.at(i).pari==it->first)
                        {
                                temp.col(order++)=Denmat.at(i).state;
                        }
                }
                _PMat.insert(pair<string, MatrixXd>(it->first, temp));
        }
        double sum(0);
        double sum1(0);
        for(int i=0; i<Denmat.size(); ++i)
        {
                sum+=Denmat.at(i).lambda;
                if(i<min)sum+=Denmat.at(i).lambda;
        }
        err=sum-sum1;
        
        
}


void OP::TruncU(const OP& A)
{
        _PDim=A._PDim;
        for(auto it=_PRL.begin(); it!=_PRL.end(); ++it)
        {
                _PMat.at(it->first)=A._PMat.at(it->second).transpose()*_PMat.at(it->first)*A._PMat.at(it->first);
        }
}
