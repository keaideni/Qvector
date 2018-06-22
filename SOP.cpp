/*************************************************************************
    > File Name: SOP.cpp
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年06月05日 星期二 11时05分26秒
 ************************************************************************/

#include "SOP.h"

SOP::SOP(const Parameter& para, const OpType& type)
{
	switch (type)
	{
		case Creation:
		{
				
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			SpMat tempp(_PDim.at("negative"), _PDim.at("positive"));
			SpMat tempn(_PDim.at("positive"), _PDim.at("negative"));
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				if(i<_PDim.at("negative"))
				tempp.insert(i, i)=sqrt(2*i+1);
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				if(i+1<_PDim.at("positive"))
				tempn.insert(i+1, i)=sqrt(2*(i+1));
			}

			_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempn));
			break;
		}
		case Annihilation:
		{

			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			SpMat tempp(_PDim.at("negative"), _PDim.at("positive"));
			SpMat tempn(_PDim.at("positive"), _PDim.at("negative"));
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				if(i<_PDim.at("negative"))
				tempp.insert(i, i+1)=sqrt(2*(i+1));
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				if(i+1<_PDim.at("positive"))
				tempn.insert(i, i)=sqrt(2*i+1);
			}

			_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempn));
			break;
		}
		case Iden:
		{

			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", para.nmax()/2+1));
			_PDim.insert(pair<string, int>("negative", (para.nmax()+1)/2));
			
			SpMat tempp(_PDim.at("positive"), _PDim.at("positive"));
			SpMat tempn(_PDim.at("negative"), _PDim.at("negative"));
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				//if(i<_PDim.at("negative"))
				tempp.insert(i, i)=1;
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				//if(i+1<_PDim.at("positive"))
				tempn.insert(i, i)=1;
			}

			_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempn));
			break;
		}
		case SigmaZ:
		{
			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			SpMat tempp(1,1);
			SpMat tempn(1,1);
			tempp.insert(0,0)=-1;
                        tempn.insert(0,0)=1;
                        
			
			_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempn));
			break;



		}
		case SigmaP:
		{
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			SpMat tempp(1,1);
			//SpMat tempn(1,1);
			tempp.insert(0,0)=1;
			
			_PMat.insert(pair<string, SpMat>("positive", tempp));
			//_PMat.insert(pair<string, SpMat>("negative", tempn));
			break;


		}

		case SigmaM:
		{
			_PRL.insert(pair<string, string>("positive", "negative"));
			_PRL.insert(pair<string, string>("negative", "positive"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			SpMat tempp(1,1);
			//SpMat tempn(1,1);
			tempp.insert(0,0)=1;
			
			//_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempp));
			break;


		}

		case SigmaI:
		{
			_PRL.insert(pair<string, string>("positive", "positive"));
			_PRL.insert(pair<string, string>("negative", "negative"));
			_PDim.insert(pair<string, int>("positive", 1));
			_PDim.insert(pair<string, int>("negative", 1));

			SpMat tempp(1,1);
			SpMat tempn(1,1);
			tempp.insert(0,0)=1;
                        tempn.insert(0,0)=1;
			
			_PMat.insert(pair<string, SpMat>("positive", tempp));
			_PMat.insert(pair<string, SpMat>("negative", tempn));

			break;

		}
	}
}



const SOP& SOP::Kron(const SOP& a, const SOP& b)
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
	

	SpMat tempp(_PDim.at(_PRL.at("positive")), _PDim.at("positive"));
	SpMat tempn(_PDim.at(_PRL.at("negative")), _PDim.at("negative"));
	bool tell1, tell2;
	if(_PRL.at("positive")=="positive")
	{
		if(a._PRL.at("positive")=="positive")
		{
			SpMat temp;
			{
				auto ita=a._PMat.find("positive");
				auto itb=b._PMat.find("positive");
				tell1=ita!=(a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				Block(tempp, temp, 0, 0);}
			}

			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				Block(tempp, temp, pp, pp);}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("positive", tempp));

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				Block(tempn, temp, 0,0);}
			}
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				Block(tempn, temp, pn, pn);}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("negative", tempn));
		}else
		{
			SpMat temp;
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				Block(tempp, temp, 0,pp);}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				Block(tempp, temp, pp, 0);}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("positive", tempp));

			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell1=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				Block(tempn, temp, 0,pn);}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell2=(ita!=a._PMat.end())&&(itb!=b._PMat.end());
				if(tell2)
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				Block(tempn, temp, pn, 0);
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("negative", tempn));

		}
	}else
	{
		if(a._PRL.at("positive")=="positive")
		{
			
			SpMat temp;
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				Block(tempp, temp, 0,0);}
			}
			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				Block(tempp, temp, pn, pp);}
			}
			
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("positive", tempp));

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				Block(tempn, temp, 0,0);}
			}
			{
				
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				Block(tempn, temp, pp, pn);}
			}

			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("negative", tempn));
		}else
		{
			SpMat temp;
			{	
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("negative");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("negative"));
				Block(tempp, temp, 0,pp);}
			}
			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("positive");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("positive"));
				Block(tempp, temp, pn, 0);}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("positive", tempp));


			{
				auto ita=a._PMat.find("negative"), itb=b._PMat.find("positive");
				tell1=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell1){
				MatrixKron(temp, a._PMat.at("negative"), b._PMat.at("positive"));
				Block(tempn, temp, 0,pn);}
			}

			{
				auto ita=a._PMat.find("positive"), itb=b._PMat.find("negative");
				tell2=ita!=a._PMat.end()&itb!=b._PMat.end();
				if(tell2){
				MatrixKron(temp, a._PMat.at("positive"), b._PMat.at("negative"));
				Block(tempn, temp, pp, 0);}
			}
			if(tell1|tell2)
			_PMat.insert(pair<string, SpMat>("negative", tempn));

		}

	}
}



const SOP& SOP::add(const SOP& a)
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


const SOP SOP::operator+(const SOP& a)const
{
	SOP aa(*this);

	return aa.add(a);

	//aa.show();

	//return aa;

}




const bool SOP::Matexist(const SOP& a, const string& b)
{
	auto it=a._PMat.find(b);
	return it!=a._PMat.end();
}



const SOP& SOP::operator=(const SOP& a)
{
	_PDim=a._PDim;
	_PRL=a._PRL;
	_PMat=a._PMat;

	return *this;
}



const SOP& SOP::time(const SOP& a, const SOP& b)
{
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
			_PMat.insert(pair<string, SpMat>(itb->first, ita->second*itb->second));
	}
	for(auto itb=b._PRL.begin(); itb!=b._PRL.end(); ++itb)
	{
		_PRL.insert(pair<string, string>(itb->first, a._PRL.at(itb->second)));
	}
	return *this;
}






const SOP SOP::operator*(const SOP& a)const
{
	SOP b;
	b.time(*this, a);
	return b;
}
	






const SOP& SOP::time(const double& d)
{
	for(auto it=_PMat.begin(); it!=_PMat.end();++it)
		it->second*=d;

	return *this;
}
	


const SOP SOP::operator*(const double& d)const
{
	SOP temp(*this);
	return temp.time(d);
}



const SOP operator*(const double& d, const SOP& a)
{
	return a*d;
}









