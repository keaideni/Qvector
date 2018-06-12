/*************************************************************************
    > File Name: OP.cpp
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年06月05日 星期二 11时05分26秒
 ************************************************************************/

#include "OP.h"

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
			for(int i=0; i<_PDim.at("positive"); ++i)
			{
				if(i<_PDim.at("negative"))
				tempp(i, i+1)=sqrt(2*(i+1));
			}
			for(int i=0; i<_PDim.at("negative"); ++i)
			{
				if(i+1<_PDim.at("positive"))
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
		case SingmaZ:
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
		case SingmaP:
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

		case SingmaM:
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

		case SingmaI:
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
		outfile<<it->first<<"\t"<<endl;
	}

	outfile.precision(20);
    outfile<<_PMat.size()<<endl;
	for(auto it=_PMat.begin(); it!=_PMat.end();++it)
	{
		outfile<<it->first<<endl
			<<it->second<<endl;
	}
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
        MatrixXd temp(MatrixXd::Zero(_PDim.at(Rname), _PDim.at(_PRL.at(Rname))));

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
