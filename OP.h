/*************************************************************************
    > File Name: OP.h
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年06月05日 星期二 09时58分06秒
 ************************************************************************/
#ifndef OP_H
#define OP_H
#include "Parameter.h"

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>

#include<sstream>
#include<string>
#include <unordered_map>

using namespace std;
using namespace Eigen;
enum OpType{
	Creation, Annihilation, Iden, SingmaZ, SingmaP, SingmaM, SingmaI
};

enum 
class OP
{
private:
	unordered_map<string, MatrixXd> _PMat;
	unordered_map<string, string> _PRL;
	unordered_map<string, int> _PDim;
	void MatrixKron(MatrixXd& ab,const MatrixXd& a, const MatrixXd& b)
	{
		ab=MatrixXd::Zero(a.rows()*b.rows(),a.cols()*b.cols());

		int sizer(b.rows()),sizec(b.cols());
	        for(int i=0; i<a.rows(); ++i)
		{
			for(int j=0; j<a.cols(); ++j)
			{
				int startr(i*b.rows()),startc(j*b.cols());

	                        ab.block(startr,startc,sizer, sizec)=a(i,j)*b;
		        }
		}
	}
	
	const OP& Kron(const OP& a, const OP& b);
	const bool Matexist(const OP& a, const string& b);//whther b is exist or not.
public:
	const unordered_map<string, MatrixXd>& PMat()const{
	return _PMat;
	};

	const unordered_map<string, string>& PRL()const{
	return _PRL;
	};
	const unordered_map<string, int>& PDim()const{
	return _PDim;
	};

	OP(const Parameter& para, const OpType& type);
	OP(const OP& a, const OP& b)
	{
		Kron(a, b);
	};

	OP(const OP& a):
	_PMat(a._PMat),
	_PRL(a._PRL),
	_PDim(a._PDim)
	{

	};

	OP(){

	};

	const OP& add(const OP& a);

	const OP operator+(const OP& a)const;
	const OP& operator=(const OP& a);
	const OP& time(const OP& a, const OP& b);






	void show()const
	{
		if(_PRL.at("positive")=="positive")
		{
			cout<<"This is a parity no change operator"<<endl;
		}else
		{
			
			cout<<"This is a parity changed operator"<<endl;
		}

		cout<<"===========The dimension for each parity=========="<<endl;
		cout<<"The positive => "<<_PDim.at("positive")<<endl;
		cout<<"THe negative => "<<_PDim.at("negative")<<endl;
		cout<<"===========The matrix for different parity========="<<endl;
		auto it=_PMat.find("positive");
		if(it!=_PMat.end())cout<<"The positive => "<<endl<<_PMat.at("positive")<<endl;
		auto itt=_PMat.find("negative");
		if(itt!=_PMat.end())
		cout<<"The negative => "<<endl<<_PMat.at("negative")<<endl;


		
	}

	

};

#endif
