/*************************************************************************
    > File Name: SOP.h
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年06月05日 星期二 09时58分06秒
 ************************************************************************/
#ifndef SOP_H
#define SOP_H
#include "Parameter.h"

#include <Eigen/Sparse>
#include <cmath>

#include <unordered_map>

using namespace std;
using namespace Eigen;
enum SOPType{
	Creation, Annihilation, Iden, SigmaZ, SigmaP, SigmaM, SigmaI
};
typedef SparseMatrix<double> SpMat;
 
class SOP
{
private:
	unordered_map<string, SpMat> _PMat;
	unordered_map<string, string> _PRL;
	unordered_map<string, int> _PDim;
	void MatrixKron(SpMat& ab, const SpMat& a, const SpMat& b)
        {
        ab.setZero();
        ab.resize(a.rows()*b.rows(), a.cols()*b.cols());

        for (int k=0; k<a.outerSize(); ++k)
                for(SpMat::InnerIterator it(a,k); it; ++it)
                {
                        for(int l=0; l<b.outerSize(); ++l)
                        {
                                for(SpMat::InnerIterator itt(b,l); itt; ++itt)
                                {
                                        ab.insert(it.row()*b.rows()+itt.row(), 
                                                it.col()*b.cols()+itt.col())
                                        =it.value()*itt.value();
                                }
                        }
                };
        }

        void Block(SpMat& a, const SpMat& b, const int& startr, const int& startc)
        {
                for(int k=0; k<b.outerSize(); ++k)
                {
                        for(SpMat::InnerIterator it(b,k); it; ++it)
                        {
                                a.inset(startr+it.row(), startc+it.col());
                        }
                }
        }

	const bool Matexist(const SOP& a, const string& b);//whther b is exist or not.
public:
	const unordered_map<string, SpMat>& PMat()const{
	return _PMat;
	};

	const unordered_map<string, string>& PRL()const{
	return _PRL;
	};
	const unordered_map<string, int>& PDim()const{
	return _PDim;
	};

	SOP(const Parameter& para, const SOPType& type);
	SOP(const SOP& a, const SOP& b)
	{
		Kron(a, b);
	};

	SOP(const SOP& a):
	_PMat(a._PMat),
	_PRL(a._PRL),
	_PDim(a._PDim)
	{

	};

	SOP(){

	};
	const SOP& Kron(const SOP& a, const SOP& b);

	const SOP& add(const SOP& a);

	const SOP operator+(const SOP& a)const;
	const SOP& operator=(const SOP& a);
	const SOP& time(const SOP& a, const SOP& b);
	const SOP operator*(const SOP& a)const;
	const SOP& time(const double& d);
	const SOP operator*(const double& d)const;

	

	void show()const
	{
		if(_PRL.at("positive")=="positive")
		{
			cout<<"This is a parity no change SOPerator"<<endl;
		}else
		{
			
			cout<<"This is a parity changed SOPerator"<<endl;
		}

		cout<<"===========The dimension for each parity=========="<<endl;
		cout<<"The positive => "<<_PDim.at("positive")<<endl;
		cout<<"THe negative => "<<_PDim.at("negative")<<endl;
		cout<<"===========The matrix for different parity========="<<endl;
		auto it=_PMat.find("positive");
		if(it!=_PMat.end())cout<<"The positive => "<<endl<<_PMat.at("positive").rows()<<"X"<<_PMat.at("positive").cols()<<endl;
		auto itt=_PMat.find("negative");
		if(itt!=_PMat.end())
		cout<<"The negative => "<<endl<<_PMat.at("negative").rows()<<"X"<<_PMat.at("negative").cols()<<endl;


		
	}

	

};

const SOP operator*(const double& d, const SOP& a);
#endif
