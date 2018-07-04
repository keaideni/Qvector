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
#include<Eigen/Eigenvalues>
using namespace std;
using namespace Eigen;
enum Parity
{
        Positive, Negative
};

 
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
	const OP& Kron(const OP& a, const OP& b);

	const OP& add(const OP& a);

	const OP operator+(const OP& a)const;
	const OP& operator=(const OP& a);
	const OP& time(const OP& a, const OP& b);
	const OP operator*(const OP& a)const;
	const OP& time(const double& d);
	const OP operator*(const double& d)const;
        //To calculate an OP time wave; left |xy> to |z>.
        const OP& LWavetime(const OP& a, const OP& wave);
        //To calculate a wave*TruncU time OP; right |xy> to |z>.
        const OP& RWavetime(const OP& wave, const OP& a);
        //To calculate wave*TruncU.transpose. right |z> to |xy>.
        const OP& RWavetime2(const OP& wave, const OP& a);
        //To calculate TruncU*wave. left |z> to |xy>.
        const OP& LWavetime2(const OP& wave, const OP& a);

	
	void save(ofstream& outfile)const;
	void read(ifstream& infile);
        
        void DenTruncL(const OP& A, const int& D, double& err);//OP must be a wave OP in QWave.h.
        void DenTruncR(const OP& A, const int& D, double& err);
        void TruncU(const OP& A);
        void TruncSave(const int& orbital)const;//To save the truncation operator;
        void TruncRead(const int& orbital);


	void show()const
	{
		if(_PRL.at("positive")=="positive")
		{
			cout<<"This is a parity no change operator"<<endl;
		}else
		{
			
			cout<<"This is a parity changed operator"<<endl;
		}

		//cout<<"===========The dimension for each parity=========="<<endl;
		//cout<<"The positive => "<<_PDim.at("positive")<<endl;
		//cout<<"THe negative => "<<_PDim.at("negative")<<endl;
		cout<<"===========The matrix for different parity========="<<endl;
		auto it=_PMat.find("positive");
		if(it!=_PMat.end())cout<<"The positive => "<<endl<<_PMat.at("positive").rows()<<"X"<<_PMat.at("positive").cols()<<endl;
		auto itt=_PMat.find("negative");
		if(itt!=_PMat.end())
		cout<<"The negative => "<<endl<<_PMat.at("negative").rows()<<"X"<<_PMat.at("negative").cols()<<endl;


		
	}

        friend class QWave;	

};

const OP operator*(const double& d, const OP& a);
#endif
