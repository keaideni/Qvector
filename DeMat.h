/*************************************************************************
	> File Name: DeMat.h
	> Author: 
	> Mail: 
	> Created Time: 2018 Мог 21, Пү 10:24:03
 ************************************************************************/

#ifndef _DEMAT_H
#define _DEMAT_H
#include<Eigen/Dense>
//using namespace std;
using namespace Eigen;
class DeMat: public MatrixXd
{
public:
        DeMat(const int& i, const int& j):MatrixXd(i, j){};

        ~DeMat(){};


        double& insert(const int& i, const int& j)
        {
                return(*this)(i, j);
        }
};
/*ostream& operator<<(ostream& os, const DeMat& a)
{
        operator<<(os, cast<MatrixXd>(a));
        return os;
}*/
#endif
