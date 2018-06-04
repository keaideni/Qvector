//for the M and N single site. Because the matrix is sparse, so construct a class especially.
#ifndef SINGLE_SUB_H
#define SINGLE_SUB_H
#include "Parameter.h"
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat;

class SingleSub
 {

 private:
        SpMat _System;
        SpMat _SysA;
        SpMat _SysAdag;
        void Kron(SpMat& ab, const SpMat& a, const SpMat& b);

 public:
        const SpMat& System()const{return _System;};
        const SpMat& SysA()const{return _SysA;};
        const SpMat& SysAdag()const{return _SysAdag;};
         SingleSub(){};
         ~SingleSub(){};
         SingleSub(const Parameter& para);
         
 }; 


#endif // SINGLE_SUB_H
