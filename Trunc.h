#ifndef TRUNC_H
#define TRUNC_H
struct Eigstruct
{
        double lamda;
        VectorXd state;
};



bool comp(const Eigstruct& a, const Eigstruct& b);
bool comp(const Eigstruct& a, const Eigstruct& b)
{
        return (a.lamda > b.lamda);
}
const MatrixXd TruncL(const MatrixXd& Wave, const int& D);
const MatrixXd TruncR(const MatrixXd& Wave, const int& D);
const MatrixXd DenTruncL(const MatrixXd& Wave, const int& D);
const MatrixXd DenTruncR(const MatrixXd& Wave, const int& D);
const MatrixXd DenTruncL(const MatrixXd& Wave, const int& D, double& entropy);
const MatrixXd DenTruncR(const MatrixXd& Wave, const int& D, double& entropy);



const MatrixXd TruncL(const MatrixXd& Wave, const int& D)
{

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
        MatrixXd truncU(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncU.col(i)=denmat.at(i).state;
        }

        return truncU;

}

const MatrixXd TruncR(const MatrixXd& Wave, const int& D)
{
        vector<Eigstruct> denmat;

        if(Wave.cols()*Wave.rows()>16)
        {
                BDCSVD<MatrixXd> svd(Wave, ComputeFullV);
                for(int i=0; i<svd.singularValues().rows(); ++i)
                {
                        Eigstruct base;
                        base.lamda=svd.singularValues()(i);
                        base.state=svd.matrixV().col(i);

                        denmat.push_back(base);
                }
        }else
        {
                JacobiSVD<MatrixXd> svd(Wave, ComputeFullV);
                for(int i=0; i<svd.singularValues().rows(); ++i)
                {
                        Eigstruct base;
                        base.lamda=svd.singularValues()(i);
                        base.state=svd.matrixV().col(i);

                        denmat.push_back(base);
                }
        }

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.cols());
        int ncol(D<denmat.size()?D:denmat.size());
        MatrixXd truncV(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncV.col(i)=denmat.at(i).state;
        }

        return truncV;
}

const MatrixXd DenTruncL(const MatrixXd& Wave, const int& D)
{
        vector<Eigstruct> denmat;
        
        SelfAdjointEigenSolver<MatrixXd> es(Wave*Wave.transpose());
        for (int i = 0; i<es.eigenvalues().size(); i++)
        {
                Eigstruct temp = {es.eigenvalues()(i), es.eigenvectors().col(i) };
                denmat.push_back(temp);
        }
        
        

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.rows());
        int ncol(D<denmat.size()?D:denmat.size());
        MatrixXd truncU(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncU.col(i)=denmat.at(i).state;
        }

        return truncU;
}

const MatrixXd DenTruncL(const MatrixXd& Wave, const int& D, double& entropy)
{
        vector<Eigstruct> denmat;
        //double sum(0);
        entropy=0;
        
        SelfAdjointEigenSolver<MatrixXd> es(Wave*Wave.transpose());
        for (int i = 0; i<es.eigenvalues().size(); i++)
        {
                Eigstruct temp = {es.eigenvalues()(i), es.eigenvectors().col(i) };
                denmat.push_back(temp);
                //sum+=es.eigenvalues()(i);
        }
        
        

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.rows());
        int ncol(D<denmat.size()?D:denmat.size());
        MatrixXd truncU(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncU.col(i)=denmat.at(i).state;
                entropy-=denmat.at(i).lamda*log(denmat.at(i).lamda);
        }

        return truncU;
}


const MatrixXd DenTruncR(const MatrixXd& Wave, const int& D)
{
        vector<Eigstruct> denmat;
        
        SelfAdjointEigenSolver<MatrixXd> es(Wave.transpose()*Wave);
        for (int i = 0; i<es.eigenvalues().size(); i++)
        {
                Eigstruct temp = {es.eigenvalues()(i), es.eigenvectors().col(i) };
                denmat.push_back(temp);
        }
        
        

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.cols());
        int ncol(D<denmat.size()?D:denmat.size());
        MatrixXd truncV(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncV.col(i)=denmat.at(i).state;
        }

        return truncV;
}


const MatrixXd DenTruncR(const MatrixXd& Wave, const int& D, double& entropy)
{
        vector<Eigstruct> denmat;
        entropy=0;
        
        SelfAdjointEigenSolver<MatrixXd> es(Wave.transpose()*Wave);
        for (int i = 0; i<es.eigenvalues().size(); i++)
        {
                Eigstruct temp = {es.eigenvalues()(i), es.eigenvectors().col(i) };
                denmat.push_back(temp);
        }
        
        

        sort(denmat.begin(), denmat.end(), comp);

        int nrow(Wave.cols());
        int ncol(D<denmat.size()?D:denmat.size());
        MatrixXd truncV(MatrixXd::Zero(nrow, ncol));

        for(int i=0; i<ncol; ++i)
        {
                truncV.col(i)=denmat.at(i).state;
                entropy-=denmat.at(i).lamda*log(denmat.at(i).lamda);
        }

        return truncV;
}


#endif // TRUNC_H
