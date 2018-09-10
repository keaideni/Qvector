#ifndef CALCU_H
#define CALCU_H
#define pi 3.1415926
#include <cmath>
#include "DMRG.h"

double ParticleNo(const Parameter& para, ofstream& outfile);
double ParticleNo(const Parameter& para, ofstream& outfile)
{
        Sub Sys, m(para,1);
        vector<OP> ParticleMat;
        ParticleMat.push_back(m.SysAdag()*m.SysA());

        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                OP tempTrunc;
                tempTrunc.TruncRead(i);
                for(int j=0; j<i-1; ++j)
                {
                        OP tempMat1(ParticleMat.at(j), m.SysEye());
                        ParticleMat.at(j)=tempMat1;
                        if(i==para.LatticeSize()/2)continue;
                        ParticleMat.at(j).TruncU(tempTrunc);

                }
                Sys.Read(i-1);

                OP tempppp(Sys.SysEye(), m.SysAdag()*m.SysA());
                ParticleMat.push_back(tempppp);
                if(i==para.LatticeSize()/2)break;
                ParticleMat.at(i-1).TruncU(tempTrunc);
        }
        ParticleMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                OP tempTrunc;
                tempTrunc.TruncRead(para.LatticeSize()-i+1);
                for(int j=0; j<i-1; ++j)
                {
                        OP tempMat2(ParticleMat.at(para.LatticeSize()/2+j), m.SysEye());
                        ParticleMat.at(para.LatticeSize()/2+j)=tempMat2;
                        if(i==para.LatticeSize()/2)continue;
                        ParticleMat.at(para.LatticeSize()/2+j).TruncU(tempTrunc);

                }
                Sys.Read(para.LatticeSize()+2-i);

                OP temppp(Sys.SysEye(), m.SysAdag()*m.SysA());
                ParticleMat.push_back(temppp);
                if(i==para.LatticeSize()/2)break;
                ParticleMat.at(para.LatticeSize()/2+i-1).TruncU(tempTrunc);
        }


        

        OP wave;
        wave.TruncRead(10000);

        vector<double> Particle;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OP tempaverage;
                Particle.push_back(tempaverage.AverageL(wave, ParticleMat.at(i)));
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                OP tempaverage;
                Particle.push_back(tempaverage.AverageR(wave, ParticleMat.at(i)));
        }
        //ofstream outfile("./result/ParticleNo");
        outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i=\t"<<(i+1)<<"\t ParticleNo=\t"<<Particle.at(i)<<endl;
                cout<<"i=\t"<<(i+1)<<"\t ParticleNo=\t"<<Particle.at(i)<<endl;
        }
        //outfile.close();

        return Particle.at(para.LatticeSize()/2);


}





double SigmaParticleNo(const Parameter& para, ofstream& outfile);
double SigmaParticleNo(const Parameter& para, ofstream& outfile)
{
			
	Sub Sys, m(para, 1);
	OP tempsigmaP(para, SigmaP), tempSigmaM(para, SigmaM);
        OP tempsigma(tempsigmaP*tempSigmaM);
        OP tempeye(para, Iden);

	OP Sigma(tempeye, tempsigma);


	vector<OP> SigmaParticleMat;
	SigmaParticleMat.push_back(Sigma);
	for(int i=2; i<=para.LatticeSize()/2; ++i)
	{
		for(int j=0; j<i-1; ++j)
		{
                        OP tempMat1(SigmaParticleMat.at(j), m.SysEye());
			SigmaParticleMat.at(j)=tempMat1;
			if(i==para.LatticeSize()/2)continue;
			OP tempTrunc;
			tempTrunc.TruncRead(i);
			SigmaParticleMat.at(j).TruncU(tempTrunc);
		}
		Sys.Read(i-1);
                OP temppp(Sys.SysEye(), Sigma);
		SigmaParticleMat.push_back(temppp);
		if(i==para.LatticeSize()/2)break;
		OP tempTrunc;
		tempTrunc.TruncRead(i);
		SigmaParticleMat.at(i-1).TruncU(tempTrunc);

	}
	SigmaParticleMat.push_back(Sigma);
	for(int i=2; i<=para.LatticeSize()/2;++i)
	{

		for(int j=0; j<i-1; ++j)
		{
                        OP tempppp(SigmaParticleMat.at(para.LatticeSize()/2+j), m.SysEye());
							
			SigmaParticleMat.at(para.LatticeSize()/2+j)=tempppp;
				
			if(i==para.LatticeSize()/2)continue;
			OP tempTrunc;
			tempTrunc.TruncRead(para.LatticeSize()-i+1);
			SigmaParticleMat.at(para.LatticeSize()/2+j).TruncU(tempTrunc);

		}
		Sys.Read(para.LatticeSize()+2-i);
                OP tempMat2(Sys.SysEye(), Sigma);
		SigmaParticleMat.push_back(tempMat2);
		if(i==para.LatticeSize()/2)break;
		OP tempTrunc;
		tempTrunc.TruncRead(para.LatticeSize()-i+1);
		SigmaParticleMat.at(para.LatticeSize()/2+i-1).TruncU(tempTrunc);
	}

        OP wave;
        wave.TruncRead(10000);

        vector<double> SigmaParticle;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OP tempaverage;
                SigmaParticle.push_back(tempaverage.AverageL(wave, SigmaParticleMat.at(i)));
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                OP tempaverage;
                SigmaParticle.push_back(tempaverage.AverageR(wave, SigmaParticleMat.at(i)));
        }
        outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i=\t"<<(i+1)<<"\tParticleNo=\t"<<SigmaParticle.at(i)<<endl;
                cout<<"i=\t"<<(i+1)<<"\t ParticleNo=\t"<<SigmaParticle.at(i)<<endl;
        }

        return SigmaParticle.at(para.LatticeSize()/2);
}









void Correlation(const Parameter& para, ofstream& outfile);
void Correlation(const Parameter& para, ofstream& outfile)
{
        Sub m(para, 0);
        vector<OP> CorrMat;
        OP Adag1(m.SysAdag());

        CorrMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OP tempmat(CorrMat.at(j), m.SysEye());
                        CorrMat.at(j)=tempmat;
                        
                        if(i==para.LatticeSize()/2)continue;
                        
                        OP tempTrunc;
                        tempTrunc.TruncRead(i);
                        CorrMat.at(j).TruncU(tempTrunc);
                        
                }
                OP tempMat1(Adag1, m.SysA());
                CorrMat.push_back(tempMat1);
                OP tempAdag1(Adag1, m.SysEye());
                Adag1=tempAdag1;
                if(i==para.LatticeSize()/2)continue;
                OP tempTrunc;
                tempTrunc.TruncRead(i);
                CorrMat.at(i-1).TruncU(tempTrunc);
                Adag1.TruncU(tempTrunc);
        }

        vector<double> Corr;
        OP wave;
        wave.TruncRead(10000);
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OP tempaverage;
                Corr.push_back(tempaverage.AverageL(wave, CorrMat.at(i)));
        }
        Sub Sys;
        Sys.Read(para.LatticeSize()/2+2);
        OP A(Sys.SysA1(), m.SysEye());

        OP averagemat;
        Corr.push_back(averagemat.Average(wave, Adag1, A));

        //ofstream outfile("./result/Correlation");
        for(int i=0; i<Corr.size(); ++i)
        {
                outfile<<"R=\t"<<i<<"\t Corr(R)=\t"<<Corr.at(i)<<endl;
                cout<<"R=\t"<<i<<"\t Corr(R)=\t"<<Corr.at(i)<<endl;
        }
        //outfile.close();
}



void SigmaCorrelation(const Parameter& para, ofstream& outfile);
void SigmaCorrelation(const Parameter& para, ofstream& outfile)
{
        Sub m(para, 0);
        vector<OP> CorrMat;
	OP tempsigmaplu(para, SigmaP ), tempsigmamin(para, SigmaM);
	
        OP tempeye(para, Iden);
        OP Adag1(tempeye, tempsigmaplu), A(tempeye, tempsigmamin);

        CorrMat.push_back(Adag1*A);
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OP tempMat(CorrMat.at(j), m.SysEye());
                        CorrMat.at(j)=tempMat;
                        
                        if(i==para.LatticeSize()/2)continue;
                        
                        OP tempTrunc;
                        tempTrunc.TruncRead(i);
                        CorrMat.at(j).TruncU(tempTrunc);
                        
                }
                OP tempMat1(Adag1, A);
                CorrMat.push_back(tempMat1);

                OP tempMat2(Adag1, m.SysEye());
                Adag1=tempMat2;
                if(i==para.LatticeSize()/2)continue;
                OP tempTrunc;
                tempTrunc.TruncRead(i);
                CorrMat.at(i-1).TruncU(tempTrunc);
                Adag1.TruncU(tempTrunc);
        }

        vector<double> Corr;
        OP wave;
        wave.TruncRead(10000);
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OP tempaverage;
                Corr.push_back(tempaverage.AverageL(wave, CorrMat.at(i)));
        }
        Sub Sys;
        Sys.Read(para.LatticeSize()/2+2);
        OP A1(Sys.SysA1(), m.SysEye());
        
        OP averagemat;
        Corr.push_back(averagemat.Average(wave, Adag1, A1));


        //Corr.push_back((wave.adjoint()*Adag1*wave*A1.transpose()).trace());

        //ofstream outfile("./result/SigmaCorrelation");
        for(int i=0; i<Corr.size(); ++i)
        {
                outfile<<"R=\t"<<i<<"\t Corr(R)=\t"<<Corr.at(i)<<endl;
                cout<<"R=\t"<<i<<"\t Corr(R)=\t"<<Corr.at(i)<<endl;
                
        }
        outfile.close();
}


/*
double secondcorrelation(const Parameter& para);
double secondcorrelation(const Parameter& para)
{
	Sub a,m(para,0);
	a.Read(para.LatticeSize()/2-1);
	MatrixXd wave;
	ReadTruncM(wave, 10000);
	MatrixXd Adag(Kron(a.SysEye(),m.SysAdag())), A(Kron(a.SysEye(),m.SysA()));
	double fenzi((wave.transpose()*Adag*Adag*A*A*wave).trace());
	double fenmu(pow((wave.transpose()*Adag*A*wave).trace(),2));

	return fenzi/fenmu;

}
*/




#endif // CALCU_H
