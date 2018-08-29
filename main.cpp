#include "DMRG.h"
//#include "test.h"
#include "Calcu.h"


//int Sub::nmax;

int main(void)
{
	//test()
        Parameter para;
        //Sub::nmax=para.nmax();
        //test(para);

        ofstream outfile, outfile1;
        DMRG haha(para, Positive);

        outfile.open("./result/ResultP");
        outfile.precision(20);
        outfile<<"gr=\t"<<para.gr()<<"\t gcr=\t"<<para.gcr()<<"\t Jr=\t"<<para.Jr()
        <<"\t Jcr=\t"<<para.Jcr()<<"\t Energy=\t"<<haha.FEnergy()<<"\t";


        outfile1.open("./result/ParticleP");
        outfile<<"\t AParticleNo=\t"<<ParticleNo(para, outfile1)<<"\t";

        outfile1.close();
        outfile1.open("./result/SigmaParticleP");
        outfile<<"\t SigmaParticleNo=\t"<<SigmaParticleNo(para, outfile1)<<endl;
       
        outfile1.close();
        
        outfile1.open("./result/CorrelationP");
        Correlation(para, outfile1);
        outfile1.close();

        outfile1.open("./result/SigmaCorrelationP");
        SigmaCorrelation(para, outfile1);
        outfile1.close();
        outfile.close();

        /*DMRG hahaha(para, Negative);

        outfile.open("./result/ResultN");
        outfile.precision(20);
        outfile<<"gr=\t"<<para.gr()<<"\t gcr=\t"<<para.gcr()<<"\t Jr=\t"<<para.Jr()
        <<"\t Jcr=\t"<<para.Jcr()<<"\t Energy=\t"<<hahaha.FEnergy()<<"\t";

        outfile1.open("./result/ParticleN");
        outfile<<"\t AParticleNo=\t"<<ParticleNo(para, outfile1)<<"\t";

        outfile1.close();
        outfile1.open("./result/SigmaParticleN");
        outfile<<"\t SigmaParticleNo=\t"<<SigmaParticleNo(para, outfile1)<<endl;;
       
        outfile1.close();
        
        outfile1.open("./result/CorrelationN");
        Correlation(para, outfile1);
        outfile1.close();

        outfile1.open("./result/SigmaCorrelationN");
        SigmaCorrelation(para, outfile1);
        outfile1.close();
        outfile.close();*/
}
