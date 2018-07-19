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


        //DMRG haha(para, Positive);

        /*ofstream outfile("./result/ResultP");
        outfile.precision(20);
        outfile<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        <<" ,Jcr= "<<para.Jcr()<<" ,Energy= "<<haha.FEnergy()<<endl;
        outfile.close();*/

        //DMRG hehe(para, Negative);

        //outfile.open("./result/ResultN");
        //outfile.precision(20);
        //outfile<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        //<<" ,Jcr= "<<para.Jcr()<<" ,Energy= "<<hehe.FEnergy()<<endl;
        //outfile.close();
        /*<<" ,ExcitedEnergy= "<<haha.Excited()<<" ,Entropy= "
        <<haha.Entropy()<<endl<<" ,AParticleNo= "<<ParticleNo(para)<<" ,SigmaParticleNo= "
	<<SigmaParticleNo(para)<<" ,<A>= "<<OrderParameter(para)
	<<" ,SecondCorrelation= "<<secondcorrelation(para)
	<<" ,Parity= "<<Parity(para)<<endl;*/
        ofstream outfile1("./result/ParticleP");
	std::cout<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        <<" ,Jcr= "<<para.Jcr()<<" ,AParticleNo= "<<ParticleNo(para, outfile1)<<endl;
        outfile1.close();
        outfile1.open("./result/SigmaParticleP");
	std::cout<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        <<" ,Jcr= "<<para.Jcr()<<" ,AParticleNo= "<<SigmaParticleNo(para, outfile1)<<endl;;
        outfile1.close();
        
        outfile1.open("./result/CorrelationP");
        Correlation(para, outfile1);
        outfile1.close();

        outfile1.open("./result/SigmaCorrelationP");
        SigmaCorrelation(para, outfile1);
        outfile1.close();

	/*<<" ,SigmaParticleNo= "<<SigmaParticleNo(para)<<" ,<A>= "<<OrderParameter(para)
	<<" ,SecondCorrelation= "<<secondcorrelation(para)
	<<" ,Parity= "<<Parity(para)<<endl;*/

        /*Correlation(para);
	SigmaCorrelation(para);*/
}
