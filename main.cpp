#include "DMRG.h"
//#include "test.h"
#include "Calcu.h"


int Sub::nmax;

int main(void)
{
        Parameter para;
        Sub::nmax=para.nmax();
        DMRG haha(para);
        ofstream outfile("./result/Result");
        outfile.precision(20);
        outfile<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        <<" ,Jcr= "<<para.Jcr()<<" ,Energy= "<<haha.FEnergy()<<" ,ExcitedEnergy= "<<haha.Excited()<<" ,Entropy= "
        <<haha.Entropy()<<endl<<" ,AParticleNo= "<<ParticleNo(para)<<" ,SigmaParticleNo= "
	<<SigmaParticleNo(para)<<" ,<A>= "<<OrderParameter(para)
	<<" ,SecondCorrelation= "<<secondcorrelation(para)
	<<" ,Parity= "<<Parity(para)<<endl;
	/*std::cout<<"gr= "<<para.gr()<<" ,gcr= "<<para.gcr()<<" ,Jr= "<<para.Jr()
        <<" ,Jcr= "<<para.Jcr()<<" ,AParticleNo= "<<ParticleNo(para)
	<<" ,SigmaParticleNo= "<<SigmaParticleNo(para)<<" ,<A>= "<<OrderParameter(para)
	<<" ,SecondCorrelation= "<<secondcorrelation(para)
	<<" ,Parity= "<<Parity(para)<<endl;*/

        Correlation(para);
	SigmaCorrelation(para);
}
