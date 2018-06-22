#ifndef TEST_H
#define TEST_H

#include "QWave.h"

void test(const Parameter& para);
void test(const Parameter& para)
{
        /*Sub haha(para, 1);
        Sub hehe(para, haha, haha, 2);


        //hehe.Show();

        hehe.Save();

        Sub heha;

        heha.Read(2);
        //heha.Show();
        //
        MatrixXd temp;
        temp=hehe.SysA().PMat().at("negative")-heha.SysA().PMat().at("negative");

        //cout<<emp.sum()<<endl;*/
        SingleSub haha(para);
        
        Sub hehe(para, 1);

        QWave wave(hehe, haha, haha, hehe);

        vector<double> f;
        wave.Wave2f(f, Positive);

        VectorXd ff(MatrixXd::Random(f.size(), 1));

        for(int i=0; i<f.size(); ++i)
        {
                f.at(i)=ff(i);
        }

        wave.f2Wave(f, Positive);

        vector<double> fff;
        wave.Wave2f(fff, Positive);
        double sum(0);
        for(int i=0; i<fff.size(); ++i)
        {
                sum+=abs(f.at(i)-fff.at(i));
        }
        cout<<sum<<endl<<f.at(2)<<endl<<fff.at(2)<<endl;
        
        QWave wave2(wave);
        //wave.SysOPWave(hehe.System(), Positive);

        wave2.Hamiltanian(hehe, haha, haha, hehe, wave, para, Positive);
        QWave wave3(hehe, haha, haha, hehe, wave, para, Positive);
        
}






#endif // TEST_H
