#ifndef TEST_H
#define TEST_H

#include "SuperEnergy.h"

void test(Parameter& para);
void test(Parameter& para)
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
        //cout<<para.Energy<<endl;
        //OP a(para, SigmaI);
        //OP b(para, Annihilation);
        //OP c(b, a);
        //c.show();

        //SingleSub ha(para);
        
        Sub hehe(para, 1);//hehe.Show();
        Parity pari(Positive);
        Sub haha(para, hehe, hehe, 2);
        haha.Show();
        //Super sup(para, hehe, haha, haha, hehe, pari);
        //SuperEnergy supp(para, sup, pari);
        //OP U;
        //supp.wave.Wave2SMEN(U, pari);
        //OP UU;
        //double err;
        //UU.DenTruncL(U, para.D(), err);

        //Sub a(para, hehe, hehe, 2);
        //a.Trunc(UU);

        //a.Show();
}






#endif // TEST_H
