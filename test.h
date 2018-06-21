#ifndef TEST_H
#define TEST_H

#include "SOP.h"

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
        SOP haha(para, Creation);
        haha.show();
}






#endif // TEST_H
