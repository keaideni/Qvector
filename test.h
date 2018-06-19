#ifndef TEST_H
#define TEST_H

#include "Sub.h"

void test(const Parameter& para);
void test(const Parameter& para)
{
        Sub haha(para, 1);
        Sub hehe(para, haha, haha, 2);


        //hehe.Show();

        hehe.Save();

        Sub heha;

        heha.Read(2);
        //heha.Show();
        //
        MatrixXd temp;
        temp=hehe.SysA().PMat().at("negative")-heha.SysA().PMat().at("negative");

        //cout<<emp.sum()<<endl;
}






#endif // TEST_H
