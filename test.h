#ifndef TEST_H
#define TEST_H

#include "OP.h"

void test();
void test()
{
	Parameter para;
	OP haha(para, Annihilation);
	//haha.show();
	OP hehe(para, Creation);
	//hehe.show();

	//hehe.add(haha);
	//hehe.show();
	OP heha(hehe*haha*12);
        	

	heha.show();

        ofstream outfile("heha");
        heha.save(outfile);

        ifstream infile("heha");
        OP heheha;
        heheha.read(infile);
        heheha.show();
}






#endif // TEST_H
