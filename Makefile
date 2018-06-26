## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=g++ -std=c++11 -g
#########


## Flags to give the compiler for "release mode"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -I/home/xuejian/WORK/Lib/spectra/include/ -I/home/xuejian/WORK/Lib/eigen-git-mirror/





obj=main.o OP.o SOP.o Sub.o QWave.o #Super.o DMRG.o
main:$(obj)
	$(CCCOM) -o main $(obj)  $(LIBSPECTRA)
main.o:main.cpp  test.h SingleSub.h #DMRG.h Calcu.h
	$(CCCOM) -c main.cpp -O2 $(LIBSPECTRA)
OP.o:OP.cpp OP.h
	$(CCCOM) -c OP.cpp -O2 $(LIBSPECTRA)
SOP.o:SOP.cpp SOP.h
	$(CCCOM) -c SOP.cpp -O2 $(LIBSPECTRA)
Sub.o:Sub.cpp Sub.h Parameter.h
	$(CCCOM) -c Sub.cpp -O2 $(LIBSPECTRA)
QWave.o:QWave.cpp QWave.h
	$(CCCOM) -c QWave.cpp -O2 $(LIBSPECTRA)
#Super.o:Super.cpp Super.h
#	$(CCCOM) -c Super.cpp -O2 $(LIBSPECTRA)
#DMRG.o:DMRG.cpp DMRG.h SuperEnergy.h Trunc.h
#	$(CCCOM) -c DMRG.cpp -O2 $(LIBSPECTRA)
.PHONY:clean
clean:
	rm -f main $(obj)















