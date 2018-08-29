## Set which compiler to use by defining CCCOM:
##GNU GCC compiler
#CCCOM=g++ -m64 -std=c++11 
##Clang compiler (good to use on Mac OS)
#CCCOM=clang++ -std=c++1y
##Intel C++ compiler (good to use with Intel MKL if available)
CCCOM=g++ -std=c++11 -g
#########
TGT = main
SRCS = OP.cpp SOP.cpp Sub.cpp QWave.cpp SuperEnergy.cpp DMRG.cpp main.cpp
CUR_DIR = $(shell pwd)
## Flags tcpp give the ccppmpiler fcppr "release mcppde"



#LIBFLAGS = -larmadillo
LIBSPECTRA = -O2 -I/home/keaideni/WORK/Lib/spectra/include/ -I/home/keaideni/WORK/Lib/eigen-git-mirror/



$(TGT) : $(SRCS:.cpp=.o)
	$(CCCOM) $^ -o $@ $(LIBSPECTRA) 
%.o:%.cpp
	$(CCCOM) -c $<   $(LIBSPECTRA) 

%.d:%.cpp
	@$(CCCOM) -MM $< > $@ 

sinclude $(SRCS:.cpp=.d)


clean:
	rm -f $(TGT) $(SRCS:.cpp=.o) $(SRCS:.cpp=.d)

.PHONY:clean














