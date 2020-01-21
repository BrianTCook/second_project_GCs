#  There are a number of code options that need to be switched on/off
#  on compile-time.
#
#  From the following list, please select the desired options:
OPT1     =  -DWARN        # to enable or disable warnings
OPT2     =  -DDIM6        # Does calculations in 6 dimensional space
#OPT2    =  -DDIM3         # Does calculations in 3 dimensional space
#OPT2    =  -DDIMO         # specify dimension in allvars.h (#define ND 4) 
#OPT3    =  -DMEDIAN       # reduces memory, if SPH not needed
#OPT4    =  -DPERIODIC     # Does calculations in 6 dimensional space


OPTIMIZE =    -O3 -Wall 
# Extra for intel machines
#OPTIMIZE =    -O3 -Wall -ffast-math -msse2
# Extra for MAC OS X
#OPTIMIZE =    -O3 -Wall -fast  

OPTIONS =  $(OPTIMIZE) $(OPT1) $(OPT2) $(OPT3) $(OPT4) $(OPT5) $(OPT6)\
           $(OPT7) $(OPT8) $(OPT9)  $(OPT10)

CXX       = g++  #   
CXXFLAGS =   $(OPTIONS)
LIBS   =  -lm 

OBJS   = main.o  begrun.o   init.o   io.o allvars.o functions.o  \
	allocate.o  read_ic.o  nr.o    tree.o tree_search.o  ngb_search.o density_nd.o

INCL   = allvars.h proto.h  nr.h functions.h ngb_search.h Makefile 


EXEC   = Enbid

$(EXEC): $(OBJS) 
	$(CXX) $(OBJS) $(LIBS) -o  $(EXEC)  
	mv $(EXEC) ..

$(OBJS): $(INCL) 

.PHONY : clean
clean:
	rm -f $(OBJS) 

