

PROJECT = int2
HOST = osx-gnu
HOST = linux-gnu

ifeq ($(HOST),linux-gnu)
  FC = gfortran
  FFLAGS = -O2 -g -w -fallow-argument-mismatch
  CC = gcc
  IDIRS = -I../include -I../src -I$(HOME)/.local/include
  CFLAGS = -std=c99 -Wall -O2 -funroll-loops -march=native
  LDFLAGS = -lopenblas -llapack -L$(HOME)/.local/lib
  LINK = gfortran -o $(PROJECT) 
endif

# ifeq ($(HOST),osx-gnu)
#   FC = gfortran-10
#   FFLAGS = -O2 -g -w -fallow-argument-mismatch
#   CC = gcc-10
#   IDIRS = -I../../utils -I../src -I$(FMM2D)/include -I/usr/local/include
#   CFLAGS = -std=c99 -Wall -O2 -funroll-loops -march=native
#   LDFLAGS = -framework Accelerate -L/usr/local/lib -lid
#   LINK = gfortran-10 -o $(PROJECT) 
# endif

# ifeq ($(HOST),osx-intel)
#   FC = ifort
#   FFLAGS = -O2 -g
#   CC = icc
#   IDIRS = -I../../utils -I../src -I$(FMM2D)/include -I/usr/local/include
#   CFLAGS = -std=c99 -Wall -O2 -funroll-loops -march=native
#   #CFLAGS = -std=c99 -O2
#   LDFLAGS = -mkl -nofor_main -L/usr/local/lib -lid
#   LINK = ifort -o $(PROJECT) 
# endif




csrcs = test_wrs2d.c ../src/cprini.c ../src/quadtree.c ../src/wrs2d.c \
  ../src/kernels.c

xsrcs = 

f90srcs = 

fsrcs = 


objs = $(csrcs:.c=.o) $(xsrcs:.cpp=.o) $(fsrcs:.f=.o) $(f90srcs:.f90=.o)
deps = $(csrcs:.c=.d) $(xsrcs:.cpp=.d)


# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)


.PHONY: all
all: $(deps) $(objs)
	rm -f $(PROJECT)
	$(LINK) $(objs) $(LDFLAGS)


%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.d: %.c
	$(CC) $(IDIRS) -MM $^ -MT $(@:.d=.o) -MF $@

%.o: %.c
	$(CC) $(CFLAGS) $(IDIRS) -c $< -o $@



.PHONY: clean
clean:
	rm -f $(objs)
	rm -f $(deps)
	rm -f $(PROJECT)
