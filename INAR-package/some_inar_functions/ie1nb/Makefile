# integrate in R is linked to a source file integrate.c
# integrate.c essentially consists of the fortran code
#  in the SRC files after applying f2c
# integrate.c makes use of SEXP (S expressions) for passing a function
#   from R to C

# from dags.f, documentation for limit, lenw, iwork() and work()
#         dimensioning parameters
#            limit - integer
#                    dimensioning parameter for iwork
#                    limit determines the maximum number of subintervals
#                    in the partition of the given integration interval
#                    (a,b), limit.ge.1.
#                    if limit.lt.1, the routine will end with ier = 6.
#
#            lenw  - integer
#                    dimensioning parameter for work
#                    lenw must be at least limit*4.
#                    if lenw.lt.limit*4, the routine will end
#                    with ier = 6.


SRC=d1mach.f i1mach.f dqagi.f dqagie.f dqags.f dqagse.f dqelg.f dqk15i.f dqk21.f dqpsrt.f xerror.f xerrwv.f 

libquadpk.a: $(SRC)
	gfortran -c $(SRC)
	ar r libquadpk.a *.o
	ranlib libquadpk.a

OBJ=d1mach.o i1mach.o dqagi.o dqagie.o dqags.o dqagse.o dqelg.o dqk15i.o \
 dqk21.o dqpsrt.o xerror.o xerrwv.o

#OBJ = $(SRC:.f=.o) 

%.o: %.f
	 gfortran -c $< -o $@;

# convert .f files to .o object files
objects: $(SRC)
	gfortran -c -m64 $(SRC) 
	
# convert .f files to .o object files with fpic
objectslinux: $(SRC)
	gfortran -c -m64 -fpic $(SRC) 

# checks for inverse chf to match R code in initialchecks.r and thinop-invchf.r
# also check for I2  operator  sum_{j=0}^xprev K_j(x;alpha,gamma)
invchf: invchf.f90 $(OBJ)
	gfortran -o invchf invchf.f90 $(OBJ)

# example   invchf < invchf.in.txt > invchf.out.txt

# function ebinom to be called from fortran90 or C (see mnex1.f90, mnex2.c)
I2invchf.o: I2invchf.f90 
	gfortran -c I2invchf.f90 

FCLIB=/usr/local/gfortran/lib

# sample main program calling ebinom from C, linking f90
# fortran to C library may depend on version of Unix
ie1nb: ie1nb.c lgamma.c I2invchf.f90 $(SRC) 
	rm $(OBJ) lgamma.o I2invchf.o
	gfortran -c $(SRC)
	gfortran -c I2invchf.f90 
	gcc -c ie1nb.c lgamma.c 
	gcc -DMAIN -L$(FCLIB) -o ie1nb ie1nb.c lgamma.o I2invchf.o $(OBJ) -lgfortran -lm

Rie1nb: ie1nb.c lgamma.c I2invchf.f90 $(SRC)
	rm $(OBJ) lgamma.o I2invchf.o
	gfortran -fpic -c $(SRC) 
	gfortran -fpic -c I2invchf.f90 
	gcc -fpic -c ie1nb.c lgamma.c
	gcc -L$(FCLIB) -shared -o ie1nb.so ie1nb.o lgamma.o I2invchf.o $(OBJ) -lgfortran -lm


#======================================================================

zipfile:
	zip ie1nb.zip Makefile *.f90 *.f *.c *.r *.txt


#======================================================================
# bounded intervals, main program in f77
tstfin: tstfin.f $(OBJ) 
	gfortran -o tstfin tstfin.f $(OBJ)

# infinite intervals, main program in f77
tstinf: tstinf.f $(OBJ)
	gfortran -o tstinf tstinf.f $(OBJ)

# bounded intervals, main program in f90
testfin: testfin.f90 
	gfortran -o testfin testfin.f90 $(OBJ)

# infinite intervals, main program in f90
testinf: testinf.f90 
	gfortran -o testinf testinf.f90 $(OBJ)

#tarfile:
#	tar cvzf quadpackf90.tgz Makefile *.f *.f90

#zipfile:
#	zip -l quadpackf90.zip Makefile *.f *.f90
