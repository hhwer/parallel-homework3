
ALL: ex2
CFLAGS = 
FFLAGS = 
CPPFLAGS = 
FPPFLAGS = 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

ex2: ex2.o chkopts
	${CLINKER}-o ex2 ex2.o ${PETSC_LIB}
	${RM} ex2.o
	
