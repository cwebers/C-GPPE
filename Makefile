COMP	=g++ -arch x86_64 -I ../Eigen
OPTIONS	=-c
EXECUTABLE	=test
OBJETS		=testcov.o Covfunc.o Gppe.o

${EXECUTABLE}:${OBJETS}
	${COMP} -o ${EXECUTABLE} ${OBJETS}
	
testcov.o:testcov.cpp Covfunc.h Gppe.h
	${COMP} ${OPTIONS} testcov.cpp
	
Covfunc.o :Covfunc.cpp Covfunc.h
	${COMP} ${OPTIONS} Covfunc.cpp

Gppe.o :Gppe.cpp Gppe.h
	${COMP} ${OPTIONS} Gppe.cpp
	
	
	
	#-arch x86_64