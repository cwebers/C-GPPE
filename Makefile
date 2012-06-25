COMP	=g++ -I ../Eigen
OPTIONS	=-c
EXECUTABLE	=test
OBJETS		=testcov.o Covfunc.o Gppe.o Tool.o

${EXECUTABLE}:${OBJETS}
	${COMP} -o ${EXECUTABLE} ${OBJETS}
	
testcov.o:testcov.cpp Covfunc.h Gppe.h Tool.h
	${COMP} ${OPTIONS} testcov.cpp
	
Covfunc.o :Covfunc.cpp Covfunc.h
	${COMP} ${OPTIONS} Covfunc.cpp

Gppe.o :Gppe.cpp Gppe.h
	${COMP} ${OPTIONS} Gppe.cpp
	
Tool.o :Tool.cpp Tool.h
	${COMP} ${OPTIONS} Tool.cpp