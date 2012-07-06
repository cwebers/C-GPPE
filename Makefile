COMP	=g++ -I ../Eigen -I ../ #-O3
OPTIONS	=-c
EXECUTABLE	=test
OBJETS		=test.o Covfunc.o Gppe.o Tool.o Learn.o

${EXECUTABLE}:${OBJETS}
	${COMP} -o ${EXECUTABLE} ${OBJETS}
	
test.o:test.cpp Covfunc.h Gppe.h Tool.h Learn.o
	${COMP} ${OPTIONS} test.cpp
	
Covfunc.o :Covfunc.cpp Covfunc.h
	${COMP} ${OPTIONS} Covfunc.cpp

Gppe.o :Gppe.cpp Gppe.h
	${COMP} ${OPTIONS} Gppe.cpp
	
Tool.o :Tool.cpp Tool.h
	${COMP} ${OPTIONS} Tool.cpp
	
Learn.o :Learn.cpp Learn.h
	${COMP} ${OPTIONS} Learn.cpp

clean:
	rm -rf *.o test
