CPP	   = g++
INCL       = -I ../Eigen -I../dlib
OPTIONS	   = # -O3
EXECUTABLE = test
OBJECTS	   = test.o Covfunc.o Gppe.o Tool.o Learn.o

${EXECUTABLE}: ${OBJECTS}
	${CPP} -o ${EXECUTABLE} ${OBJECTS}
	
test.o: test.cpp Covfunc.h Gppe.h Tool.h Learn.o
	${CPP} ${INCL} ${OPTIONS} -c test.cpp
	
Covfunc.o: Covfunc.cpp Covfunc.h
	${CPP} ${INCL} ${OPTIONS} -c Covfunc.cpp

Gppe.o: Gppe.cpp Gppe.h
	${CPP} ${INCL} ${OPTIONS} -c Gppe.cpp
	
Tool.o: Tool.cpp Tool.h
	${CPP} ${INCL} ${OPTIONS} -c Tool.cpp
	
Learn.o: Learn.cpp Learn.h
	${CPP} ${INCL} ${OPTIONS} -c Learn.cpp

clean:
	rm -rf ${OBJECTS} ${EXECUTABLE}
