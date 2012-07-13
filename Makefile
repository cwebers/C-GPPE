CPP	   = g++
INCL       = -I ../Eigen -I../dlib
OPTIONS	   =  -fopenmp
OPTIONS   += -O3
EXECUTABLE = test
OBJECTS	   = test.o Covfunc.o CGppe.o Tool.o CLearner.o

${EXECUTABLE}: ${OBJECTS}
	${CPP} -o ${EXECUTABLE} ${OBJECTS}
	
test.o: test.cpp Covfunc.h CGppe.h Tool.h CLearner.o
	${CPP} ${INCL} ${OPTIONS} -c test.cpp
	
Covfunc.o: Covfunc.cpp Covfunc.h
	${CPP} ${INCL} ${OPTIONS} -c Covfunc.cpp

CGppe.o: CGppe.cpp CGppe.h
	${CPP} ${INCL} ${OPTIONS} -c CGppe.cpp
	
Tool.o: Tool.cpp Tool.h
	${CPP} ${INCL} ${OPTIONS} -c Tool.cpp
	
CLearner.o: CLearner.cpp CLearner.h
	${CPP} ${INCL} ${OPTIONS} -c CLearner.cpp

clean:
	rm -rf ${OBJECTS} ${EXECUTABLE}
