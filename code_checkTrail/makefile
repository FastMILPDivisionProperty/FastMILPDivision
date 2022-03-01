CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++17 -fopenmp -Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

IFLAGS = -I $(GUROBI_HOME)/include/ -I /usr/include/x86_64-linux-gnu/
LFLAGS = -L $(GUROBI_HOME)/lib/ -lgurobi_g++5.2 -lgurobi90 -lgmpxx -lgmp -ltbb

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@

main :  main.o customCallback.o Matrix.o aux_function.o led.o divtrail.o Polynome.o DivTable.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o main main.o customCallback.o Matrix.o aux_function.o led.o divtrail.o Polynome.o DivTable.o $(LFLAGS)

timing :  timing.o customCallback.o Matrix.o aux_function.o skinny64.o divtrail.o Polynome.o DivTable.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o timing timing.o customCallback.o Matrix.o aux_function.o skinny64.o divtrail.o Polynome.o DivTable.o $(LFLAGS)

testMinor :  testMinor.o Matrix.o aux_function.o 
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o testMinor testMinor.o Matrix.o aux_function.o $(LFLAGS)

generateAllMatrixFiles :  generateAllMatrixFiles.o genMatrixFiles.o Matrix.o aux_function.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o generateAllMatrixFiles generateAllMatrixFiles.o genMatrixFiles.o Matrix.o aux_function.o $(LFLAGS)

clean :
	rm -rf *.o
