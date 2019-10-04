# Variables:
CXX= g++
CXXFLAGS= -g -std=c++11 -Wall -I/usr/include -I$(INC) -O3 -march=native
SRC= src/
INC= include/
OBJ= obj/
BIN= bin/

# Reglas para construir el programa principal:
all: $(BIN)QAP_AG $(BIN)QAP_AM

$(BIN)QAP_AG: $(OBJ)PBasedALG.o $(OBJ)random.o $(OBJ)mainAG.o
	$(CXX) $(CXXFLAGS) -o $(BIN)QAP_AG $^

$(BIN)QAP_AM: $(OBJ)PBasedALG.o $(OBJ)random.o $(OBJ)mainAM.o
	$(CXX) $(CXXFLAGS) -o $(BIN)QAP_AM $^

$(OBJ)PBasedALG.o: $(SRC)PBasedALG.cpp $(INC)PBasedALG.h
	$(CXX) $(CXXFLAGS) -o $(OBJ)PBasedALG.o -c $(SRC)PBasedALG.cpp

$(OBJ)random.o: $(SRC)random.cpp $(INC)random.h
	$(CXX) $(CXXFLAGS) -o $(OBJ)random.o -c $(SRC)random.cpp

$(OBJ)mainAG.o: $(SRC)mainAG.cpp
	$(CXX) $(CXXFLAGS) -o $(OBJ)mainAG.o -c $(SRC)mainAG.cpp

$(OBJ)mainAM.o: $(SRC)mainAM.cpp
	$(CXX) $(CXXFLAGS) -o $(OBJ)mainAM.o -c $(SRC)mainAM.cpp  


# Limpiar
clean:
	rm $(OBJ)*.o
	rm $(BIN)QAP_AG
	rm $(BIN)QAP_AM

cleanobj:
	rm $(OBJ)*.o
