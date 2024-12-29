CC = g++ -std=c++17
SRC = src
OBJ = obj

all: $(OBJ)/fit.o

$(OBJ)/fit.o: $(SRC)/fit.cpp $(SRC)/matrix_calc.cpp
	$(CC) $(SRC)/fit.cpp $(SRC)/matrix_calc.cpp -o $(OBJ)/fit.o

clean:
	rm $(OBJ)/*