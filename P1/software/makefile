# Macros de directorios y opciones de compilación

BIN=./bin
SRC=./src
IN=./input
OUT=./output

OPT=-Wall -O2



# Macros

all : compile_all
	./executeAlgorithms.sh

compile_all : $(BIN)/greedy $(BIN)/localSearch $(BIN)/hybrid

examples : exampleGreedy exampleLocalSearch exampleHybrid



# Compilar

$(BIN)/greedy :
	g++ $(OPT) $(SRC)/greedy.cpp -o $(BIN)/greedy

$(BIN)/localSearch :
	g++ $(OPT) $(SRC)/localSearch.cpp -o $(BIN)/localSearch

$(BIN)/hybrid :
	g++ $(OPT) $(SRC)/hybrid.cpp -o $(BIN)/hybrid



# Ejemplos de ejecución

exampleGreedy : $(BIN)/greedy
	$(BIN)/greedy < $(IN)/MDG-a_01_n500_m50.txt
	$(BIN)/greedy < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/greedy < $(IN)/MDG-c_10_n3000_m400.txt

exampleLocalSearch : $(BIN)/localSearch
	$(BIN)/localSearch < $(IN)/MDG-a_01_n500_m50.txt
	$(BIN)/localSearch < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/localSearch < $(IN)/MDG-c_10_n3000_m400.txt

exampleHybrid : $(BIN)/hybrid
	$(BIN)/hybrid < $(IN)/MDG-a_01_n500_m50.txt
	$(BIN)/hybrid < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/hybrid < $(IN)/MDG-c_10_n3000_m400.txt



# Eliminar binarios

clean :
	rm -f $(BIN)/*
	rm -f $(OUT)/*
