# Macros de directorios y opciones de compilación

BIN=./bin
SRC=./src
IN=./input
OUT=./output

OPT=-Wall -O2



# Macros

all : compile_all
	./executeAlgorithms.sh

compile_all : $(BIN)/ES $(BIN)/BMB $(BIN)/ILS $(BIN)/ILS_ES



# Compilar

$(BIN)/ES :
	g++ $(OPT) $(SRC)/ES.cpp -o $(BIN)/ES

$(BIN)/BMB :
	g++ $(OPT) $(SRC)/BMB.cpp -o $(BIN)/BMB

$(BIN)/ILS :
	g++ $(OPT) $(SRC)/ILS.cpp -o $(BIN)/ILS

$(BIN)/ILS_ES :
	g++ $(OPT) $(SRC)/ILS_ES.cpp -o $(BIN)/ILS_ES



# Ejemplos de ejecución

exampleES : $(BIN)/ES
	$(BIN)/ES < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/ES < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/ES < $(IN)/MDG-c_10_n3000_m400.txt

exampleBMB : $(BIN)/BMB
	$(BIN)/BMB < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/BMB < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/BMB < $(IN)/MDG-c_10_n3000_m400.txt

exampleILS : $(BIN)/ILS
	$(BIN)/ILS < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/ILS < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/ILS < $(IN)/MDG-c_10_n3000_m400.txt

exampleILS_ES : $(BIN)/ILS_ES
	$(BIN)/ILS_ES < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/ILS_ES < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/ILS_ES < $(IN)/MDG-c_10_n3000_m400.txt



# Eliminar binarios

clean :
	rm -f $(BIN)/*
	rm -f $(OUT)/*
