# Macros de directorios y opciones de compilación

BIN=./bin
SRC=./src
IN=./input
OUT=./output

OPT=-Wall -O2



# Macros

all : compile_all
	./executeAlgorithms.sh

compile_all : $(BIN)/AGGuniform $(BIN)/AGGposition $(BIN)/AGEuniform $(BIN)/AGEposition $(BIN)/AM1 $(BIN)/AM2 $(BIN)/AM3



# Compilar

$(BIN)/AGGuniform :
	g++ $(OPT) $(SRC)/AGGuniform.cpp -o $(BIN)/AGGuniform

$(BIN)/AGGposition :
	g++ $(OPT) $(SRC)/AGGposition.cpp -o $(BIN)/AGGposition

$(BIN)/AGEuniform :
	g++ $(OPT) $(SRC)/AGEuniform.cpp -o $(BIN)/AGEuniform

$(BIN)/AGEposition :
	g++ $(OPT) $(SRC)/AGEposition.cpp -o $(BIN)/AGEposition

$(BIN)/AM1 :
	g++ $(OPT) $(SRC)/AM1.cpp -o $(BIN)/AM1

$(BIN)/AM2 :
	g++ $(OPT) $(SRC)/AM2.cpp -o $(BIN)/AM2

$(BIN)/AM3 :
	g++ $(OPT) $(SRC)/AM3.cpp -o $(BIN)/AM3



# Ejemplos de ejecución

exampleAGGuniform : $(BIN)/AGGuniform
	$(BIN)/AGGuniform < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AGGuniform < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AGGuniform < $(IN)/MDG-c_10_n3000_m400.txt

exampleAGGposition : $(BIN)/AGGposition
	$(BIN)/AGGposition < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AGGposition < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AGGposition < $(IN)/MDG-c_10_n3000_m400.txt

exampleAGEuniform : $(BIN)/AGEuniform
	$(BIN)/AGEuniform < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AGEuniform < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AGEuniform < $(IN)/MDG-c_10_n3000_m400.txt

exampleAGEposition : $(BIN)/AGEposition
	$(BIN)/AGEposition < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AGEposition < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AGEposition < $(IN)/MDG-c_10_n3000_m400.txt

exampleAM1 : $(BIN)/AM1
	$(BIN)/AM1 < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AM1 < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AM1 < $(IN)/MDG-c_10_n3000_m400.txt

exampleAM2 : $(BIN)/AM2
	$(BIN)/AM2 < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AM2 < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AM2 < $(IN)/MDG-c_10_n3000_m400.txt

exampleAM3 : $(BIN)/AM3
	$(BIN)/AM3 < $(IN)/MDG-a_10_n500_m50.txt
	$(BIN)/AM3 < $(IN)/MDG-b_21_n2000_m200.txt
	$(BIN)/AM3 < $(IN)/MDG-c_10_n3000_m400.txt



# Eliminar binarios

clean :
	rm -f $(BIN)/*
	rm -f $(OUT)/*
