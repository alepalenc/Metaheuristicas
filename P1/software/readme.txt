////////// CONTENIDO DEL DIRECTORIO SOFTWARE //////////

- makefile: fichero para compilar y ejecutar los algoritmos
- executeAlgorithms.sh: script para ejecutar los algoritmos
- src: códigos fuente
- bin: binarios
- input: instancias del problema (inicialmente vacío)
- output: resultados obtenidos al ejecutar el script


////////// COMANDOS DEL FICHERO MAKEFILE //////////

- all: compila los tres algoritmos y ejecuta el script executeAlgorithms.sh
- bin/<algoritmo>: compila el código fuente con el algoritmo especificado (greedy, localSearch, hybrid)
- example<algoritmo>: compila el algoritmo especificado (greedy, localSearch, hybrid) y lo ejecuta con tres 
                      instancias del problema con distintos tamaños, mostrando los resultados por la terminal
- examples: compila los tres algoritmos y realiza sus ejemplos de ejecución
- compile_all: compila los tres algoritmos
- clean: elimina el contenido de los directorios bin y output


////////// PASOS PARA REPETIR LOS EXPERIMENTOS //////////

1. Introducir los ficheros .txt con las instancias en el directorio input
2. Desde una terminal, escribir: make all
