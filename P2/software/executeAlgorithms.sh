#!/bin/bash


for algorithm in `ls bin`
do
    rm -f output/${algorithm}.dat
    touch output/${algorithm}.dat
    echo "ALGORITHM: $algorithm"
    for input_file in `ls input`
    do
        echo -e "\t$input_file"
        bin/${algorithm} < input/${input_file} >> output/${algorithm}.dat
    done
done

exit 0
