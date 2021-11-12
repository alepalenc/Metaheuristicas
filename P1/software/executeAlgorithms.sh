#!/bin/bash

output_file_greedy="output/greedy.dat"
output_file_local_search="output/localSearch.dat"
output_file_hybrid="output/hybrid.dat"

rm -f ${output_file_greedy}
rm -f ${output_file_local_search}
rm -f ${output_file_hybrid}
touch ${output_file_greedy}
touch ${output_file_local_search}
touch ${output_file_hybrid}

echo "ALGORITHM: GREEDY"
for input_file in `ls input`
do
    echo -e "\t$input_file"
    bin/greedy < input/${input_file} >> ${output_file_greedy}
done

echo "ALGORITHM: LOCAL SEARCH"
for input_file in `ls input`
do
    echo -e "\t$input_file"
    bin/localSearch < input/${input_file} >> ${output_file_local_search}
done

echo "ALGORITHM: HYBRID"
for input_file in `ls input`
do
    echo -e "\t$input_file"
    bin/hybrid < input/${input_file} >> ${output_file_hybrid}
done

exit 0
