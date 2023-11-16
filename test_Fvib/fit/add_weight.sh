#!/bin/bash
FILE1='casm_learn_input'
while IFS= read -r line
do
    awk '$27= "1.0"' casm_learn_input > tmp.txt && mv tmp.txt file.txt
done <"$FILE1"
