#!/bin/bash
filename="fit.json"
rm -rf ../.casm/tmp
rm -rf ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.__tmp
rm ${filename%.*}_*

casm-learn -s $filename
casm-learn -s $filename --checkhull > fit_log.txt
casm-learn -s $filename --select 0
casm-learn -s $filename --hall --indiv 0 --format json > ${filename%.*}-eci.json
casm query -k  'comp(a)' 'formation_energy' 'clex(formation_energy)' 'hull_dist(ALL,comp)'  'clex_hull_dist(ALL,comp)' 'comp_n(Ta)'  -c casm_learn_input   -o data.dat

cp data.dat ${filename%.*}_fit.dat
cp fit_fit.dat hull.dat
cp ../cluster_expansions/clex.formation_energy/calctype.default/ref.default/bset.default/eci.__tmp/eci.json .
echo ${filename%.*}
