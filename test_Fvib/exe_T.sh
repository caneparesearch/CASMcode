#!/bin/bash
mkdir chempot_${1}
cd chempot_${1}
cat > monte.json <<EOF
{
    "comment" : "Built from example",
    "debug" : false,
    "ensemble" : "grand_canonical",
    "method" : "metropolis",
    "model" : {
      "formation_energy" : "formation_energy"
    },
    "supercell" : [
      [10, 0, 0],
      [0, 10, 0],
      [0, 0, 10]
    ],
    "data" : {
      "sample_by" : "pass",
      "sample_period" : 1,
      "min_pass" : 1000,
      "max_pass" : 10000,
      "confidence" : 0.95,
      "measurements" : [
        {
          "quantity" : "formation_energy",
		 "precision" : 1e-3
        },
        {
          "quantity" : "potential_energy",
		 "precision" : 1e-3
        },
        {
          "quantity" : "clex_hull_dist(ALL,comp)",
		 "precision" : 1e-3
        },
        {
          "quantity" : "atom_frac"
        },
        {
          "quantity" : "site_frac"
        },
        {
          "quantity" : "comp",
          "precision" : 1e-3
        },
        {
          "quantity" : "comp_n"
        }
      ],
      "storage" : {
        "write_observations" : false,
        "write_trajectory" : false,
        "output_format" : ["csv", "json"]
      }
    },
    "driver" : {
      "mode" : "incremental",
      "motif" : {
        "configname" : "restricted_auto"
      },
      "initial_conditions" : {
        "param_chem_pot" : {
          "a" : ${1},
          "b" : 0
        },
        "temperature" : ${5},
        "tolerance" : 0.001
      },
      "final_conditions" : {
        "param_chem_pot" : {
          "a" : ${1},
          "b" : 0
        },
        "temperature" : ${6},
        "tolerance" : 0.001
      },
      "incremental_conditions" : {
        "param_chem_pot" : {
          "a" : 0,
          "b" : 0
        },
        "temperature" : ${7},
        "tolerance" : 0.001
      }
    }
  }
EOF

echo `date` "Starting job ..." > stdout.txt
echo `date` "Execute Job${2} chempot = ${1} using Slot${3} from T = ${5} to ${6} with dT = ${7}"  >> stdout.txt

casm monte -s monte.json >> stdout.txt
mv stdout.txt stdout_${1}.txt
python ../../Analysis.py
