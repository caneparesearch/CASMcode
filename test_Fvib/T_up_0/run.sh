#!/bin/bash
rm -rf conditions.* stdout.txt results*

casm monte -s monte.json > stdout.txt
