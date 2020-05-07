#!/usr/bin/env bash

#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 main.py

mv enbid_*.txt ../data/enbid_files
mv times_*.txt ../data/simulation_times
mv tree_*.txt ../data/tree_data/
mv *.pdf ../figures
mv *.png ../figures
mv *.gif ../figures

