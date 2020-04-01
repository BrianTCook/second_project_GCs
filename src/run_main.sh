#!/usr/bin/env bash

#should be in ~/second_projects_GCs/src/

echo 'run script'
python testing_nemesis.py >> testing_nemesis.txt

#echo 'create gifs'
'''
convert -delay 2 'snapshot_Nbody_SingleStar_*.png' -loop 0 Nbody_singlestar.gif
convert -delay 2 'snapshot_Nbody_SingleCluster_*.png' -loop 0 Nbody_singlecluster.gif
convert -delay 2 'snapshot_tree_SingleStar_*.png' -loop 0 tree_singlestar.gif
convert -delay 2 'snapshot_tree_SingleCluster_*.png' -loop 0 tree_singlecluster.gif
convert -delay 2 'snapshot_nemesis_SingleStar_*.png' -loop 0 nemesis_singlestar.gif
convert -delay 2 'snapshot_nemesis_SingleCluster_*.png' -loop 0 nemesis_singlecluster.gif
'''

#rm -rf snapshot_*.png

mv enbid_*.txt ../data/enbid_files
mv times_*.txt ../data/simulation_times
mv tree_*.txt ../data/tree_data
#mv Nbody_*.txt ../data/Nbody_data
mv nemesis_*.txt ../data/nemesis_data
mv *.npy.gz ../data/zipped_files
#mv *.csv ../data
mv *.pdf ../figures
mv *.png ../figures
mv *.gif ../figures

