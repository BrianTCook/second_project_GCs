#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 testing_nemesis.py

echo 'create gifs'
convert -delay 10 'snapshot_Nbody_SingleStar_*.png' -loop 0 Nbody_singlestar.gif
convert -delay 10 'snapshot_tree_SingleStar_*.png' -loop 0 tree_singlestar.gif
convert -delay 10 'snapshot_Nbody_SingleCluster_*.png' -loop 0 Nbody_singlecluster.gif
convert -delay 10 'snapshot_tree_SingleCluster_*.png' -loop 0 tree_singlecluster.gif

#rm -rf snapshot_*.png
#mv *.txt ../data
#mv *.npy ../data
#mv *.csv ../data
#mv *.pdf ../figures
#mv *.png ../figures
#mv *.gif ../figures

