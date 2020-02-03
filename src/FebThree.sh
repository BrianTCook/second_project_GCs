#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 testing_nemesis.py

echo 'convert .png files to .gif'

convert -delay 10 'phase_space_map_frame=*_SingleStar_tree.png' -loop 0 map_tree_singlestar.gif
convert -delay 10 'phase_space_map_frame=*_BinaryCluster_tree.png' -loop 0 map_tree_binarycluster.gif
convert -delay 10 'phase_space_map_frame=*_SingleCluster_tree.png' -loop 0 map_tree_singlecluster.gif
convert -delay 10 'phase_space_map_frame=*_SingleStar_Nbody.png' -loop 0 map_Nbody_singlestar.gif
convert -delay 10 'phase_space_map_frame=*_BinaryCluster_Nbody.png' -loop 0 map_Nbody_binarycluster.gif
convert -delay 10 'phase_space_map_frame=*_SingleCluster_Nbody.png' -loop 0 map_Nbody_singlecluster.gif

rm -rf phase_space_map_frame=*.png

mv *.txt ../data
mv *.npy ../data
mv *.png ../figures
mv *.gif ../figures
