#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 testing_nemesis.py

echo 'convert .png files to .gif'

convert -delay 10 'phase_space_map_frame=*_tree_singlestar.png' -loop 0 map_tree_singlestar.gif
convert -delay 10 'phase_space_map_frame=*_tree_singlestar.png' -loop 0 map_tree_binarycluster.gif
convert -delay 10 'phase_space_map_frame=*_tree_singlestar.png' -loop 0 map_tree_singlecluster.gif
convert -delay 10 'phase_space_map_frame=*_Nbody_singlestar.png' -loop 0 map_tree_singlestar.gif
convert -delay 10 'phase_space_map_frame=*_Nbody_singlestar.png' -loop 0 map_tree_binarycluster.gif
convert -delay 10 'phase_space_map_frame=*_Nbody_singlestar.png' -loop 0 map_tree_singlecluster.gif

rm -rf phase_space_map_frame=*.png

mv *.txt ../data
mv *.npy ../data
mv *.png ../figures
mv *.gif ../figures
