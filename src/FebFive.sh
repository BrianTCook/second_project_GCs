#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 testing_nemesis.py

#echo 'create BinaryCluster gifs'
#convert -delay 10 'snapshot_Nbody_*.png' -loop Nbody.gif
#convert -delay 10 'snapshot_tree_*.png' -loop tree.gif

mv *.txt ../data
mv *.npy ../data
mv *.pdf ../figures
mv *.png ../figures
mv *.gif ../figures
