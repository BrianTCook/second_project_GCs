#'should be run in ~/amuse directory'
#'should put into ~/second_project_GCs for safekeeping'

echo 'run script'
./amuse.sh ~/Desktop/second_project_GCs/src/main.py 

echo 'convert collection of .png files into .gif'
convert -delay 10 'frame_*_Nemesis.png' -loop 0 GCs_in_background_nem.gif
convert -delay 10 'frame_*_Brute.png' -loop 0 GCs_in_background_bru.gif
rm -rf frame_*.png

mv *.png ~/Desktop/second_project_GCs/figures/
mv *.gif ~/Desktop/second_project_GCs/figures/
mv *.hdf5 ~/Desktop/second_project_GCs/src/
cp main_GCs.sh ~/Desktop/second_project_GCs/src/
