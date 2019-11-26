#'should be run in ~/amuse directory'
#'should put into ~/second_project_GCs for safekeeping'

echo 'run script'
#rm -rf logfile_main.txt
./amuse.sh ~/Desktop/second_project_GCs/src/main.py #> logfile_main.txt

echo 'convert collection of .png files into .gif'
convert -delay 10 'frame_*.png' -loop 0 GCs_in_background.gif
rm -rf frame_*.png

mv *.png ~/Desktop/second_project_GCs/figures/
mv GCs_in_background.gif ~/Desktop/second_project_GCs/figures/
cp main_GCs.sh ~/Desktop/second_project_GCs/src/
