#should be in ~/second_projects_GCs/src/

echo 'run script'
python3 testing_nemesis.py

echo 'no gifs, just run the simulation and get the energy/\av{w} plots'

mv *.txt ../data
mv *.npy ../data
mv *.png ../figures
mv *.gif ../figures
