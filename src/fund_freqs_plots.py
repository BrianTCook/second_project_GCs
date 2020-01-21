from amuse.io import read_set_from_file

directory = '/home/brian/Desktop/second_project_GCs/data/'
filename = 'Brute_data.hdf5'
bodies = read_set_from_file(directory+filename, "hdf5")
print(bodies)
