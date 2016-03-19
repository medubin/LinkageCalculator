import sys
import csv
import glob

def write_to_file(l): # writes to a new file with a name equal to the given file (minus the 3 character extension) + the mod string
	with open("All_Averages_Metadata_D.csv", "wb") as f:
		writer = csv.writer(f)
		writer.writerows(l)



def create_average_list(): #Calculates
	the_list = [] 
	list_to_file = []
	transposed_list=[]

	for files in glob.glob("*_averages.csv"):
		the_list.append(files[0:-13])
		print files
		first_line =  True
		for row in csv.reader(open(files, 'rU'), delimiter = ','):
			if first_line: 
				first_line = False
			else:
				the_list.append(row[2])

		list_to_file.append(the_list)
		the_list = []

	distance_list = range(1,2001)
	distance_list.insert(0, "Distance")
	list_to_file.insert(0, distance_list)
	transposed_list = zip(*list_to_file)
	write_to_file(transposed_list)












create_average_list()