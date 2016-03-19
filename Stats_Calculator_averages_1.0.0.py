import sys
import csv
import glob

def write_to_file(l): # writes to a new file with a name equal to the given file (minus the 3 character extension) + the mod string
	with open("Calculations_from_average_values.csv", "wb") as f:
		writer = csv.writer(f)
		writer.writerows(l)

def calculations_background(l):
	average = 0.0
	counting = 0.0
	for row in l[-500:]:
		#print row
		average = average + float(row[1])
		counting += 1

	return average/counting

def calculation_distance(l,av): # Calculates the distance until the value hits the the 110% of the background
	counter = 0
	av = float(av * 1.1)
	for row in l[1:]:
		if float(row[1]) < av:
			counter += 1
		else:
			counter = 0
		if counter == 3:
			return row[0]
	else:	
		return 3000.0
	
		


def calculate_stats(): #Calculates
	the_list = [] 
	list_to_file = [["File", "Starting_value", "Background_value", "Distance_to_background", "Average n"]]
	for files in glob.glob("*averages.csv"):
		print files
		for row in csv.reader(open(files, 'rU'), delimiter = ','):
			the_list.append(row)

		average_value = calculations_background(the_list)
		print average_value
		list_to_file.append([files[0:-13],the_list[1][1],average_value, calculation_distance(the_list, average_value), calculate_average_n(the_list) ])
		the_list = []

	
	write_to_file(list_to_file)


def calculate_average_n(l):
	n_total = 0.0
	n_counter = 0.0

	for row in l[1:]:
		n_total = n_total + float(row[3])
		n_counter += 1
	return n_total/n_counter










calculate_stats()