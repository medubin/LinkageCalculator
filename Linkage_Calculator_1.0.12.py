import sys
import csv
from operator import itemgetter
import re
import math


"""
Linkage_Calculator Version 1.0.12

Arguments:
p=	'all'		default
	'trim'
	'ldcalc'
	'avcalc'
	'regional'    calc:average
	'regionald'   calc:average   uses Dprime
	'regionalc'   calc:average:comparison   compares to Mcguagh
	'avcaclcc'     calc:comparison
	Sets the proccess to run
d=	2000		default
	takes a value for maximum comparison distance
s 	removes singletons		default = false
f 	removes snps with less than 4/4 nucleotide combinations default = false
n= 	5 			default
	number of Ns allowed per loci default
lr=	0
	number of lines removed off the end during trim. Doesn't use them for calculation either

New in version 1.0.12
Added avcalcc comparing calcuation data to Mcguagh data
Added avecalc10kb comparing calculation data to premade bin sizes (currently 50kb, in spite of the 10kb name)

New in version 1.0.11
Added comparing regional data to Mcguagh

New in version 1.0.10
Adding a new method: regionald. Same as regional, but calculates the value for Dprime. 

New in version 1.0.9
The line removal function was broken. I fixed it I will need to rerun EVERTYTHING on 2/13/15

New in version 1.0.8
Fixed the calculate regions method. Previously it was skipping some of the data.

New in version 1.0.7
Fixed a with AvgBP not calculating properly.
Modifed regional calculations to recalculate AvgBP instead of taking the incorrect one in the files

New in version 1.0.6
Line by line file opening for trim and avcalc. Still can't figure out how to do it for ldcalc



New in version 1.0.5
Sorting doesn't work for large files. Created new function to deal with large files
Fixed a bug in the original average_calculator function which caused it to  not append the file with the last average (the one equal to comparison_distance)



New in version 1.0.4
Optimization for larger files:
	Fixed an issue where the ldcalc iterator moves through each possible comparison (even though it doesn't make comparsions past distance)
	Set the second row in the linkage_disequilibrium_iterator to have a range end at bp1 + comparison_distance. This significantly sped up processing of large files

New in version 1.0.3
Changed format of file input. Now requires first line to be base number
Added an argument that removes the last x lines during trim
Added Dprime calculations into average

New in version 1.0.2
Adding more comments to code
Fixed the all process, in which the calculate average function was not called
Added support for tab delimited txt files

New in version 1.0.1
Added arguments [p, d,s,f, n].
Added the ability to run all 3 basic process in a single run

New in version 1.0.0
Added version tracking. Now functional in all basic processes




"""
#global default arguments 	Defaults
process_type = 				'all' # Which of the 3 functions it will run. The default is all three of them (trim -> ldcalc -> avcalc)
comparison_distance = 		2000 # How are away the maximum distance for bp comparison. Used in the ldcalc process. Default is 2000 
singleton_removal = 		False # Whether singletons (loci that only have a single instance of polymorphism) are removed during the trim process
full_only = 				False # Removed all comparisons in which not all 4 nucleotide combinations are present. Used in ldcalc process
ns_allowed =				5 # The number of Ns (or blanks or not ATCG) allowed. Used during the trim process
lines_removed = 			0 # Removes lines from the end of the file during the trim function.


def delimiter_uncoverer():
	if str(sys.argv[1])[-3:] == 'csv':
		return ','
	elif str(sys.argv[1])[-3:] == 'txt':
		return '\t'

def file_input(): # returns the file
	return (open(sys.argv[1], 'rU')) # sets new_file equal to the inputed file

def modified_file_input(modstring): #used exclusively for the all processes
	file_name = str(sys.argv[1])
	file_name = file_name[0:-4] + modstring #Since the all process sequentially makes files and opens them this is required
	return (open(file_name, 'rU'))


def print_file(f): #prints the file as a series of rows
	try:
	    reader = csv.reader(f)  # creates the reader object
	    for row in reader:   # iterates the rows of the file in orders
	        print row    # prints each row
	finally:
		f.close()      # closing

def print_dict(d): # prints a dictionary alternating keys and values line per line
	for keys,values in d.items():
		print (keys)
		print (values)

def print_list(l): # prints a list with each sublist on a new line
	for row in l:
		print row		


def read_to_list(l, file_extension = str(sys.argv[1])[-3:]): # reads the file and returns a list containing rows
	
	li = [] # sets up a new list to fill with rows
	try:
		if file_extension == 'csv': 
			reader = csv.reader(l, delimiter = ',') #reads the file as comma deliminated
			
		elif file_extension == 'txt': 
			reader = csv.reader(l, delimiter = '\t') #reads the file as tab deliminated
		for row in reader:
				li.append(row)
	finally:
		l.close()  # closing
	return li


def write_to_file(l, modstring): # writes to a new file with a name equal to the given file (minus the 3 character extension) + the mod string
	with open(str(sys.argv[1])[0:-4] + modstring, "wb") as f:
		writer = csv.writer(f)
		writer.writerows(l)


def append_to_file(l, modstring): # appends to an existing file with a name equal to the given file (minus the 3 character extension) + the mod string
	with open(str(sys.argv[1])[0:-4] + modstring, "a") as f: # The "a" appends instead of writes over
		writer = csv.writer(f)
		writer.writerows(l)

	
def trim_non_poly(): # Removes non informative markers from the list l
	global ns_allowed # how many Ns are allowed per loci
	global singleton_removal # boolean of whether singletons (loci with a single version of a polymorphic snp) are removed
	global lines_removed
	file_name = file_input() # sets new_file equal to the inputed file
	nl = []
	heading = ''.join(file_name.readline())
	headlist = [heading.split(delimiter_uncoverer())] 
	headlist =[headlist[0][:-lines_removed]]# All of this is just to remove a few stupid lines from the header
	write_to_file(headlist, '_trimmed.csv') # creates new list
	x = 0 # sets up the counter, to keep track of where I am in the list (is there a better way to do this?)
	counter = 0
	for row in csv.reader(file_name, delimiter = delimiter_uncoverer()): # goes through the rest of the lines 1 at a time: # iterates through each line in the original list
		if lines_removed !=0: row = row[:-lines_removed]
		A, T, G, C, N, blank = 0, 0, 0, 0, 0,0 # for each new line set the count of bp at 0
		for bp in row[1:]: # iterates through each individual at a loci
			if bp == 'A': A += 1 # next 6 lines: adds up each possible bp for the individuals
			elif bp == 'T': T += 1
			elif bp == 'G': G += 1
			elif bp == 'C': C += 1
			elif bp == 'N': N += 1
			else: blank += 1
		totallength = T + G + A + C # calculates number of non N/blank individuals
		if all([(N+ blank) <= ns_allowed, A != totallength, T != totallength, G != totallength, C != totallength, trim_too_many_nucleotides(A, T, G, C)]): # checks if N is 5 or more of the individuals, and if any of the nucleotides are 100% fixed and if there are more than 2 nucleotides
			if singleton_removal == False or all([singleton_removal, A != 1, T != 1, G != 1, C != 1]): # checks if singleton removal is true, if it is, checks for singletons and denies appending if it finds them
				nl.append(row) #if the marker is informative, pass it to the new list
				counter +=1
				if counter % 10000 == 0: 
					print "trimmed row " + str(counter)
					append_to_file(nl, '_trimmed.csv')
					nl =[]
		x = x + 1 # add to the counter

	append_to_file(nl, '_trimmed.csv')
	

	


def trim_too_many_nucleotides(A,T,G,C): # returns True if there are only 2 nucleotide types
	counter = 0 # counts number of unique nucleotides
	if A > 0: counter += 1 # adds to the counter if that nucleotide is present
	if T > 0: counter += 1
	if G > 0: counter += 1
	if C > 0: counter += 1
	if counter > 2: # if there are more than 2 nucleotide types, return false
		return False
	else:
		return True # otherwise return true


def trim_process(): # Runs the rest of the functions for removing uninformative snps
	print sys.argv[1]
	trim_non_poly()

def linkage_disequilibrium_iterator(l):
	global comparison_distance
	file_list = []
	write_to_file([['BP1','BP2','Distance', 'AvgBP', 'D','Dmax',"Dprime",'r2','n']], "_calculated.csv")
	list_header = l[0] # saves the header for later
	bp1 = 0 # gets assigned the location of the first row being examined
	bp2 = 0 # gets assigned the location of the second row being examined
	del l[0] # deletes header off of list	
	row_counter = 0 # keeps track of what the first row is
	for row in l: # goes through each row of the list
		bp1 = int(row[0]) # assigns bp1 to the current first row
		row_counter += 1
		if row_counter % 10000 == 0: #Every 10000 rows
			print "on row " + str(row_counter) + " of " + str(len(l)) # Every 1000 rows prints what row we are on
			append_to_file(file_list, "_calculated.csv") # Every 1000 rows appends the rows we have calculated to the data file
			file_list = [] # Every 1000 rows epmties the working list
		for secondrow in l[(row_counter):(row_counter + comparison_distance) ]: # goes through all the remaining rows of the list, starting at row + 1
			bp2 = int(secondrow[0]) # assigns bp2 to the current second row
			#if bp2 - bp1 > comparison_distance: print "comparing " + str(bp1) + ' and ' + str(bp2)
			if bp2 - bp1 > comparison_distance: break #if the comparison loci is further than the comparison distance allows, end all comparisons to bp1
			AA, AG, AT, AC = 0, 0, 0, 0 # mass assignment of all possible nucleotide combinations
			GA, GG, GT, GC = 0, 0, 0, 0
			TA, TG, TT, TC = 0, 0, 0, 0
			CA, CG, CT, CC = 0, 0, 0, 0
			for bp in range(1,len(row)): # iterates through all bp of both rows
				if row[bp] == 'A': # figures out, and adds to the counter for, each pair of nucleotides
					if secondrow[bp] == 'A': AA += 1 
					elif secondrow[bp] == 'G': AG += 1 
					elif secondrow[bp] == 'T': AT += 1 
					elif secondrow[bp] == 'C': AC += 1 
				elif row[bp] == 'G': 
					if secondrow[bp] == 'A': GA += 1
					elif secondrow[bp] == 'G': GG += 1
					elif secondrow[bp] == 'T': GT += 1
					elif secondrow[bp] == 'C': GC += 1
				elif row[bp] == 'T':
					if secondrow[bp] == 'A': TA += 1
					elif secondrow[bp] == 'G': TG += 1
					elif secondrow[bp] == 'T': TT += 1
					elif secondrow[bp] == 'C': TC += 1
				elif row[bp] == 'C':
					if secondrow[bp] == 'A': CA += 1
					elif secondrow[bp] == 'G': CG += 1
					elif secondrow[bp] == 'T': CT += 1
					elif secondrow[bp] == 'C': CC += 1
		
			temp_list = calculator_iterator(AA, AG, AT, AC, GA, GG, GT, GC, TA, TG, TT, TC, CA, CG, CT, CC, bp1,bp2) # builds a temporary list to hold these values
			if temp_list: # checks to make sure the list is not empty (which would happen if the calculator_iterator detected loss of polymporphism due to Ns)
				file_list.append(temp_list) # adds the calculated values to the full list
	append_to_file(file_list, "_calculated.csv") #at the end there will be some values left in the file list, so we dumped the rest into the file


		

def calculator_iterator(AA, AG, AT, AC, GA, GG, GT, GC, TA, TG, TT, TC, CA, CG, CT, CC, bp1, bp2): #takes a combination of two loci genotypes and passes them to the various statistical calculators
	global full_only
	xx,xy,yx,yy = 0.0,0.0,0.0,0.0 # the generic ints to hold the 4 values used in calculating D. They start at 0, and stay at 0 unless a value is added to them
	total_value = 0.0
	nuc_list = [['AA',AA],[ 'AG',AG],[ 'AT',AT],[ 'AC', AC],[ 'GA',GA],[ 'GG',GG],[ 'GT',GT],['GC',GC],[ 'TA',TA],['TG', TG],[ 'TT',TT],['TC', TC],['CA', CA],['CG', CG],['CT', CT],['CC',CC]] # builds a list of lists, containing all possible Genotypes and their values
	real_list = [] # this is the list that will contan all the values that actually have a non zero number of items
	for nucpair in nuc_list: # iterates through all the nucleotides in the list
		if nucpair[1] > 0: real_list.append(nucpair) #if the nucleotide count is non zero, add it to the new list
	
	if len(real_list) > 1: #checks to make sure there are at least 2 values in real list (the presence of Ns can cause only 1, which leads to an uninformative d value of 0)
		if full_only == False or all([full_only, len(real_list) >= 4]): # if full only, requires 4 nucleotides instead just 1	
			xx = float(real_list[0][1]) # first sets xx to be the first value in the list
			for nuc in real_list[1:]: # compares xx to others
				if nuc[0][0] == real_list[0][0][0]: xy = float(nuc[1]) #if the first nucleotide of next in list matches xx, it's xy 
				elif nuc[0][1] == real_list[0][0][1]: yx = float(nuc[1]) #else if the second nucleotide of the next in list matches xx, it's yx
				else: yy = float(nuc[1]) #if none of those match, it's yy


			if ((xx+xy) * (yx+yy) * (xx+yx) * (yy+xy)) != 0: #preliminary test of ratios. Makes sure we didn't lose polymorphism at one of the sites
				# All calculations are done here
				total_value = float(xx + xy + yx + yy) #total number of nucletide pairs (the n)
				xx,xy,yx,yy = float(xx/total_value), float(xy/total_value), float(yx/total_value), float(yy/total_value) # Calculates proportions for the pairs
				d =  d_calculator(xx,xy,yx,yy) # calculates d
				dmax = d_max_calculator(xx,xy,yx,yy, d) # calculates d max
				r = r_calculator(xx,xy,yx,yy, d) # calculates r2
				dprime = d/dmax # calculates dprime
				distance = bp2-bp1 # calculates distance
				avgbp = float((bp1+bp2)/2.0) #calculates the midpoint of the bps
		

				return [bp1,bp2, distance, avgbp, abs(d), dmax, abs(dprime), r, total_value] # returns all of the calculated values to the ld interator

		


def d_calculator(xx,xy,yx,yy): # d calculator
	d = (xx*yy) - (xy*yx) # calculates d based on xx*yy - xy*yx = d

	return d

def d_max_calculator(xx,xy,yx,yy, d): #if d > 0: min[x1y2, y1x2] if d < 0: min[x1x2, y1y2]
	x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0 # initializes the proportions of x,y for the first bp and second bp
	dmax1, dmax2 = 0.0, 0.0 # these will contain the two possible dmaxes
	x1 = xx + xy # proportions are calculated
	y1 = yx + yy
	x2 = xx + yx
	y2 = yy + xy

	if d > 0: # if d is +
		dmax1 = x1 * y2 # the first dmax
		dmax2 = y1 * x2 # the second dmx
		return min([dmax1, dmax2]) # returns the smaller of the two
	elif d < 0: # if d is -
		dmax1 = x1 * x2 # the first dmax
		dmax2 = y1 * y2 # the second dmax
		return min([dmax1, dmax2]) # returns the smaller of the two
	else: return 1 # if d is 0
		
def r_calculator(xx,xy,yx,yy, d): #Calculates r^2 = (d/(x1*x2*y1*y2)^.5)^2
	x1, y1, x2, y2 = 0.0, 0.0, 0.0, 0.0 # initializes the proportions of x,y for the first bp and second bp
	x1 = (xx + xy) # proportions are calculated
	y1 = (yx + yy)
	x2 = xx + yx
	y2 = yy + xy
	r = (d)/((x1*y1*x2*y2)**.5)
	r = r**2
	return r


def ld_process(): # highest level in the ld calculator process
	datafile = read_to_list(file_input()) 
	linkage_disequilibrium_iterator(datafile)


def calculate_average_per_distance_process(): #used to calculate the average r per distance measured
	average_calculator_big_data()


def average_calculator_big_data(modstring = '.csv'): # Replaces the average calculator. Reads file 1 line at a time instead of all at once.
	counterizer = 0 # Counts how many rows of data we have done
	av_list = [] # creates a list to store averages
	av_dict = {} # a dictionary to hold the values as they are being added
	file_name = modified_file_input(modstring) # sets new_file equal to the inputed file
	file_name.readline() # reads the header line
	for row in csv.reader(file_name, delimiter = ','): # goes through the rest of the lines 1 at a time
		if row[2] in av_dict: # checks the distance metric in calculation file to see if it already exists in the dictionary if it does:
			av_dict[row[2]][0] += 1.0 #add 1 to the counter
			av_dict[row[2]][1] += float(row[6]) #add the value of Dprime to the running Dprime total
			av_dict[row[2]][2] += float(row[7]) #add the value of r2 to the running r2 total
		else: #if it isn't in the dictionary
			av_dict[row[2]] = [1.0 ,float(row[6]), float(row[7])] #create a new dictionary key with the values of the 1, dprime, and r2
		counterizer += 1 # counterizer keeps track of how far we've gone in the program
		if counterizer % 10000 == 0: print "on row " + str(counterizer)  #lets the user know we've made pogress
	for key in av_dict: #finally, adds the dictionary into a list
		av_list.append([key, (float(av_dict[key][2])/float(av_dict[key][0])), (float(av_dict[key][1])/float(av_dict[key][0])),  float(av_dict[key][0]) ])


	av_list = sorted(av_list, key=lambda av_list: int(av_list[0])) #sorts the list, shouldn't be a problem. List will be 2000
	av_list.insert(0,['Distance','Average r2','Average Dprime', 'n'] ) # insert the header
	write_to_file(av_list, "_averages.csv") #write the list to a file






def average_calculator(l):
	del l[0] #removes the header
	av_list = [['Distance','Average r2','Average Dprime', 'n'],] # creates a list to store averages

	l = sorted(l, key=lambda l: int(l[2])) # sorts the list by distance 
	current_r_total = 0.0 # calculates the additive total for all r2 in that distance (which will later be divided by the number of items in this)
	current_d_total = 0.0 # calculates the additive total for all dprime in that distance (which will later be divided by the number of items in this)
	current_count = 0.0 # count of number of instances of this distance
	current_distance = int(l[0][2]) # starting distance. The first distance in
	for row in l:
		if int(row[2]) == current_distance:
			current_r_total = current_r_total + float(row[7])
			current_d_total = current_d_total + float(row[6])
			current_count += 1
		else:
			av_list.append([current_distance, (current_r_total/current_count),(current_d_total/current_count), current_count])
			current_r_total = float(row[7])
			current_d_total = float(row[6])
			current_count = 1.0
			current_distance = int(row[2])
	av_list.append([current_distance, (current_r_total/current_count),(current_d_total/current_count), current_count])
	
	write_to_file(av_list, "_averages.csv")


def iterate_through_arguments(): # is the first thing called, iterates through the extra arguments (not the program or the data file) and sets the global parameters equal to the arguments given
	global process_type # defines which calculations the program runs
	global comparison_distance # defines the distance which
	global ns_allowed # defines number of Ns allowed in program
	global singleton_removal  # boolean whether singletons are removed or not
	global full_only # boolean whether only 4/4 nucletide comparisons are allowed
	global lines_removed
	global average_per_location
	for argu in sys.argv[2:]:
		argulist = argu.split('=', 1)
		if argulist[0] == 'p': process_type = argulist[1]
		elif argulist[0] == 'd': comparison_distance = int(argulist[1])
		elif argulist[0] == 'n': ns_allowed = int(argulist[1])
		elif argulist[0] == 's': singleton_removal = True
		elif argulist[0] == 'f': full_only = True
		elif argulist[0] == 'lr': lines_removed = int(argulist[1])



def process_all(): # Runs all of the process
	trim_process()
	datafile = read_to_list(modified_file_input('_trimmed.csv'), 'csv') 
	linkage_disequilibrium_iterator(datafile)
	average_calculator_big_data('_calculated.csv') 





def control_process(): # The only function run in the main body. Controls which processes take place
	global process_type # which proccess takes place
	
	iterate_through_arguments() # Iterates through the arguments passed to the program, passing them off to the appropriate global variable
	if process_type == 'trim': trim_process()
	elif process_type == 'ldcalc': ld_process()
	elif process_type == 'avcalc': calculate_average_per_distance_process()
	elif process_type == 'all': process_all()
	elif process_type == 'regional': calculate_regions()
	elif process_type == 'regionald': calculate_regions(2,6, '_regionald.csv')
	elif process_type == 'regionalc': calculate_regions_set_region()
	elif process_type == 'avcalcc' : average_calculator_big_data_comparison()
	elif process_type == 'avecalc10kb': average_calculator_big_data_10kb()







def calculate_regions(ave_int = 1, calc_int = 7, output_string = '_regional.csv'):
	print sys.argv[1].split(":")[0]
	region_size = 1000 # size of the region bin
	data_holder = [['kb','average', 'n']] # will be appended with data. Currently holds just the header rows for the output file
	calc_reader = open(sys.argv[1].split(":")[0], 'rU') # reads the first half of the argument. The calculation file
	ave_reader = open(sys.argv[1].split(":")[1], 'rU') # reads the second half of the argument. The average file
	dict_holder = {} # holds the dictionary containing the data
	ave_dict = {} # holds the average values for each distance
	operational_avg_bp = 0.0 #holds the bp used for what region the r2 will go to
	for row in csv.reader(ave_reader, delimiter = ','): #populates the ave_dict with values from the given average file
		ave_dict[row[0]] = row[ave_int] 
	
	calc_reader.readline() #skips the header row in the calculation file
	for i,row in enumerate(csv.reader(calc_reader, delimiter = ',')): # moves through the calculation file
		if i % 100000 == 0: print "on row " + str(i)
		operational_avg_bp = float((float(row[0])+float(row[1]))/2.0) #sets the operational_avg_bp to the midpoint between the two bps in the file
		operational_avg_bp = float(int(region_size) * round(int(operational_avg_bp)/region_size)) #rounds down the the nearest region_size

		if operational_avg_bp in dict_holder: # Checks to see if the bp bin is already in the dictionary
			dict_holder[operational_avg_bp][0] += ( float(row[calc_int]) / float(ave_dict[str((int(row[1]) - int(row[0])))]) ) # if it is add the r2 value to the bins r2 holder
			dict_holder[operational_avg_bp][1] += 1 # and add one to the counter for the bin
		else:
			dict_holder[operational_avg_bp] = [( float(row[calc_int]) / float(ave_dict[str((int(row[1]) - int(row[0])))]) ), 1.0] #if not, make a key in the dictionary for the bin and add the r2 and 1 to the counter


	for key, values in dict_holder.items(): # go through the dictionary, adding to the list if thre are at least 100 values in the bin

		if int(values[1]) >= 100: data_holder.append([float(key)/float(region_size), float(values[0])/float(values[1]), values[1]])

	with open(sys.argv[1].split(":")[0][0:-4] + output_string, "wb") as f: # writes the data_holder list to a new file
		writer = csv.writer(f)
		writer.writerows(data_holder)




def calculate_regions_set_region(ave_int = 1, calc_int = 7, output_string = '_regional_comparion.csv'):
	print sys.argv[1].split(":")[0]
	data_holder = [] # will be appended with data. Currently holds just the header rows for the output file
	operational_avg_bp = 0.0 #holds the bp used for what region the r2 will go to
	ave_dict = {} # holds the average values for each distance


	calc_reader = open(sys.argv[1].split(":")[0], 'rU') # reads the first third of the argument. The calculation file
	ave_reader = open(sys.argv[1].split(":")[1], 'rU') # reads the second third of the argument. The average file
	region_size = open(sys.argv[1].split(":")[2], 'rU') # reads the last third of the argument. The regions from the Mcguagh supplementary data
	
	for row in csv.reader(ave_reader, delimiter = ','): #populates the ave_dict with values from the given average file
		ave_dict[row[0]] = row[ave_int]
		region_size

	for row in csv.reader(region_size, delimiter = ','):
		data_holder.append([row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], 0.0, 0.0]) # appends the region data into the data_holder I MIGHT NEED TO APPEND MORE!!!!
	del data_holder[0]
	calc_reader.readline()
	for i,row in enumerate(csv.reader(calc_reader, delimiter = ',')): # moves through the calculation file
		if i % 100000 == 0: print "on row " + str(i)
		operational_avg_bp = float((float(row[0])+float(row[1]))/2.0) #sets the operational_avg_bp to the midpoint between the two bps in the file
		for new_row in data_holder: # Here it attempts to find the right bin to place the r2 into.
			#print new_row
			if operational_avg_bp >= float(new_row[0]) and operational_avg_bp <= float(new_row[1]):
				new_row[8] += ( float(row[calc_int]) / float(ave_dict[str((int(row[1]) - int(row[0])))]) )
				new_row[9] += 1.0
				break

	for row in data_holder:
		row[8] = float(row[8])/float(row[9])
	data_holder.insert(0, ["Marker1", "Marker2", "#of_crossovers", "total", "bp_between_markers", "cM/MB",	"Low_cM/MB", "High cM/Mb", "average", "n"])

	with open(sys.argv[1].split(":")[0][0:-4] + output_string, "wb") as f: # writes the data_holder list to a new file
		writer = csv.writer(f)
		writer.writerows(data_holder)



def average_calculator_big_data_comparison(l=[]): # Replaces the average calculator. Reads file 1 line at a time instead of all at once.
	counterizer = 0 # Counts how many rows of data we have done
	average_list = [] 
	operational_avg_bp = 0.0 #holds the bp used for what region the r2 will go to
	temporary_file_list = []


	calc_reader = open(sys.argv[1].split(":")[0], 'rU') # reads the first half of the argument. The calculation file
	if not l:
		region_size = open(sys.argv[1].split(":")[1], 'rU') # reads the second half of the argument. The regions from the Mcguagh supplementary data
		region_size.readline() #skips header
		for row in csv.reader(region_size, delimiter = ','): #creates the list of Mgaugh data and an empty dictionary.
			average_list.append([row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], {}]) #empty dictionary will contain the bins 1-2000 of distance data
	else:
		
		average_list = l
	


	calc_reader.readline() # reads the header line

	for i,row in enumerate(csv.reader(calc_reader, delimiter = ',')):
		if i % 100000 == 0: print "on row " + str(i)
		operational_avg_bp = float((float(row[0])+float(row[1]))/2.0) #sets the operational_avg_bp to the midpoint between the two bps in the file
		for region in average_list:
			if operational_avg_bp >= float(region[0]) and operational_avg_bp <= float(region[1]):

				if row[2] in region[8]: # checks the distance metric in calculation file to see if it already exists in the dictionary if it does:
					region[8][row[2]][0] += 1.0 #add 1 to the counter
					region[8][row[2]][1] += float(row[6]) #add the value of Dprime to the running Dprime total
					region[8][row[2]][2] += float(row[7]) #add the value of r2 to the running r2 total
				else: 
					region[8][row[2]] = [1.0 ,float(row[6]), float(row[7])] #create a new dictionary key with the values of the 1, dprime, and r2
				break

	
	for row in average_list:
		temporary_file_list = []
		for key in row[8]:
			#print "key is " + str(key)
			#print str(row[8][key][2])
			temporary_file_list.append([key, (float(row[8][key][2])/float(row[8][key][0])), (float(row[8][key][1])/float(row[8][key][0])),  float(row[8][key][0]) ])
		temporary_file_list = sorted(temporary_file_list, key=lambda temporary_file_list: int(temporary_file_list[0])) #sorts the list, shouldn't be a problem. List will be 2000
		temporary_file_list.insert(0,['Distance','Average r2','Average Dprime', 'n'] ) # insert the header
		with open(sys.argv[1].split(":")[0][0:-15] + str(row[0]) + ".csv", "wb") as f: # writes the data_holder list to a new file
			#print "hello"
			writer = csv.writer(f)
			writer.writerows(temporary_file_list)


def average_calculator_big_data_10kb():
	lis = []
	for row in range(0,35000000, 500000):
		lis.append([row, (row + 500000), 0, 0 , 0, 0 , 0, 0,{}])

	average_calculator_big_data_comparison(lis)
		
















control_process()








