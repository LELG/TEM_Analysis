##############################################################
##############################################################
######			 Folder Handling Data Anlysis			######
##############################################################
##############################################################

##
## This file require a pathq as input to strt the analysis
# $ python Data_Analysis.py "/Volumes/bioinf head/Data/V2/Thu_Jul__7_09_36_55_2016"



import sys
import os.path
from TEM_File_Parser import *

if __name__ == '__main__':

	#Use argv to get input from command line
	# argv[0] will contain name of python code
	# argv[1] will contain path to riun the analysis
	Data_Repository_Path = sys.argv[1]
	Simulations = []
	Data_Paths = []
	print Data_Repository_Path
	if ( os.path.exists( Data_Repository_Path ) ):
		Simulations = os.listdir( Data_Repository_Path ) # returns list
	else:
		print "Invalid path"

	print "Simulations: ", Simulations

	#print Data_Repository_Path  + "/" + Simulations[0]

	for simulation in Simulations: # Ooen IT Folder iteration
		inside_path =  Data_Repository_Path + "/" + simulation
		for iteration in os.listdir( inside_path ): # Open ID Iteration
			ID_path = inside_path + "/" + iteration
			for _id in os.listdir(ID_path): # Check for Files
				Data_Paths.append(ID_path + "/" + _id)
				if (not os.path.isfile( ID_path + "/" + _id + "/Evolution/V2_Final_Population.txt") and os.path.isfile( ID_path + "/" + _id + "/Evolution/Tumour_Evolution.txt") and os.path.isfile( ID_path + "/" + _id + "/Growth/Tumour_Growth.txt" ) ):
					print "Filde not found in: ", ID_path + "/" + _id
	

	# Send the info to the Parser Object

	par = Parser()
	for path in Data_Paths:
		if ("Growth" in path):
			growth_file = path + "/Tumour_Growth.txt"
			if( par.Valid_Growth_File(growth_file) ):
				print "Processing... ", growth_file
				
				Parser().Read_and_Parse_Growth_File(growth_file)
				Parser().Plot_Simulated_Tumour_Growth(False)
				Parser().Plot_Simulated_Number_of_Clones(False)

	#if( Parser().Valid_Growth_File() ):
	#	print "Growth"
	#	Parser().Read_and_Parse_Growth_File()
	#	Parser().Plot_Simulated_Tumour_Growth(False)
	#	Parser().Plot_Simulated_Number_of_Clones(False)