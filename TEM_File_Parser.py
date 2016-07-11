########################################################
########################################################
###### Python Parser Reader of Tumour Evolution ######## 
########################################################
########################################################
import pandas as pd
import sys
import os.path
import numpy as np
import matplotlib.pyplot as plt
import itertools
import re
from ete3 import Tree, TreeStyle, TextFace
from pandas.tools.plotting import parallel_coordinates
from pandas.tools.plotting import andrews_curves
from pandas.tools.plotting import radviz
from pandas.tools.plotting import autocorrelation_plot
from pandas.tools.plotting import lag_plot
from natsort import natsorted, ns

if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO




DataStructs = []
TimeSeries = []
Tumour_Evolution = dict()

Final_Population_df = pd.DataFrame()



x = []
y = []
newborn =[]
clones = []
s_pressure = []
av_MR = []
av_PR = []

Tumour_Burden = []
Years = []
Total_Clones =[]

class Vividict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class Parser(object):

	"""
		This Function Check for the validity of the 
		intpur file. Does not check if the file is in the correct format
	"""
	def Get_TEM_File(self):
		valid_file = False
		if ( os.path.isfile(sys.argv[1]) ):
			valid_file = True
		else:
			valid_file = False

		return valid_file

	"""
		This function extracts the data from the snapshot file
		It saves each spanshot as a Pandas df and each header as
		as list of strings separated by tabs.
	"""
	def Read_Tumour_Evolution_File(self, Test = True, stop = 100):
		file = open(sys.argv[1], 'r')
		data = ""
		First_Time = True
		i = 0
		
		for line in file:
			if (line[0:1] == "T"):
				TimeSeries.append(line)
				if(First_Time):
					data = line
					First_Time = False
				else:
					Pandas_df = pd.read_csv(StringIO(data), sep = "\t", engine = 'python')
					Pandas_df.columns = ['ID', 'Size', 'MR', 'PR']
					DataStructs.append(Pandas_df)
					data = line
					i = i+1
			else:
				data += line

			if (i == stop and Test):
				break

		Pandas_df = pd.read_csv(StringIO(data), sep = "\t", engine = 'python')
		Pandas_df.columns = ['ID', 'Size', 'MR', 'PR']
		DataStructs.append(Pandas_df)

	"""
		Pandas df final population size dataframe
	"""
	def Read_Final_Population_File(self):
		Final_Population_df = pd.read_csv(sys.argv[3], sep="\t", engine ='python').sort_values(by='Clone_Size', ascending=0)
		return Final_Population_df
		

	def Read_Final_Population_File_Top_values(self, Top_Percent_Values = 1):
		Final_Population_df = pd.read_csv(sys.argv[3], sep="\t", engine ='python').sort_values(by='Clone_Size', ascending=0)
		
		if( int(len(Final_Population_df) * Top_Percent_Values) > 5000 ):
			#Final_Population_df = Final_Population_df.nlargest( 5000, 'Clone_Size' ) 
			Final_Population_df = Final_Population_df.nlargest( int (len(Final_Population_df) * Top_Percent_Values), 'Clone_Size' ) 
		else:
			Final_Population_df = Final_Population_df.nlargest( int (len(Final_Population_df) * Top_Percent_Values), 'Clone_Size' ) 
		return Final_Population_df

	def Normalise_df_axis(self):
		IDs = Final_Population_df['Generation_ID']
		CSs = Final_Population_df['Clone_Size']
		max_CS = float(max(CSs))
		n_CSs = []
		MRs = Final_Population_df['Mutation_Rate']
		n_MRs = []
		max_MR = float( max(MRs)  )
		PRs = Final_Population_df['Proliferation_Rate']
		n_PRs = []
		max_PR = float( max(PRs) )

		Muts = Final_Population_df['Mutations']
		n_Muts = []
		max_Muts = float( max(Muts) )
		Bur = Final_Population_df['Burden']
		n_Bur = []
		max_Bur = float( max(Bur) )
		Alive = []
		
		#print "CS", CSs
		for _cs, _mr, _pr, _mut, _burd in zip(CSs, MRs, PRs, Muts, Bur):
			alive_f = True
			if(_cs == 0):
				alive_f = False
				Alive.append(1)
				n_CSs.append(float(_cs))
			else:
				clone_s = float(float(_cs)/max_CS)
				n_CSs.append(clone_s)

				#if(clone_s >= 0.1 ):
					#Alive.append(2)
					#alive_f = False
				
			
			if(_mr == 0):
				n_MRs.append(float(_mr))
			else:
				n_MRs.append( float( float(_mr)/max_MR ) )

			if(_pr == 0):
				alive_f = False
				Alive.append(1)
				n_PRs.append( float(_pr) )
			else:
				n_PRs.append( float(_pr)/max_PR )

			if(_burd == 0):
				n_Bur.append(0)
			else:
				n_Bur.append( float( float(_burd)/max_Bur ) )

			n_Muts.append( float(_mut)/max_Muts )

			if(alive_f):
				Alive.append(0)

			
		return pd.DataFrame( zip(n_CSs,n_MRs,n_PRs,n_Muts, n_Bur, Alive), columns=['Clone_Size', 'Mutation_Rate', 'Proliferation_Rate', 'Mutstions','Burden','Extinct' ] )
		


	
	"""
		Provided a snapshot df, this function generates 
		a clonal dictionary containing the evolution of
		PR, MR and Sizes
	"""
	def Generate_Clonal_Evolution_Dict(self, df, headerData):
		IDs   = df['ID']
		Sizes = df['Size']
		MRs   = df['MR']
		PRs   = df['PR']

		for id_, size, mr, pr in zip (IDs, Sizes, MRs, PRs):
			if id_ in Tumour_Evolution:
	 				Tumour_Evolution[id_][0].append(size)
	 				Tumour_Evolution[id_][1].append(mr)
	 				Tumour_Evolution[id_][2].append(pr)
	 				Tumour_Evolution[id_][3].append( float(headerData[0].replace("T ", ""))/8760.0 )
	 		else:
	 			Tumour_Evolution[id_] = [[size],[mr],[pr], [float(headerData[0].replace("T ", ""))/8760.0]]

	"""
		This generates the full snapshot average of clonal values.
		The average of all active clones given a minimum delta size 
		fraction of the total population.
	"""
	def Filtered_Clonal_Value(self, delta=0.1):
		i = 0
		for df, t in zip(DataStructs, TimeSeries):
			sizes = df['Size']
			population_size = float( sum(df['Size']) )
			
			df = df.drop( df[df.Size <= population_size * delta ].index )
			
			y.append( float( sum(df['Size'])) )
			av_MR.append( float( sum(df['MR']/float( len(df['MR'])) ) ) )
			av_PR.append( float( sum(df['PR']/float( len(df['PR'])) ) ) )
			
			headerData = t.split('\t', 4)
			
			x.append( float(headerData[0].replace("T ",""))/8760.0 )
			newborn.append( float(headerData[1].replace("N ","") ) )
			clones.append ( float(headerData[2].replace("C ","") ) )
			s_pressure.append( float(headerData[3].replace("S ","").replace("\n","") ) )
			DataStructs[i] = df
			self.Generate_Clonal_Evolution_Dict( df, headerData )
			i += 1

			
	"""
		The same as Filtered_Clonal_Values but without
		a filter. This can be computational intensive for some 
		values.
	"""
	def Raw_Clonal_Values(self):
		# Raw struct no filtering
		for df, t in zip(DataStructs, TimeSeries):
			y.append( float( sum(df['Size'])) )
			av_MR.append( float( sum(df['MR']/float( len(df['MR'])) ) ) )
			av_PR.append( float( sum(df['PR']/float( len(df['PR'])) ) ) )
			headerData = t.split('\t', 4)
			x.append( float(headerData[0].replace("T ",""))/8760.0 )
			newborn.append( float(headerData[1].replace("N ","") ) )
			clones.append ( float(headerData[2].replace("C ","") ) )
			s_pressure.append( float(headerData[3].replace("S ","").replace("\n","") ) )
			self.Generate_Clonal_Evolution_Dict( df, headerData )


	def Plot_Clonal_Evolution_Growth(self, plot=True, Log = False):
		fig = plt.figure()
		fig.suptitle('Clonal Evolution Growth', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Growth')

		if(Log):
			for df in Tumour_Evolution:
				ax.plot( np.asarray(Tumour_Evolution[df][3]), np.log10(np.asarray(Tumour_Evolution[df][0])) )

		else:
			for df in Tumour_Evolution:
				ax.plot( np.asarray(Tumour_Evolution[df][3]), np.asarray(Tumour_Evolution[df][0]))

		if(plot):
			plt.show()
		else:
			plt.savefig('Clonal_Growth.eps', format='eps', dpi=1000)


	def Plot_Clonal_Evolution_PR(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Clonal Evolution PR Plot', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('years')
		ax.set_ylabel('Proliferation Rate')

		for df in Tumour_Evolution:
			ax.plot( np.asarray(Tumour_Evolution[df][3]), np.asarray(Tumour_Evolution[df][2]))

		if(plot):
			plt.show()
		else:
			plt.savefig('Clonal_Evolution_PR.eps', format='eps', dpi=1000)


	def Plot_Clonal_Evolution_MR(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Clonal Evolution MR', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Mutation Rate')

		for df in Tumour_Evolution:
			ax.plot( np.asarray(Tumour_Evolution[df][3]), np.asarray(Tumour_Evolution[df][1]))

		if(plot):
			plt.show()
		else:
			plt.savefig('Clonal_Evolution_MR.eps', format='eps', dpi=1000)


	def Plot_Estimated_Growth(self, plot = True, Log = False ):
		fig = plt.figure()
		fig.suptitle('Estimated Number of Cells', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Number Of Cells')

		if(Log):
			ax.plot( np.asarray(x), np.log10( np.asarray(y) ) )
		else:
			ax.plot( np.asarray(x), np.asarray(y) )
		
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Growth.eps', format='eps', dpi=1000)




	def Plot_Estimated_Newborn_Clones(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Estimated Number of Newborn Clones', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Number Of Cells')
		ax.plot(np.asarray(x), np.asarray(newborn))
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Newborn_Clones.eps', format='eps', dpi=1000)

	def Plot_Estimated_Clones_in_Tumour(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Clones in Tumour', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Total Clones in Tumour')
		ax.plot(np.asarray(x), np.asarray(clones))
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Active_Clones_in_Tumor.eps', format='eps', dpi=1000)

	def Plot_Estimated_Selective_Pressure(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Clones in Tumour', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Selective Pressure')
		ax.plot(np.asarray(x), np.asarray(s_pressure))
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Selective_Pressure.eps', format='eps', dpi=1000)

	def Plot_Estimated_Clonal_Average_MR(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Estimated Clonal Average Mutation Rate', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Mutation Rate')
		ax.plot(np.asarray(x), np.asarray(av_MR))
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Clonal_Avg_MR.eps', format='eps', dpi=1000)


	def Plot_Estimated_Clonal_Average_PR(self, plot=True):
		fig = plt.figure()
		fig.suptitle('Estimated Clonal Average Proliferation Rate', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Proliferation Rate')
		ax.plot(np.asarray(x), np.asarray(av_PR))
		if(plot):
			plt.show()
		else:
			plt.savefig('Estimated_Clonal_Avg_PR.eps', format='eps', dpi=1000)


	def Valid_Final_Population_File(self):
		valid_file = False
		if(os.path.isfile(sys.argv[3])):
			valid_file = True
		else:
			valid_file = False

		return valid_file

	def Valid_Growth_File(self, Path =""):
		valid_file = False
		if ( os.path.isfile(Path)  ):
			valid_file = True
		else:
			valid_file = False

		return valid_file



	def Read_and_Parse_Growth_File(self, Path=""):
		file = open(Path, 'r')
		data = ""
		for line in file:
			headerData = line.split('\t', 3)
			Tumour_Burden.append( float(headerData[0]) )
			Years.append( float(headerData[1])/8760.0 )
			Total_Clones.append( float(headerData[2].replace("\n","")) )


	def Plot_Simulated_Tumour_Growth(self, plot = True):
		fig = plt.figure()
		fig.suptitle('Tumour Growth', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Tumour Growth')
		ax.plot(np.asarray(Years), np.asarray(Tumour_Burden))
		if(plot):
			plt.show()
		else:
			plt.savefig('Simulated_Tumour_Growth.eps', format='eps', dpi=1000)

	def Plot_Simulated_Number_of_Clones(self, plot = True):
		fig = plt.figure()
		fig.suptitle('Number of Clones', fontsize=14, fontweight='bold')
		ax = fig.add_subplot(111)
		ax.set_xlabel('Years')
		ax.set_ylabel('Number of Clones')
		ax.plot(np.asarray(Years), np.asarray(Total_Clones))
		if(plot):
			plt.show()
		else:
			plt.savefig('Simulated_Number_of_Clones.eps', format='eps', dpi=1000)

	def Surviving_Phylogenetic_Tree(self):
		C_IDs = Tumour_Evolution.keys()
		C_Times = []
		ID_lenghts  = []
		Phylogeny = []

		C_IDs.remove("P-0:0")
		main_branches = []
		
		#while loop

		regex_str = 'PC[0-9]+'
		branch = ',[0-9]+'
		entering_flag = True
		Phylo_Struct = []
		i = 0
		
		while(  C_IDs ):
			if(entering_flag):
				entering_flag = False
				searchRegex = re.compile(regex_str + '-.*')
				matches = [m.group(0) for l in C_IDs for m in [searchRegex.search(l)] if m]
				
				for m in matches:
					C_IDs.remove(m)
			

				Phylo_Struct.append(matches)
			else:
				regex_str = regex_str + branch
				searchRegex = re.compile(regex_str + '-.*')
				matches = [m.group(0) for l in C_IDs for m in [searchRegex.search(l)] if m]
	
				for m in matches:
					C_IDs.remove(m)

				Phylo_Struct.append(matches)
				
		#print "PS",Phylo_Struct
		# trabnsalte this into a tree
		main_branches = Phylo_Struct[0] 
		
		branch_ids = []
		for clone in main_branches:
			branch_ids.append(clone[0:3])

		Phylogeny = []

	
		initial_step = Phylo_Struct.pop(0)
		initial_step = natsorted(initial_step, key=lambda y: y.lower())
		#print "I", initial_step
		

		ID_time = dict()
		regex_str = 'PC[0-9]+'
		for clone in initial_step:
			parent_str = re.search(regex_str,clone)
			ID_time[parent_str.group(0)] = int(clone.split("-", 1)[1].split(":")[0])
			Phylogeny.append( ("P", clone, 100 - int(clone.split("-", 1)[1].split(":")[0]) ) ) ## Year length

		#print ID_time
		#print "Remaining ", Phylo_Struct

		## Generating Phylogenetic Tree
		regex_str = 'PC[0-9]+'
		branch = ',[0-9]+'
		for step in Phylo_Struct:
			step = natsorted(step, key=lambda y: y.lower())
			for clone in step:
				clone_str = re.search(regex_str,clone)
				for _parent in initial_step:
					parent_str = re.search(regex_str,_parent)
					if( clone_str.group(0) == parent_str.group(0)):
						main_parent = re.search('PC[0-9]+', clone)
						clone_year = int(clone.split("-", 1)[1].split(":")[0])
						Phylogeny.append( (_parent, clone, abs(ID_time[main_parent.group(0)] - clone_year)  ) )     ## year lengt normalised
			
			regex_str = regex_str + branch
			initial_step = step
			
		


		t = Tree.from_parent_child_table( Phylogeny  )#a=np.unique(t).tolist()
		ts = TreeStyle()
		ts.show_leaf_name = True
		

		#ts.rotation = 90

		
		ts.mode = "c"
		ts.arc_start = -180 # 0 degrees = 3 o'clock
		ts.arc_span = 180

		t.show(tree_style=ts)



	

		#print "B"
		#for branch in Phylo_Struct:
		#	for e in branch:
		#		for _id in branch_ids: # by regex instead
		#			if (e[0:3] == _id):
		#				if( _id in phy ):
		#					phy[_id].append(e)
		#				else:
		#					phy[_id] = [e]

		#print phy
		#keys =  phy.keys()

		#print "K"

		

		
		#print "ALL, IDs", C_IDs

		#for _id in C_IDs:
		#	if ( _id.find(",") == -1 ):
		#		main_branches.append(_id[0:3])
		#		C_IDs.remove(_id)

		#print "PT", C_IDs
		#print "MB", main_branches




	# Phylogeny tree plot
	def Get_Clonal_IDs(self):
		C_IDs = Tumour_Evolution.keys()
		C_Times = []
		ID_lenghts  = []
		Phylogeny = []

		C_IDs.remove("P-0:0")
		main_branches = []
		
		print "ALL, IDs", C_IDs

		for _id in C_IDs:
			if ( _id.find(",") == -1 ):
				main_branches.append(_id[0:3])
				C_IDs.remove(_id)

		print "Main, branches ", main_branches

		for branch in main_branches:
			Phylogeny.append( ("P", branch) )

		print "ph ", Phylogeny

		#print "Phylogeny", Phylogeny
		print "ids ", C_IDs

		for clone in C_IDs:
			
			generations = clone.split(",")
			head = generations[0]
			last_and_hour = generations[ len(generations)-1 ].split("-")
			last = last_and_hour[0]
			time = last_and_hour[1]
			generations.remove( generations[0] )
			generations.remove( generations[ len( generations )-1 ] )
			if (not generations):
				Phylogeny.append( (head, head+","+last+"-"+time) )
			else:
				print clone
				carry = head
				for element in generations:
					Phylogeny.append( (head, head+","+element) )
					carry += "," + element
				Phylogeny.append( (carry, carry+","+last+"-"+time) )


		print Phylogeny

		p = sorted(list(set(Phylogeny)))
		print "p", p
		t = Tree.from_parent_child_table( p  )#a=np.unique(t).tolist()
		ts = TreeStyle()
		ts.show_leaf_name = True
		#ts.rotation = 90
		t.show(tree_style=ts)


	def Clonal_Evolution_Multidimensional_Data(self):
		i = 0.0
		Clonal_Evolution_df = pd.DataFrame()
		for df in DataStructs:
			if( i == 0 ):
				t = [i] * len(df)
				Clonal_Evolution_df = df 
				Clonal_Evolution_df['t'] = pd.Series(t, index=Clonal_Evolution_df.index)
			else:
				t = [i] * len(df)
				df['t'] = pd.Series(t, index=df.index)
				Clonal_Evolution_df = pd.concat([Clonal_Evolution_df, df], ignore_index = True)
			
			i = i + 1.0

		C = Clonal_Evolution_df[ 'ID' ]
		S = Clonal_Evolution_df['Size']
		M = Clonal_Evolution_df[ 'MR' ]
		P = Clonal_Evolution_df[ 'PR' ]
		T = Clonal_Evolution_df[ 't'  ]

		Normalised_df = pd.DataFrame(  zip( T/max(T), S/max(S),  P/max(P), M/max(M), C ), columns = ['t','Size','PR','MR','ID']  )
		
		plt.figure()
		parallel_coordinates(Normalised_df, 'ID', colormap='jet' ).set_title("PC Plot")
		plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
		plt.savefig('Clonal_Evolution_Parallel_Coords_Plot.eps', format='eps', dpi=1000)

		plt.figure()
		andrews_curves(Normalised_df, 'ID', colormap='jet')
		plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
		plt.savefig('Clonal_Evolution_Andrews_Curves_Plot.eps', format='eps', dpi=1000)

		plt.figure()
		radviz(Normalised_df, 'ID', colormap='jet')
		plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
		plt.savefig('Clonal_Evolution_RadViz_Plot.eps', format='eps', dpi=1000)



if __name__ == '__main__':

	Growth_Flag         = True 
	Evolution_Flag      = False
	CE_Frequency_Filter = 0.01
	Clonal_Evol_Flag    = False
	Surviving_Phylogeny = True



	if( Parser().Valid_Growth_File(sys.argv[2]) and Growth_Flag ):
		print "Growth"
		Parser().Read_and_Parse_Growth_File(sys.argv[2])
		Parser().Plot_Simulated_Tumour_Growth(False)
		Parser().Plot_Simulated_Number_of_Clones(False)
 
	if( Parser().Get_TEM_File() and Evolution_Flag ) :
		print "Evolution"
		Parser().Read_Tumour_Evolution_File(False, 3000)
		#arser().Raw_Clonal_Values()
		Parser().Filtered_Clonal_Value(CE_Frequency_Filter)
		# To visualise
		Parser().Plot_Estimated_Growth(False, True)
		Parser().Plot_Estimated_Clones_in_Tumour(False)
		Parser().Plot_Estimated_Selective_Pressure(False)
		Parser().Plot_Estimated_Clonal_Average_MR(False)
		Parser().Plot_Estimated_Clonal_Average_PR(False)
		Parser().Plot_Clonal_Evolution_PR(False)	
		Parser().Plot_Clonal_Evolution_MR(False)	
		Parser().Plot_Clonal_Evolution_Growth(False, True)

	#For ploting clonal evolution  Parallel coordinates
	if( Parser().Get_TEM_File() and Clonal_Evol_Flag ):
		print "Clonal Evolution Multidimensional"
		Parser().Read_Tumour_Evolution_File(False, 3000)
		Parser().Raw_Clonal_Values()
		#Parser().Filtered_Clonal_Value(CE_Frequency_Filter)
		Parser().Clonal_Evolution_Multidimensional_Data()


	if( Parser().Get_TEM_File() and Surviving_Phylogeny ):
		print "Phylogeny "
		Parser().Read_Tumour_Evolution_File(False, 2000)
		#Parser().Filtered_Clonal_Value(CE_Frequency_Filter)
		Parser().Raw_Clonal_Values()
		Parser().Surviving_Phylogenetic_Tree()


	#Phylogenetic treee reconstructions
	if( Parser().Get_TEM_File() and False  ):
		Parser().Read_Tumour_Evolution_File(False, 2000)
		Parser().Filtered_Clonal_Value(0.1)
		#print Tumour_Evolution.keys()
		#print DataStructs[len(DataStructs)-1]
		#P_H =  Tumour_Evolution['P-0:0'][3]
		#print P_H[0], P_H[len(P_H)-1]
		print "PC Ploting"
		



	# Ploting final population parallel coordinates
	if( Parser().Valid_Final_Population_File() and False ):
		Final_Population_df = Parser().Read_Final_Population_File_Top_values(1)
		#Final_Population_df = Parser().Read_Final_Population_File()
		Final_Population_df = Parser().Normalise_df_axis()
		#print Final_Population_df
		print "Ploting", len(Final_Population_df)
		plt.figure()
		#parallel_coordinates(Final_Population_df, 'Extinct', color=['blue','black','red']).set_title("PC Plot")
		andrews_curves(Final_Population_df, 'Extinct', colormap='jet')
		plt.show()
		print "Done"
		


	

