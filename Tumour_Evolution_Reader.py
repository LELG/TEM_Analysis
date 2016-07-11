
########################################################
###### Python Parser Reader of Tumour Evolution ######## 
########################################################
import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import sys

if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

print "File: ", sys.argv[1]
file = open(sys.argv[1], 'r')

data = ""
DataStructs = []

First_Time = True
i = 0
Test = True
stop = 100
TimeData = []

for line in file:
	if (line[0:1] == "T"):
		TimeData.append(line)
		if(First_Time):
			data = line
			First_Time = False
		else:
			Pandas_df = pd.read_csv(StringIO(data), sep = "\t", engine = 'python')
			Pandas_df.columns = ['ID', 'Size', 'MR','PR']
			DataStructs.append(Pandas_df)
			data = line
			i = i+1
	else:
		data += line

	if( i == stop and Test):
		break

Pandas_df = pd.read_csv(StringIO(data), sep = "\t", engine = 'python')
Pandas_df.columns = ['ID', 'Size', 'MR','PR']
DataStructs.append(Pandas_df)


x = []
y = []
newborn = []
clones = []
s_pressure = []
av_MR = []
av_PR = []

TE = dict()

j=0

for df, t in zip(DataStructs, TimeData):
	 y.append( float(sum(df['Size'])) )
	 av_MR.append( float(sum(df['MR'])/float(len(df['MR'])) ) )
	 av_PR.append( float(sum(df['PR'])/float(len(df['PR'])) ) )
	 headerData = t.split('\t', 4)
	 x.append( float(headerData[0].replace("T ", ""))/8760.0)
	 newborn.append( float(headerData[1].replace("N ", "") ) )
	 clones.append( float(headerData[2].replace("C ", "") ) )
	 s_pressure.append( float(headerData[3].replace("S ","").replace("\n", "") ) )


	 # This could be separte to avoid memory leaks
	 IDs = df['ID']
	 Sizes = df['Size']
	 MRs = df['MR']
	 PRs = df['PR']

	 for id_,size,mr,pr in zip(IDs, Sizes, MRs, PRs):
	 	if id_ in TE:
	 		TE[id_][0].append(size)
	 		TE[id_][1].append(mr)
	 		TE[id_][2].append(pr)
	 		TE[id_][3].append( float(headerData[0].replace("T ", ""))/8760.0 )
	 	else:
	 		TE[id_] = [[size],[mr],[pr], [float(headerData[0].replace("T ", ""))/8760.0]]

	# TE = df.set_index('ID').to_dict()
	# print TE

	 j= j +1
	 if(j == stop-5 and Test):
	 	break;


#fig = plt.figure()
#fig.suptitle('Clonal Size Plot', fontsize=14, fontweight='bold')
#ax = fig.add_subplot(111)
#ax.set_xlabel('years')
#ax.set_ylabel('Number Of Cells')

#for df in TE:
	#print TE[df][0]
#	ax.plot( np.asarray(TE[df][3]), np.asarray(TE[df][0]))



#fig = plt.figure()
#fig.suptitle('Clonal MR Plot', fontsize=14, fontweight='bold')
#ax = fig.add_subplot(111)
#ax.set_xlabel('years')
#ax.set_ylabel('Mutation Rate')

#for df in TE:
	#print TE[df][0]
#	ax.plot( np.asarray(TE[df][3]), np.asarray(TE[df][1]))



fig = plt.figure()
fig.suptitle('Clonal PR Plot', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Proliferation Rate')

for df in TE:
	#print TE[df][0]
	ax.plot( np.asarray(TE[df][3]), np.asarray(TE[df][2]))

#print x
#print y
#print newborn
#print clones
#print s_pressure

fig = plt.figure()
fig.suptitle('Estimated Number of Cells', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Number Of Cells')
ax.plot(np.asarray(x), np.asarray(y))

fig = plt.figure()
fig.suptitle('Estimated Number of Newborn Clones', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Number of Newborn Clones')
ax.plot(np.asarray(x), np.asarray(newborn))

fig = plt.figure()
fig.suptitle('Clones in Tumour', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Total Clones in Tumour')
ax.plot(np.asarray(x), np.asarray(clones))

fig = plt.figure()
fig.suptitle('Selective Pressure', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Selective Pressure')
ax.plot(np.asarray(x), np.asarray(s_pressure))

fig = plt.figure()
fig.suptitle('Estimated Clonal Average Mutation Rate', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Mutation Rate')
ax.plot(np.asarray(x), np.asarray(av_MR))

fig = plt.figure()
fig.suptitle('Estimated Clonal Average Proliferation Rate', fontsize=14, fontweight='bold')
ax = fig.add_subplot(111)
ax.set_xlabel('years')
ax.set_ylabel('Proliferation Rate')
ax.plot(np.asarray(x), np.asarray(av_PR))

plt.show()

#print TE['P-0:0']


######### Data For Open pandas dataframes
#for df in dataFrames:
#	file_ = open('temp_df.txt','w')
#	file_.write( df )
#	file_.close()
#	Pandas_df = pd.read_csv('temp_df.txt', sep="\t", engine='python')
#	Pandas_df.columns = ['ID', 'Size', 'MR','PR']
#	DataStructs.append(Pandas_df)



#print DataStructs[len(TimeData)-1]
#print DataStructs[0]
#print list(DataStructs[0])
