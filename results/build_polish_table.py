#!/usr/bin/python
import sys, optparse, os, csv, math
order = []
data = {}

if len(sys.argv) < 3:
  print("Please enter csv files in order LP_UB, minelib, grasp random, merge (no split), merge, polished merge");
  exit();

read = False;

for in_file in sys.argv[1:]:

	inFile = open(str(in_file))
	reader = csv.DictReader(inFile)


	for row in reader:
		if not read:
			order.append(row["Instance"])

		if row["Instance"]:
			try:
				data[row["Instance"]][in_file].append(float(row["Best obj"]))
			except: 
				try:
					data[row["Instance"]][in_file] = [float(row["Best obj"])]
				except:
					data[row["Instance"]] = {}
					data[row["Instance"]][in_file] =[float(row["Best obj"])]
		
	read = True
	inFile.close()

for inst in order:

	print("\\texttt{"+inst.replace("_","\\_").lower()+"} & " +str('%.3E'%(data[inst]["LP_UB.csv"][0])) +" & " +str('%.3E'%(data[inst]["minelib.csv"][0])))
	for file in sys.argv[3:]:
		if file == sys.argv[-1]:
			# if abs(avgs[file] - max_avg) < 1e4:
			# 	print("& "+"\\best{"+str('%.3E'%(avgs[file]))+"} \\\\")
			# else:
			# 	print("& "+str('%.3E'%(avgs[file]))+" \\\\")data[inst][file]
			print("& "+str('%.3E'%(max(data[inst][file])))+"\\\\")
		else:
			print("& "+str('%.3E'%(max(data[inst][file]))))
	print("%")