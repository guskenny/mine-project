#!/usr/bin/python
import sys, optparse, os, csv, math
order = []
data = {}
preloaded_std_devs = {}

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
		
		if row["Std Dev"]:
			preloaded_std_devs[row["Instance"]] = {}
			preloaded_std_devs[row["Instance"]][in_file] =row["Std Dev"]	


	read = True
	inFile.close()

for inst in order:
	avgs = {}
	square_diffs = {}
	std_devs = {}
	for file in sys.argv[3:-1]:
		avgs[file] = float(sum(data[inst][file])/len(data[inst][file]))
		square_diffs[file] = [math.pow(x - avgs[file],2) for x in data[inst][file]]
		if inst in preloaded_std_devs:
			if file in preloaded_std_devs[inst]:
				std_devs[file] = preloaded_std_devs[inst][file]
			else:
				std_devs[file] = '%.3E'%(math.sqrt(float(sum(square_diffs[file]))/len(square_diffs[file])))
	# std_dev = math.sqrt(float(sum(square_diff))/len(square_diff))
	key_max = max(avgs.keys(), key=(lambda k: avgs[k]))
	max_avg = avgs[key_max]

	avgs[sys.argv[-1]] = float(data[inst][sys.argv[-1]][0])

	print("\\texttt{"+inst.replace("_","\\_").lower()+"} & " +str('%.3E'%(data[inst]["LP_UB.csv"][0])) +" & " +str('%.3E'%(data[inst]["minelib.csv"][0])))
	for file in sys.argv[3:]:
		if file == sys.argv[-1]:
			# if abs(avgs[file] - max_avg) < 1e4:
			# 	print("& "+"\\best{"+str('%.3E'%(avgs[file]))+"} \\\\")
			# else:
			# 	print("& "+str('%.3E'%(avgs[file]))+" \\\\")
			print("& "+str('%.3E'%(avgs[file]))+" \\\\")
		else:
			if abs(avgs[file] - max_avg) < 1e4:
				print("& "+"\\best{"+str('%.3E'%(avgs[file]))+"} & "+str(std_devs[file]),end='')
			else:
				print("& "+str('%.3E'%(avgs[file]))+" & "+str(std_devs[file]),end='')
			print("")
	print("%")