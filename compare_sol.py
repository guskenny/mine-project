#!/usr/bin/python

BLOCK = 0;
DEST = 1;
TIME = 2;
Y_VAL= 3;

import sys;

global args
args = sys.argv;

if len(args) <= 2:
  print "Please enter two filenames to compare";
  exit();

solution = {};

match = True;

with open(args[1], 'r') as in1:
  for line in in1:
    if line.startswith("#"):
      continue;

    split_line = line.split();
  
    block = split_line[BLOCK];
    dest = split_line[DEST];
    time = split_line[TIME];
    y_val = split_line[Y_VAL];
  
    solution[float(block)] = (float(dest), float(time), float(y_val));

with open(args[2], 'r') as in2:
  for line in in2:
    if line.startswith("#"):
      continue;

    split_line = line.split();

    if float(split_line[BLOCK]) in solution:
      for i in range(len(split_line)-1):
        if float(split_line[i+1]) != solution[float(split_line[BLOCK])][i]:
          print "Error found in block:", split_line[BLOCK];
          err_field = ""
          if i+1 == DEST:
            err_field = "Destination"
          elif i+1 == TIME:
            err_field = "Time"
          elif i+1 == Y_VAL:
            err_field = "y"

          print err_field, "values", float(split_line[i+1]), "and", solution[float(split_line[BLOCK])][i], "dont match!"
          match = False;
          #exit();
if match:
  print "Solution files match perfectly!"
else:
  print "Solution files do not match!"
