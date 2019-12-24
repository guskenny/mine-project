#!/bin/bash

cmd="python draw_plot.py $1/SEED_0.csv $1/SEED_0_BEST.csv $1/TIME_0.csv $1/WALL_TIME_0.csv $1/RESTARTS_0.csv"

iter="all iterations"
iter_arg="0"

if [ "$#" == 3 ]
then
  iter="first $3 iterations"
  iter_arg="$3"
fi

if [ $# -gt 1 ]
then
  if [ $2 == "-g" ]
  then
    echo "Plotting $1 (groups) $iter"
    cmd="$cmd $1/GROUPS.csv -g $iter_arg"
  elif [ $2 == "-n" ]
  	then
  	echo "Plotting $1 (no data) $iter"
  	cmd="$cmd $1/GROUPS.csv -n $iter_arg"
  else
    echo "Plotting $1 (times) $iter"
    cmd="$cmd -t $iter_arg"
  fi
else
  echo "Plotting $1 (times)"
fi

$cmd

