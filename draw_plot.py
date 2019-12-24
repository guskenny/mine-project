#!/usr/bin/python
import matplotlib.pyplot as plt

import sys
import os

global args
args = sys.argv

# constants
LABEL = 1;
VALUE = 2;

# order of colours
COLOURS = ['b','r','g','c','m','y','k']

# double check there are files
if len(args) <= 1:
  print("Please enter csv files");
  exit();

iter_limit = int(args[len(args)-1])

save_dir = os.path.dirname(os.path.abspath(args[1]))

# container for y values
y_vals = []

# containers for metadata
lines=[1 for _ in range(len(args)-2)]
names=["" for _ in range(len(args)-2)]
pop_size=[0 for _ in range(len(args)-2)]
swaps=[0 for _ in range(len(args)-2)]
num_merges=[0 for _ in range(len(args)-2)]
restarts=[0 for _ in range(len(args)-2)]
cpu_times=[0 for _ in range(len(args)-2)]
wall_times=[0 for _ in range(len(args)-2)]
groups=[0 for _ in range(len(args)-2)]
max_sub_swaps=[0 for _ in range(len(args)-2)]
mip_timeout=[0 for _ in range(len(args)-2)]
mip_rel_gap=[0 for _ in range(len(args)-2)]

wall_times_idx = -1
times_idx = -1
restarts_idx = -1
skip_times = False

y_max = 0;

nodata = False

if args[len(args)-2] == "-n":
    nodata = True

for infile_name in args[1:len(args)-2]:
    with open(infile_name, 'r') as infile:
        count = 0
        y = []
        for line in infile:
            if iter_limit > 0 and count > iter_limit:
                break;
            count+=1
            if line.startswith("#"):
                continue
            if line.startswith("!"):
                split_line = line.split();
                if split_line[LABEL] == "LINE":
                    lines[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "CPU_TIME":
                    cpu_times[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "WALL_TIME":
                    wall_times[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "MIP_REL_GAP":
                    mip_rel_gap[args.index(infile_name)-1] = float(split_line[VALUE])
                if split_line[LABEL] == "NAME":
                    names[args.index(infile_name)-1] = split_line[VALUE]
                if split_line[LABEL] == "POP_SIZE":
                    pop_size[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "TIMES":
                    times_idx = args.index(infile_name)-1
                if split_line[LABEL] == "WALL_TIMES":
                    wall_times_idx = args.index(infile_name)-1
                if split_line[LABEL] == "RESTART":
                    restarts_idx = args.index(infile_name)-1
                if split_line[LABEL] == "SWAPS":
                    swaps[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "MAX_SUB_SWAPS":
                    max_sub_swaps[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "MIP_TIMEOUT":
                    mip_timeout[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "NUM_MERGES":
                    num_merges[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "RESTARTS":
                    restarts[args.index(infile_name)-1] = int(split_line[VALUE])
                if split_line[LABEL] == "GROUPS":
                    groups[args.index(infile_name)-1] = 1
                    skip_times = True
                continue

            try:
                y.append(float(line))
            except ValueError:
                continue;
            if args.index(infile_name) == 1:
            	y_max = max(y_max, float(line))

    y_vals.append(y)

box_text = "POP_SIZE: " + str(pop_size[0])\
    + "\nSWAPS_PER_ITER: " + str(swaps[0])\
    + "\nMAX_SUB_SWAPS: " + str(max_sub_swaps[0])
if restarts_idx > -1:
    box_text = box_text + "\nNUM_MERGES: " + str(int(len(y_vals[0])/pop_size[0] - y_vals[restarts_idx][-1]))\
    + "\nRESTARTS: " + str(int(y_vals[restarts_idx][-1]))
else:
    box_text = box_text + "\nNUM_MERGES: " + str(int(len(y_vals[0])/pop_size[0] - restarts[0]))\
    + "\nRESTARTS: " + str(restarts[0])

box_text = box_text +"\nMIP_TIMEOUT: " + str(mip_timeout[0]) + " sec"\
                    +"\nMIP_REL_GAP: " + str('%.2f'%(mip_rel_gap[0]*100)) + "%"

if wall_times_idx > -1:
    box_text = box_text + "\nWALL_TIME: " + str(y_vals[wall_times_idx][-1]) + " sec"
else:
    box_text = box_text + "\nWALL_TIME: " + str(wall_times[0]) + " sec"

if times_idx > -1:
    box_text = box_text + "\nCPU_TIME: " + str(y_vals[times_idx][-1]) + " sec"
else:
    box_text = box_text + "\nCPU_TIME: " + str(cpu_times[0]) + " sec"

fig,ax1 = plt.subplots(figsize=(12,8))

ax1.set_ylabel('Objective value')

ax2 = plt.twinx()

ax2.yaxis.set_label_position("right")
ax = ax1.twinx()

x_max = 0

ax1.set_xlabel('Iterations')
plt.title("Plot of " + names[0] + " - MAX: " + "{:.3E}".format(y_max))

for y_idx in range(len(y_vals)):
    if y_idx == wall_times_idx or y_idx == restarts_idx or (y_idx == times_idx and skip_times is True):
        continue
    line_type = ''
    line_style = 'solid'
    marker_size = 0
    if  lines[y_idx] == 1:
        line_type = 'o'
        marker_size = 4
        line_style = 'None'

    if groups[y_idx] == 1:
        ax2.plot([x * pop_size[y_idx] + pop_size[y_idx] for x in range(len(y_vals[y_idx]))],y_vals[y_idx], linewidth=lines[y_idx],color='gold', marker=line_type, mew=0, ls=line_style, ms=marker_size)

        ax2.set_ylabel('Number of groups per merge')
    elif y_idx == times_idx and skip_times is False:
         ax2.plot(range(len(y_vals[y_idx])),y_vals[y_idx], linewidth=lines[y_idx],color='gold', marker=line_type, mew=0, ls=line_style, ms=marker_size)

         ax2.set_ylabel('CPU time [sec]')
    else:
        ax1.plot(range(len(y_vals[y_idx])),y_vals[y_idx], linewidth=lines[y_idx],color=COLOURS[y_idx], marker=line_type, mew=0, ls=line_style, ms=marker_size)
        x_max = max(x_max,len(y_vals[y_idx]))

ax1.minorticks_on()
ax1.grid(True,which='major')
ax1.grid(True,which='minor', linestyle='--',alpha=0.5)

ax.yaxis.set_visible(False)
if nodata == False:
    ax.text(0.96,0.06, box_text,transform=ax.transAxes, fontsize=14, bbox={'facecolor':'white', 'alpha':0.5, 'pad':10}, ha='right', va='bottom')

ax1.set_xlim(right=x_max)

plt.tight_layout()
plt.rcParams["savefig.directory"] = save_dir
plt.show()

filename = str(save_dir) + "/" + names[0] + '_' + str(pop_size[0]) + '_' + str(num_merges[0]) 

if args[len(args)-2] == "-g":
    filename+= '_groups'
elif args[len(args)-2] == "-t":
    filename+= '_times'

# if (iter_limit > 0):
#     filename+="("+str(iter_limit)+")"
    
filename+= '.png'

fig.savefig(filename)
