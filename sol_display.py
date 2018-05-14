#!/usr/bin/python

import sys;
import os.path;
import Tkinter as tk;
import random as r;
import math;
import time;
import tkFont;

# program constants
WINDOW_WIDTH = 900;
LUM_SCALE = 0.9

def animation_loop(delay):
    # Check if there are still periods to be drawn
    if len(periods) > 0:
        t = periods[0]
        redrawPeriod(t)
        # Remove the period that was just drawn from the list of periods
        del periods[0]
        # Call the animation_loop again after a delay
        root.after(delay, animation_loop, delay)

def makePicture(t):
    img_fname = "./anim/"
    if t < 10:
        img_fname = img_fname + "image_0"+str(t)+".ps"
    else:
        img_fname = img_fname + "image_"+str(t)+".ps"
    canvas.postscript(file=img_fname, colormode='color')

def redrawPeriod(t):
    print "drawing period " + str(t)
    for b in xrange(len(blocks)):
        luminance = int(blocks[b][t] * 255);
        colour = "#%x%x%x" % (luminance,luminance,luminance);
        canvas.itemconfig(block_ids[b],fill=colour);
    canvas.itemconfig(period_disp, text=str(t))
    ones = str(ratios[t][0]) + "%"
    zeroes = str(ratios[t][1]) + "%"
    other = str(ratios[t][2]) + "%"
    canvas.itemconfig(ones_disp, text=ones)
    canvas.itemconfig(zeroes_disp, text=zeroes)
    canvas.itemconfig(other_disp, text=other)
    makePicture(t)


global args
args = sys.argv;

# test if arguments
if len(args) <= 1:
    print "No solution data found!\nExiting...\n";
    exit();

xsol_filename = str(args[1]) + ".xsol";

print "xsol_filename: " + xsol_filename;

# test if instance exists
if not os.path.isfile(xsol_filename):
    print "No solution data found: " + xsol_filename + " found!\nExiting...\n";
    exit();

def load_sol(fname):
    max_val = 0.0;
    with open(fname, 'r') as infile:
        line = "";
        blocks = [];
        line = next(infile);
        while not line.startswith('END'):
            values = line.split();
            values = [float(x) for x in values];
            blocks.append([0.0] + values);
            for v in xrange(len(values)):
                if (values[v] > max_val):
                    max_val = values[v];
            line = next(infile);

        for b in xrange(len(blocks)):
            for v in xrange(len(blocks[b])):
                if blocks[b][v] < 1.0:
                    blocks[b][v] = blocks[b][v] * (LUM_SCALE/max_val)
        return blocks;

# get info from file
blocks = load_sol(xsol_filename);

# blocks = [
#         [0.0, 0.2, 0.5, 0.8],
#         [0.3, 0.0, 0.4, 0.0],
#         [0.5, 0.7, 0.0, 1.0],
#         [0.0, 0.0, 0.3, 0.6],
#         [1.0, 0.5, 0.2, 0.9]
#       ];

nB = len(blocks);
nT = len(blocks[0]);

n_cols = int(math.sqrt(nB));
n_rows = n_cols;

if (nB % n_rows != 0):
   n_rows = n_rows + 1;

print "nB: " + str(nB);
print "nT: " + str(nT);
print "n_rows: " + str(n_rows);
print "n_cols: " + str(n_cols);

NODE_SIZE = WINDOW_WIDTH / n_cols;
WINDOW_HEIGHT = NODE_SIZE * n_rows + NODE_SIZE;

# initialise Tkinter
root = tk.Tk();

# initialise canvas
canvas = tk.Canvas(root, width=WINDOW_WIDTH+100, height=WINDOW_HEIGHT);

# draw background
canvas.create_rectangle(0,0,WINDOW_WIDTH+100,WINDOW_HEIGHT,outline="",fill="#002080")

# open canvas
canvas.pack();

# container for blocks
block_ids = [];

x = 0;
y = -1*NODE_SIZE;

# initialise nodes
for b in xrange(nB):
    if (b % n_cols == 0):
        x = 0;
        y = y + NODE_SIZE;

    luminance = int(blocks[b][0] * 255);
    colour = "#%x%x%x" % (luminance,luminance,luminance)
    id = canvas.create_rectangle(x, y, x+NODE_SIZE, y+NODE_SIZE, outline="",fill = colour);
    block_ids.append(id);
    x = x + NODE_SIZE;

ratios = [];

# compute ratios
for t in range(nT):
    ones = 0
    zeroes = 0
    other = 0
    for b in xrange(nB):
        if blocks[b][t] == 1.0:
            ones = ones + 1
        else:
            if blocks[b][t] == 0.0:
                zeroes = zeroes + 1
            else:
                other = other + 1
    ratios.append([int((ones/float(nB))*100),\
        int((zeroes/float(nB))*100),\
        int((other/float(nB))*100)])

period_font = tkFont.Font(size=50,weight='bold')
stat_font = tkFont.Font(size=30,weight='bold')
period_disp = canvas.create_text(WINDOW_WIDTH+30,50,text = '0',font=period_font,fill='red')

ones = str(ratios[0][0]) + "%"
zeroes = str(ratios[0][1]) + "%"
other = str(ratios[0][2]) + "%"

ones_disp = canvas.create_text(WINDOW_WIDTH+50,130,text = ones,font=stat_font,fill='white')
zeroes_disp = canvas.create_text(WINDOW_WIDTH+50,190,text = zeroes,font=stat_font,fill='black')
other_disp = canvas.create_text(WINDOW_WIDTH+50,250,text = other,font=stat_font,fill='gray')

periods = [0] + range(nT)

animation_loop(500)

# redrawPeriod(0)

#for t in xrange(nT):
#    root.after(5000, redrawPeriod,t);

root.mainloop();
