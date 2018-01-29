#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
import copy

def extract_csv(filename, save=False):
    # Extracting
    f = open(filename,'r')
    content = f.read().split("\n")
    sizes = [int(x.split(" ")[0]) for x in content[:-1]]
    rows = [float(x.split(" ")[-1]) for x in content[:-1]]
    return (sizes,rows)
    
if __name__ == "__main__":
    plt.figure(1)
    plt.subplot(111)
    speedup = extract_csv("speedup.csv");
    save = speedup[1][0]
    eff = [x[:] for x in speedup]
    xy = [[],[]]
    for i in range(len(speedup[0])):
        speedup[1][i] = save/float(speedup[1][i])
        eff[1][i] = speedup[1][i]/float(eff[0][i])
    for i in range(1,len(speedup[0])+1):
        xy[0].append(i)
        xy[1].append(i)
    plt.plot(xy[0],xy[1],'r')
    plt.plot(eff[0], eff[1], 'b', label='efficiency')
    plt.plot(speedup[0], speedup[1], 'g', label='speedup')
    plt.legend(loc=1, bbox_to_anchor=(1,0.87))
    plt.xlabel('nombre de processus')
    plt.ylabel('Speed up')    
    plt.title('Speedup')    
    plt.show()
        
