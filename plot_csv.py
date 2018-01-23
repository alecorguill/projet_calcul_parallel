#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv


def extract_csv(filename, save=False):
    # Extracting
    f = open(filename,'r')
    content = f.read().split("\n")
    sizes = [int(x.split(" ")[0]) for x in content[:-2]]
    rows = [float(x.split(" ")[-1]) for x in content[:-2]]
    return (sizes,rows)
    
if __name__ == "__main__":
    speedup = extract_csv("speedup.csv");
    plt.figure(1)
    plt.subplot(111)
    save = speedup[1][0]
    for i in range(len(speedup[0])):
        speedup[1][i] = save/float(speedup[1][i])
    # aos = plt.plot(aos[0], aos[1], 'b', label='aos-avx512')
    speedup = plt.plot(speedup[0], speedup[1], 'g', label='speedup')
    plt.legend(loc=1, bbox_to_anchor=(1,0.87))

    plt.xlabel('nombre de processus')
    plt.ylabel('Speed up')    
    plt.title('Speedup')    
    plt.show()
        
