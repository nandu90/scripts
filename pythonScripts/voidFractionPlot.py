#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

flag = 0
temparr = []
timestep = []
vf = []
tempstr = ""
with open("1703727.output","r") as f:
    for line in f:
        temparr = []
        flag = 0
        for word in line.split():
            temparr.append(word)
            if word == "Fraction:":
                flag = 1
        if flag == 1:
            timestep.append(int(temparr[0]))
            tempstr = temparr[3].replace("%;","")
            vf.append(float(tempstr))
            
        
            

