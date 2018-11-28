#!/usr/bin/python3

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

plt.title('Void Fraction Evolution',fontsize=14,fontweight='bold')
plt.xlabel('Time Step',fontsize=12,fontweight='bold')
plt.ylabel('% Void Fraction',fontsize=12,fontweight='bold')
plt.plot(timestep, vf,'-b',linewidth=1.0,label='no LTS, CFL = 0.9')
plt.legend(loc='upper right')
plt.grid()
#plt.show()
plt.savefig('plot.png',bbox_inches='tight',dpi=400)
            
        
            

