import numpy as np
from matplotlib import pyplot as plt


t = np.linspace(0,1,500)

print(len(t))

x = np.sin(3.14 * 2 *5 *t)
print(x)

print(len(x))

ax = plt.axes()

# set limits
plt.xlim(0,1) 
plt.ylim(-1,1)

for i in range(500):        
     # add something to axes    
     #ax.scatter(t[i], x[i]) 
     ax.plot(t[i], x[i],'bd')

     # draw the plot
     plt.draw() 
     plt.pause(0.3) #is necessary for the plot to update for some reason

     # start removing points if you don't want all shown
     #if i>2:
     #    ax.lines[0].remove()
     #    ax.collections[0].remove()