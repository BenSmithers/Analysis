

def get_nth_fibb(n):
    previous = 1
    current  = 1
    
    ratio = float(current)/previous

    iterate = 0
    while iterate < n:
        
        current = current + previous
        previous = current - previous 
        ratio = float(current)/previous 

        iterate += 1


def get_golden_to_accuracy( accuracy ):

    last_one = 0

    previous = 1
    current  = 1
    
    ratio = float(current)/previous

    iterate = 0
    while abs(last_one - ratio)>accuracy:
        last_one = ratio

        current = current + previous
        previous = current - previous 
        print("{}/{}".format(current, previous))
        ratio = float(current)/previous 

        iterate+=1
        
    return(iterate)


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np


xs = np.logspace(1,20,20)
ys = [ get_golden_to_accuracy(x**-1) for x in xs]
plt.plot(np.log10(xs),ys)
plt.show()

