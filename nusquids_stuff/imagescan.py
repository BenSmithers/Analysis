import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import numpy as np

data = "/home/benito/software/data/depo_reco.png"

# the image is effectively grayscale, so we only need one channel
# it's like, ( red, green, blue, alpha), so let's grab the reds
img = mpimg.imread(data)[:,:,0]

shape = img.shape

total_data_bins_x = 98 #14 per generation X 7 generations
total_data_bins_y = 98

print("Dims: {} by {}".format(shape[1], shape[0]))

        
# data ranges from 0->1
def color_convert(color):
    """
    Converts what's read from the image into the desired datascale 
    """
    log_scale = True
    min_data = 10**-3
    max_data = 10**0

    min_color = 0.
    max_color = 1.
    
    value = 10**( 3*color - 3 ) if log_scale else 0.
    return(value)

def x_bin_to_bin(xdata):
    """
    Takes some bin index for the data itself (after reading)
    and returns its best guess for the correlated bin in the image
    """
    x_bins = shape[1]
    plot_bin = float(xdata)*x_bins/total_data_bins_x
    return(int(plot_bin))

def y_bin_to_bin(ydata):
    """
    Same as x_bin_to_bin, but for the y axis
    """
    y_bins = shape[0]
    plot_bin = float(ydata)*y_bins/total_data_bins_y
    return(int(plot_bin))


def y_convert(y_bin):
    """
    Converts the y-bin into the mean deposited energy
    """
    
    log_scale = True
    min_data = 10**0
    max_data = 10**7
    y_bins = shape[0]

    value = 10**( (7./y_bins)*y_bin )

    return(value)

def x_convert(x_bin):
    
    min_data = 10**0
    max_data = 10**7
    x_bins = shape[1]

    value = 10**((7./x_bins)*x_bin)

    return(value)

# these are the bin indices in terms of what we actually want 
xs = np.arange(total_data_bins_x)
ys = np.arange(total_data_bins_y)

data = np.zeros(shape=(total_data_bins_x, total_data_bins_y))

for x in xs:
    for y in ys:
        plot_x = x_bin_to_bin(x)
        plot_y = y_bin_to_bin(y)

        # the little inversion happens because, graphically, the origin is placed on the top-left
        # so 0 in the y axis is actually at the top! 
        # this is different than the way we do things though. We're civilized 
        data[x][y] = color_convert( img[shape[0]-plot_y-1][plot_x] )
        

#plt.contour(x_convert(xs), y_convert(ys), data)
x_vals = [x_convert(x_bin_to_bin(x)) for x in xs]
y_vals = [y_convert(y_bin_to_bin(y)) for y in ys]
plt.pcolormesh(x_vals, y_vals, np.log10(data), cmap='gist_gray')
plt.xscale('log')
plt.yscale('log')
#plt.imshow(data,interpolation='none', origin='lower', cmap='gist_gray')
plt.show()
