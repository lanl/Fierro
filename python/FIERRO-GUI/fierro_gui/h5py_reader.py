import h5py
import numpy as np
import os
import PIL.Image
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#f = h5py.File('./fierro_gui/ASTMD638_specimen.dream3d', 'r')


data_file_path = os.path.join(os.path.dirname(__file__), "Cylinder_Synthetic.dream3d")

with h5py.File(data_file_path, 'r') as f:
    print("Keys: %s" % f.keys())
    for key in f.keys():
        print(key)
    a_group_key = list(f.keys())[1]
    print((f[a_group_key])) 
    data = list(f[a_group_key])
    print(data)
    f = np.array(f)
    print ("Array: ", f[0][0])
    #PIL.Image.fromarray(f).show()

    # Change the Size of Graph using 
    # Figsize
    fig = plt.figure(figsize=(10, 10))
    
    # Generating a 3D sine wave
    ax = plt.axes(projection='3d')
    
    # Creating array points using 
    # numpy
    x = np.arange(0, 20, 0.1)
    y = np.sin(x)
    z = y*np.sin(x)
    c = x + y
    
    # To create a scatter graph
    ax.scatter(x, y, z, c=c)
    
    # turn off/on axis
    plt.axis('off')
    
    # show the graph
    plt.show()
    # preferred methods to get dataset values:
    #ds_obj = f[a_group_key]      # returns as a h5py dataset object
    #ds_arr = f[a_group_key][()]  # returns as a numpy array
    #PIL.Image.fromarray(ds_arr).show()