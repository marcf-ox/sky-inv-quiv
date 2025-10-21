#%%
import numpy as np
from sky_inv_quiv.auxHN import int_module_in_grid_quiver, direct_sum, ind_vertex
import copy
from sky_inv_quiv.Field import Field
import matplotlib.pyplot as plt
from scipy import stats
from time import sleep
from sky_inv_quiv.mainHN import computeHN, computeHN_sub
#aezaeaz

'''
V=aux.random_grid_indec([3,4], 5,3,3)

V.display_graph("V")


HNV=computeHN_sub(V,(0,0),verbose=False)

print(HNV)
plt.show()



'''




field=Field('F_2')



V1,supp1= int_module_in_grid_quiver( [3,3], [(0,0)],[(2,0),(0,2),(1,1)],field)
V2,supp2= int_module_in_grid_quiver( [3,3], [(0,0)],[(1,0),(0,1)],field)

x=(0,0)

V3= ind_vertex(V1.vertices,V1.edges,x,field)
 
W= direct_sum(direct_sum(V1,V2),V3)




W.display_graph("1_x",verbose=False)


HNW=computeHN_sub(W,(0,0),verbose=False)

print("test:",HNW)

plt.show()

# %%
