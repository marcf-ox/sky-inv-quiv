import numpy as np
import auxHN as aux
import copy
import matplotlib.pyplot as plt
from scipy import stats
#import mainHN as main
#aezaeaz


V=aux.random_grid_indec([5,5], 4)

#V.display_graph("V")

'''

V1,supp1= aux.int_module_in_grid_quiver( [2,2], [(0,0)],[(1,1)],field)
V2,supp2= aux.int_module_in_grid_quiver( [2,2], [(0,0)],[(1,0),(0,1)],field)

x=(0,0)

V3=aux.ind_vertex(V1.vertices,V1.edges,x,field)
 
W= aux.direct_sum(aux.direct_sum(V1,V2),V3)


W.display_graph("1_x")


HNW=main.computeHN_sub(W,(0,0),verbose=True)

print(HNW)
'''
