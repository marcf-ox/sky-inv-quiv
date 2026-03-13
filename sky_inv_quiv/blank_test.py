# %%
import matplotlib.pyplot as plt
plt.ion()

# %%
import numpy as np

import copy
from scipy import stats # type: ignore
from time import sleep




from sky_inv_quiv.mainHN import computeHN, computeHN_sub
from sky_inv_quiv.Field import Field
from sky_inv_quiv.naive_code import naive_sub_HN
from sky_inv_quiv.auxHN import int_module_in_grid_quiver, direct_sum, ind_vertex
from sky_inv_quiv.testing_funs import test_skyscraper, all_tests

# from mainHN import computeHN, computeHN_sub  # type: ignore
# from Field import Field  # type: ignore
# from naive_code import naive_sub_HN  # type: ignore
# from auxHN import int_module_in_grid_quiver, direct_sum, ind_vertex  # type: ignore
# from testing_funs import test_skyscraper, all_tests  # type: ignore

from auxHN import inverse_image_vect, null_space # type: ignore

field = Field("F_2")
A = field.to_Field(np.zeros((1,1),dtype="i"))
W = field.to_Field(np.ones((1,1),dtype="i"))
print(inverse_image_vect(A,W,field).shape)
print(null_space(A,field).shape)





#test advanced
all_tests()

#test naive
field=Field('Q')
V1,supp1= int_module_in_grid_quiver( [3,3], [(0,0)],[(2,0),(0,2),(1,1)],field)
V2,supp2= int_module_in_grid_quiver( [3,3], [(0,0)],[(1,0),(0,1)],field)
x=(0,0)
V3= ind_vertex(V1.vertices,V1.edges,x,field)
W= direct_sum(direct_sum(V1,V2),V3)
W= direct_sum(V2,V3)
W.display_graph("1_x", verbose=True)
HNW=computeHN_sub(W,(0,0))
HNW2=naive_sub_HN(W,(0,0))
assert(len(HNW)==len(HNW2))

my_naive_skyscraper = lambda V, x_set=None, filtration=False, verbose=False, comp_HN_sub=computeHN_sub: computeHN(
	V,
	x_set,
	filtration=filtration,
	verbose=verbose,
	comp_HN_sub=naive_sub_HN,
)
print("Test  for naive compute:", test_skyscraper((3,3), 4, True, field, my_naive_skyscraper ))


plt.show()

# %%
