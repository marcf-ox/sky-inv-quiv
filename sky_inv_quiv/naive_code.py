
from typing import Any, List, Dict, Optional
import numpy as np

from sky_inv_quiv.Quiver import Quiver
from sky_inv_quiv.mainHN import build_spanning_maps
from sky_inv_quiv.Field import Field, build_block_diag_l
from sky_inv_quiv.auxHN import extract_basis, inverse_image_vect,spanning_subrep, subrep, quotientrep, print_frac_2darray, null_space, intersection


def naive_sub_HN(V:Quiver,x,filtration:bool=False,verbose:bool=False) -> List[np.ndarray]:
    """
    Naive code for HN filtration without any empirical optimization
    """
    field = V.field
    #trivially semistable
    if V.spaces[x] in [0,sum(list(V.spaces.values()))]:
        return[]
    #build the matrix space
    maps_span= build_spanning_maps(V, x)
    maps_span.pop(x)
    maps_span_l= [ Aj for Aj in maps_span.values() if Aj.shape[0]*Aj.shape[1]!=0] 
    n=len(maps_span_l)
    v_x=V.spaces[x]
    v_notx=sum(list(V.spaces.values()))-v_x
    # build random A
    U: Optional[np.ndarray] = None
    for d in range(1,v_x*v_notx+1):
        max_val_random = 2*d*v_x*v_notx+10
        A = np.block(
            [
                [
                field.to_Field(np.random.randint(0 , max_val_random ,size = (1,1)) ) * maps_span_l[i% n]
            for j in range (d*v_notx) 
            ]
            for i in range (d*v_x*n)
            ]
        )
        U = naive_shrunk(maps_span_l,A, field)
        if U is not None:
            break
    if U is None:
        raise ValueError("shrunk not found")
    #print("U:",U)
    span_rep= spanning_subrep(V, x, U)
    subrep_obj= subrep(V,span_rep)
    dims_subrep=np.array(list(subrep_obj.spaces.values()))
    l=[]
    if not(subrep_obj.is_zero()):
        l=naive_sub_HN(subrep_obj,x)+[dims_subrep]
        quotrep= quotientrep(V, span_rep)
        if not(quotrep.is_zero()):
            l=l+ [dims_subrep+ dims_quot for dims_quot in naive_sub_HN(quotrep,x)]
    return l



def naive_shrunk(maps_span_l:List[np.ndarray],A:np.ndarray, field:Field)->Optional[np.ndarray]:
    """
    Naive code to find a shrunk subspace without any empirical optimization
    """
    W = field.to_Field(np.zeros((A.shape[0],0),dtype="i"))
    vx = maps_span_l[0].shape[1]
    n=len(maps_span_l)
    dim_W=0
    for _ in range(min(A.shape)):
        U = np.concatenate([inverse_image_vect(A,W,field), null_space(A,field)], axis = 1)
        W_blocks = []
        for i in range (n):
            block = extract_basis(np.concatenate([ np.dot(maps_span_l[i] , U[vx*j:vx*(j+1),: ]) for j in range(A.shape[0]//vx) ], axis = 1), field)
            W_blocks.append(block)
        W_first_block = build_block_diag_l(W_blocks, field)
        W = np.kron(field.to_Field(np.eye(A.shape[0]//W_first_block.shape[0], dtype=int)), W_first_block)
        if (W.shape[1] == dim_W):
            WcapImA= intersection(W,extract_basis(A,field),field)
            if (W.shape[1] == WcapImA.shape[1]):
                return extract_basis(U[:vx,:], field)
        dim_W = W.shape[1]
    print(A.shape, extract_basis(A,field).shape,W.shape, extract_basis(W,field).shape)
    print(intersection(W,A,field).shape, extract_basis(intersection(W,A,field),field).shape)
    return None
    




    


    




