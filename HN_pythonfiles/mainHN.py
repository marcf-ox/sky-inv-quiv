# -*- coding: utf-8 -*-
import numpy as np
import scipy as scipy
import auxHN as aux
import HNcshrunk as cshrunk
from time import time,sleep
import copy
from cfractions import Fraction
epsilon=1e-10
import warnings
import traceback
import matplotlib.pyplot as plt
import Field
warnings.simplefilter("error")



x_glob=(0,0)

def computeHN_sub(V,x,filtration=False,verbose=False):    # assumes V=<V_x> 
    d=1
    field=V.field

    maps_span= build_spanning_maps(V, x)
    maps_span.pop(x)
    maps_span_l= [ Aj for Aj in maps_span.values() if Aj.shape[0]*Aj.shape[1]!=0] 

        # trivially semistable
    if V.spaces[x] in [0,sum(list(V.spaces.values()))] :
        return[]
    
    v_x=V.spaces[x]
    v_notx=sum(list(V.spaces.values()))-v_x
    #print("v_x=",v_x,"v_notx=",v_notx,"l",list(V.spaces.values()))
    U0=None
    for v_x_small in range(int(min(np.ceil(v_x/float(v_notx)),v_x)),30+2*max(2,v_x+1)):
        U0=None
        blocks_list= cshrunk.buildblock(maps_span_l,d,v_x_small)
        try:  
            for r in blocks_list:
                blocks,blocks_size=r
                size_X= sum([s[0] for s in blocks_size[:,0]]), sum([s[1] for s in blocks_size[0]])
                M= np.block([[np.eye(min(blocks.shape)),np.zeros((min(blocks.shape), blocks.shape[1]-min(blocks.shape)))],
                             [np.random.randint(1,25*max(size_X) +10, size=(blocks.shape[0]-min(blocks.shape), min(blocks.shape))),  np.zeros( (blocks.shape[0]-min(blocks.shape),  blocks.shape[1]-min(blocks.shape))) ]])
                M=field.to_Field(np.array(M,dtype="i"))
                X= (M,maps_span_l, blocks) 
                A=cshrunk.assemble_block(X)
                if v_x_small>10:
                    print(blocks.shape,"t",aux.extract_basis(A,field).shape[1],aux.extract_basis(A,field).shape[1] *5 / blocks.shape[0]  )          
                try:
                    t=time()
                    U0_temp=cshrunk.wongblockpseudo(X,field=field)
                    t1=time()
                    '''
                    U0p=aux.extract_basis(cshrunk.wongblock(X)[:v_x],field)
                    print(U0_temp.shape,U0p.shape,"e")
                    t2=time()
                    if t1-t>100:
                        print("new", np.round(t1-t,2)," s, old",np.round(t2-t1,2)," s")
                    assert(aux.intersection(U0_temp,U0p,field).shape[1]==max(U0_temp.shape[1],U0p.shape[1]))
                    '''
                except ValueError as e: 
                    #recompute in case we are unlucky with M
                    if str(e) =="nc-rk A<nc-rk B":
                        M= np.block([[np.eye(min(blocks.shape)),np.zeros((min(blocks.shape), blocks.shape[1]-min(blocks.shape)))],
                                    [np.random.randint(1,25*max(size_X) +10, size=(blocks.shape[0]-min(blocks.shape), min(blocks.shape))),  np.zeros( (blocks.shape[0]-min(blocks.shape),  blocks.shape[1]-min(blocks.shape))) ]])
                        M=field.to_Field(np.array(M,dtype="i"))
                        X= (M,maps_span_l, blocks)
                        Ap=cshrunk.assemble_block(X)
                        U0_temp=cshrunk.wongblockpseudo(X,field=field)
                    else: 
                        raise 
                if not( U0 is None or aux.intersection(U0,U0_temp,field).shape[1]==max(U0.shape[1] ,U0_temp.shape[1]) ) :
                    
                    #c_plus
                    span_rep0_temp= aux.spanning_subrep(V, x, U0_temp)
                    subrep0_temp= aux.subrep(V,span_rep0_temp)
                    dims_subrep0_temp=np.array(list(subrep0_temp.spaces.values()))
                    c_plus= (v_notx+v_x)*subrep0_temp.spaces[x]- v_x*np.sum(dims_subrep0_temp)

                    #disc_plus
                    v_notx_s_plus= aux.ceil_div(v_notx*v_x_small,v_x)
                    #v_tot_s_min=  (v_notx*v_x_small)//v_x        +    v_x_small
                    v_x_disc_plus = v_x*(    v_notx_s_plus*subrep0_temp.spaces[x]- v_x_small*(np.sum(dims_subrep0_temp) -subrep0_temp.spaces[x]))
                    disc_plus = v_x_disc_plus/float(v_x_small) 

                     #c_moins
                    span_rep0= aux.spanning_subrep(V, x, U0)
                    subrep0= aux.subrep(V,span_rep0)
                    dims_subrep0=np.array(list(subrep0.spaces.values()))
                    c_moins= (v_notx+v_x)*subrep0.spaces[x]- v_x*np.sum(dims_subrep0)
                    v_notx_s_minus= (v_notx*v_x_small)//v_x
                    v_x_disc_minus = v_x*(    v_notx_s_minus*subrep0.spaces[x]- v_x_small*(np.sum(dims_subrep0) -subrep0.spaces[x]))
                    #compute list of possible ux, u_notx
                    better_U_possible = False
                    for ux in range(U0.shape[1]+1,U0_temp.shape[1]):
                        u_notx_max_plus = (v_notx*ux-c_plus)//v_x 
                        u_notx_max_minus = (v_notx*ux-c_moins)//v_x 
                        u_notx_max = min(u_notx_max_minus,u_notx_max_plus)
                        u_notx_min_plus = aux.ceil_div((v_x-v_x_disc_plus+v_x*ux*(v_notx_s_plus)),v_x_small*v_x)
                        assert( v_x_disc_plus//v_x> 
                                v_notx_s_plus*ux- v_x_small*u_notx_min_plus and 
                                 v_x_disc_plus//v_x <= 
                                    v_notx_s_plus*ux- v_x_small*(u_notx_min_plus-1))
                        u_notx_min_minus = aux.ceil_div((-v_x_disc_minus+v_x*ux*(v_notx_s_minus)),v_x*v_x_small)
                        u_notx_min =max( u_notx_min_minus,  u_notx_min_plus )

                        #print("ux",ux,",",u_notx_min,u_notx_max,",",v_x_small)
                        better_U_possible |=  u_notx_max >= u_notx_min
                    #raise error if bound disc< floor(disc+) not good enough
                    if better_U_possible:
                        raise ValueError("weight approx incorrect")
                    if c_plus<c_moins:
                        U_temp=U0
                U0=U0_temp
            break
        except ValueError as e:
            if str(e) =="weight approx incorrect":
                if v_x>2 and v_x_small>0 and verbose:
                    print("weigt err:",v_x_small,v_x)
            elif str(e)== "nc-rk A<nc-rk B":
                if field.descr!= "F_2" or v_x_small%5==4:
                    print("needs larger blow-up: (p,p0)=",v_x_small,v_x)
                if v_x_small>v_x:
                    pass
            else: 
                raise 

    #print(aux.intersection(Ublock,U,V.field).shape[1]==U.shape[1])


    U0=Field.flatten_zero(U0, field)
    if v_x>3 and v_x_small>1:#verbose:
        pass
        #print("U0=")
        #print(aux.print_frac_2darray(U0))
    if U0.shape[1]==V.spaces[x]:
        return[]
    span_rep= aux.spanning_subrep(V, x, U0)
    subrep= aux.subrep(V,span_rep)
    dims_subrep=np.array(list(subrep.spaces.values()))
    if v_x>6 and v_x_small>3:
        print("c=",(v_notx+v_x)*subrep.spaces[x]- v_x*np.sum(dims_subrep))
    l=[]
    if verbose:
        print("srep",np.sum(dims_subrep))
    if not(subrep.is_zero()):
        l=computeHN_sub(subrep,x,filtration,verbose)+[dims_subrep]
        quotrep= aux.quotientrep(V, span_rep)
        if verbose:
            print("qrep",sum(quotrep.spaces.values()))
        if not(quotrep.is_zero()):
            try:
                l=l+ [dims_subrep+ dims_quot for dims_quot in computeHN_sub(quotrep,x,filtration,verbose)]
            except ValueError as e:
                print("e2",e)
                #V.display_graph("V")
                #subrep.display_graph("U")
                #quotrep.display_graph("V/U")
                print('V', V.spaces.values(),'subrep', subrep.spaces.values(), 'quotrep', quotrep.spaces.values())
                raise
    return l
    
def computeHN(V,x_set=False,filtration=False,verbose=False): 
    if x_set==False:
        x_set=V.vertices
    skyscraper={}
    for x in x_set:
        if V.spaces[x]==0:
            skyscraper[x]= [np.zeros(len(V.vertices),dtype="i"),np.array(list(V.spaces.values()))]
        else:
            span_maps= aux.spanning_subrep(V, x,  V.field.to_Field(np.eye(V.spaces[x])))
            span_subrep= aux.subrep(V,span_maps)
            try:
                l=computeHN_sub(span_subrep,x,filtration,verbose)
            except ValueError as e:
                #V.display_graph("V2")
                #span_subrep.display_graph("U2")
                print("e",e)
                raise
            #V.display_graph("V",verbose)
            skyscraper[x]= [np.zeros(len(V.vertices),dtype="i")]+l+[np.array(list(span_subrep.spaces.values()))]
            if sum(list(span_subrep.spaces.values()))<sum(list(V.spaces.values())):
                #span_subrep.display_graph("<V_{"+str(x)+"}>",verbose)
                skyscraper[x]+= [np.array(list(V.spaces.values()))]
    return skyscraper



def build_spanning_maps(V,x):
    field=V.field
    maps_span={x: V.field.to_Field(np.eye(V.spaces[x]))}
    file=[x]
    while len(file)>0:
        for e in V.edges_out[file.pop()]:
            if  V.edges[e][1]  in maps_span.keys():
                #check commutativity
                try:
                    assert(Field.is_all_zero_mat(maps_span[V.edges[e][1]]-np.dot( V.Ve[e], maps_span[V.edges[e][0]]), field))
                except:
                    #V.display_graph("V")
                    print(x)
                    print("r",V.edges[e],maps_span[V.edges[e][1]],np.dot( V.Ve[e], maps_span[V.edges[e][0]]), field)
                    raise
            else: 
                file.append(V.edges[e][1])
                try:
                    maps_span[V.edges[e][1]]=np.dot( V.Ve[e], maps_span[V.edges[e][0]])
                except:
                    print("e",V.Ve[e], maps_span[V.edges[e][0]])
                    raise ValueError("error in map "+str(e)+" from "+str(V.edges[e][0])+" to "+str(V.edges[e][1]))
    return maps_span


