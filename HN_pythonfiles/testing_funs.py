import Field
import numpy as np
from time import time
import copy
import auxHN as aux
from mainHN import computeHN
import scipy.linalg

epsilon=1e-10


def random_change_bases(V):
    P={}
    for x in V.vertices:
        if V.spaces[x]==0:
            P[x]=V.field.to_Field(np.ones((0, 0)))
        else:
            p= np.random.randint(1,20,(V.spaces[x],V.spaces[x]))
            while (V.field.descr in ['Q','R','C'] and
                    not scipy.linalg.orth(p,rcond=epsilon).shape[1]==V.spaces[x]) or \
                    (V.field.descr not in ['Q','R','C'] and aux.extract_basis(V.field.to_Field(p),V.field).shape[1]!=V.spaces[x]):
                p= np.random.randint(1,20,(V.spaces[x],V.spaces[x]))
            P[x]=V.field.to_Field(p)
    return aux.bases_change(V,P)
    
def test_skyscraper(grid_size,n_int,verbose=False,field=Field.Field("Q") ):
    
    V,_=aux.int_module_in_grid_quiver(grid_size, [(1,1)], [(0,0)],field)
    list_V=[]
    if n_int>2:
        V=aux.grid_indec(grid_size,field)
        n_int-=2
        list_V=[copy.deepcopy(V)]
    
    for k in range(n_int):
        s=np.random.randint(np.zeros(len(grid_size)),1#(2*np.array(grid_size))//4
                            ,2)
        t0=np.random.randint(s,np.array(grid_size),2)
        t1=np.random.randint(s,np.array(grid_size),2)
        if verbose:
            print("source,taget: ",s,t0,t1)
        V1,supp1= aux.int_module_in_grid_quiver(grid_size, [tuple(s)], [tuple(t0),tuple(t1)],field)
        #V1,supp1= aux.int_module_in_grid_quiver(grid_size, [(min(s0),min(s1))], [(max(s0),max(s1))],field)
        V=aux.direct_sum(V, V1)
        list_V.append(V1)
    #V.display_graph("V")
    
    W=random_change_bases(V)

    #W.display_graph("W")
    x_set=False
    #x_set=[(0,0)]
    #np.random.seed(2)
    HN1= compute_HN_from_support(list_V,x_set=x_set)
    HN2=computeHN(W,x_set=x_set   ,verbose=verbose)
    #HN2bis= computeHN(W,x_set=x_set   ,verbose=False)
    


    
    assert(list(HN1.keys())==list(HN2.keys()))
    success=True
    for x0 in HN1.keys():
        try:
            assert(len(HN1[x0])==len(HN2[x0]) and all([Field.is_all_zero_mat(HN1[x0][i]-HN2[x0][i],V.field) for i in range(len(HN1[x0]))] ))
        except:
            print(x0)
            print("computing from supports:")
            print(HN1[x0])
            print("computing from HN:")
            print(HN2[x0])
            #print(len(HN1[x0])==len(HN2[x0]),[aux.is_all_zero_mat(HN1[x0][i]-HN2[x0][i],V.field) for i in range(len(HN1[x0]))])
            raise ValueError("i")
            success=False
    return success

    

    
def x_span_set(vertices,x):
    return dict(zip(vertices,[all(np.array(v)>=np.array(x)) for v in vertices] ))

def intersect_set(S1,S2):
    assert(list(S1.keys())==list(S2.keys()))
    return dict(zip(S1.keys(),[S1[v] and S2[v] for v in list(S1.keys())] ))

def sum_HN(HN1,HN2,vertices):
    null_dim_vect= np.zeros(len(vertices),dtype="i")
    cards=set(list(HN1.keys())+ list(HN2.keys()))
    HN= dict(zip(cards ,  [np.copy(null_dim_vect) for _ in cards])  )
    for card in list(HN1.keys()):
        HN[card]+=  HN1[card]
    for card in list(HN2.keys()):
        HN[card]+=  HN2[card]
    return HN




def compute_HN_from_support(V,x_set=False):
    vertices=V[0].vertices
    if not(x_set):
        x_set=vertices
    skyscraper={}
    for x in x_set:
        totdims=np.sum([ np.array( list(W.spaces.values()),dtype=np.int32) for W in V ],axis=0)
        HN={}
        #interval modules
        for W in V:
            if W.spaces[x]>0:
                HN= sum_HN( HN, HN_everywhere_ss(W,x) ,vertices)      
        quotient_dims = [a[1] for a in  sorted(HN.items())]
        filtr_dims=copy.copy(quotient_dims)
        for i in range(len(filtr_dims)-1):
            filtr_dims[i+1]+=filtr_dims[i] 
        skyscraper[x]=[np.zeros_like(totdims)]+filtr_dims
        if np.any( totdims-skyscraper[x][-1]>0):
            skyscraper[x]+=[totdims]
    return skyscraper

def HN_everywhere_ss(V,x):
        support=dict(zip(V.vertices,[V.spaces[v]!=0 for v in V.vertices]))
        assert( x in support)    
        S_cap_x=intersect_set(support, x_span_set(V.vertices, x))
        dim_intersect = np.array([  int(S_cap_x[y])*V.spaces[y ] for  y in V.vertices],dtype="i")
        card= np.sum(dim_intersect)/float(V.spaces[x])
        return {card:  dim_intersect} 



def test1(n=3,n_test=4, field=Field.Field("Q"), start_test = 42):
  xmax=[n,n]
  t= time()
  success_fail=[0,0]
  print(f"field={field.descr}, n={n}, n_test={n_test}, random seed={start_test}")
  for k in range(n_test):
      np.random.seed(start_test+k)
      success_fail[test_skyscraper(xmax,n,False,field)]+=1
  if success_fail[0]==0:
    print("success")
  else:
    print(success_fail[0] ,"failures:")  
  print("time per test=",np.round((time()-t)/n_test,2),"s")






def test2 (n=4, field = Field.Field("Q")):

  xmax=[n,n]
  V1,supp1= aux.int_module_in_grid_quiver( xmax, [(0,0)],[(0,2),(1,1)],field)
  V2,supp2= aux.int_module_in_grid_quiver( xmax, [(0,0)],[(2,1),(1,2)],field)
  V3=aux.grid_indec(xmax,field)

  V=aux.direct_sum(aux.direct_sum(V1, V2), V3)
  #V.display_graph("V before change of basis")
  V=random_change_bases(V)



  HN1= compute_HN_from_support( [ V1,V2,V3],False)
  HN2=computeHN(V)
  #print("computing from HN:")
  #print(HN2[(0,0)])
  #print("computing from supports:")
  #print(HN1[(0,0)])

  if len(HN1[(0,0)])==len(HN2[(0,0)]):
    print("success")
  else:
    print("failure")


#tests
print("test 1:")
test1(field = Field.Field("F_2"))
test1(n=5 , field = Field.Field("Q"))
print("----------------")
print("test 2:")
test2(field = Field.Field("F_2"))
test2(field = Field.Field("Q"))
