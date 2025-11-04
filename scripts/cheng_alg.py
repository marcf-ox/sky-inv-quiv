# script to apply HNcompute to the text argument given in cmd line
import sys   
from sky_inv_quiv import computeHN
import sky_inv_quiv.auxHN as aux
from sky_inv_quiv import Field
from sky_inv_quiv.quiver_parser import parse_quiver_file
from time import time  
from numpy import array

#input script python cheng_alg input_file.txt field
if __name__ == "__main__":
    input_file = sys.argv[1]
    field = sys.argv[2]

    f = open(f"{input_file}_output_{field}.txt","w")

    #parse the quiver file
    Q = parse_quiver_file(input_file,field = Field(field),grid_flag=True)

    t0 = time()
    HN = computeHN(Q)
    t1=time()
    for i_v in range(len(Q.vertices)):
        x = Q.vertices[i_v]  
        if Q.spaces[x]>0:
            print(f"At {x}:", file = f)
            for i_fil in range(1,len(HN[x])):
                dim_quotient = array(HN[x][i_fil]) - array(HN[x][i_fil-1])
                slope = dim_quotient[i_v]/sum(dim_quotient)
                print(f"({dim_quotient[i_v]},{round(slope, 3)})", file = f,end=" ")
            print("", file = f)
    print(f"Time: {t1 - t0}",file=f)
    f.close()

