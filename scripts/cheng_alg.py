# script to apply HNcompute to the text argument given in cmd line
import argparse
from sky_inv_quiv import computeHN
from sky_inv_quiv import Field
from sky_inv_quiv.naive_code import naive_sub_HN
from sky_inv_quiv.quiver_parser import parse_quiver_file
from time import time  
from numpy import array

#input script python cheng_alg input_file.txt field
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute HN output for a quiver file.")
    parser.add_argument("input_file", help="Path to input quiver text file")
    parser.add_argument("field", help="Field descriptor, e.g. Q or F_2")
    parser.add_argument(
        "--naive",
        action="store_true",
        help="Use naive sub-HN algorithm (naive_sub_HN) instead of default computeHN_sub",
    )
    args = parser.parse_args()

    input_file = args.input_file
    field = args.field

    suff_naive = "_naive" if args.naive else ""
    f = open(f"{input_file}_output_{field}{suff_naive}.txt","w")

    #parse the quiver file
    Q = parse_quiver_file(input_file,field = Field(field),grid_flag=True)

    t0 = time()
    if args.naive:
        HN = computeHN(Q, comp_HN_sub=naive_sub_HN)
    else:
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

