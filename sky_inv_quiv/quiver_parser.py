import ast
import time
from typing import Any, Dict, List, Optional
import numpy as np
import math

from sky_inv_quiv.Field import Field
from sky_inv_quiv.Quiver import Quiver

def parse_quiver_file(path: str, field: Optional[Field] = None, grid_flag: bool = False) -> Quiver:
    """
    Parse a simple Quiver file
    """
    if field is None:
        field = Field("Q")

    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) < 2:
        raise ValueError("Less than two lines")

    # Find Quiver( call
    quiver_line = lines[0].strip()
    start_char = quiver_line.find("Quiver(")+len("Quiver(")
    comma = quiver_line.find(",", start_char)
    

    #parse vertices
    if grid_flag:
        grid_dims = quiver_line[start_char:comma].split(" ")
        n = int(grid_dims[0].strip())
        if len(grid_dims)==2:
            m = int(grid_dims[1].strip())
            vertices = [(i,j) for j in range(0, m ) for i in range(0, n )]
        elif (len(grid_dims)!=1) or (math.isqrt(n)**2!=n):
            print("not grid input format")
        else:
            vertices = [(i,j) for j in range(0,math.isqrt(n) ) for i in range(0, math.isqrt(n) )]
    else:
        n = int(quiver_line[start_char:comma])
        vertices = list(range(1, n + 1))
   
    #parse edges
    edges: Dict[Any, List[Any]] = {}
    str_edges = quiver_line[comma + 1 :].strip().strip(");")
    for edge_list in ast.literal_eval(str_edges):
       edges[str(edge_list[2])] = [vertices[int(edge_list[0]-1)], vertices[int(edge_list[1]-1)]]

    #line 2
    rep_line = lines[1].strip().strip(");")
    start_dimvec = rep_line.find("[") 
    startve = rep_line.find("],") + len("],")
    #parse dimvec
    dimvec_list = ast.literal_eval(rep_line[start_dimvec : startve - 1].strip())
    dimvec_list_dict = dict(zip(vertices, dimvec_list))

    #parse Ve
    Ve: Dict[Any, np.ndarray] = {}
    #initialize Ve with zero matrices
    for e in edges.keys():
        Ve[e] = field.to_Field(np.zeros((dimvec_list_dict[edges[e][1]], dimvec_list_dict[edges[e][0] ]), dtype="i"))
    #read nonzero Ve's
    Ve_list = ast.literal_eval(rep_line[startve :].strip())
    for ve in Ve_list:
        Ve[str(ve[0])] = field.to_Field(np.array(ve[1]).transpose())
    return Quiver(vertices, edges, Ve, field=field, grid=grid_flag)

