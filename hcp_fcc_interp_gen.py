#!/usr/bin/env python
# Purpose: Generate HCP and FCC lattice suitable for interpolation 
# Author: Raymond Gasper
# Adapted from FCC surface generator by Sungho Kim and Laalitha Liyanage

import os
import sys
from math import sqrt

usage="""
        Usage: ./hcp_fcc_interp_gen.py a_hcp c/a a_fcc na nb nc
        nc must multiply be a multiple of 3

        Mandatory arguments
        -------------------
        a_hcp - equilibrium lattice constant (hcp)
        c/a - equilibrium c/a ratio (hcp)
        a_fcc - equilibrium lattice constant (fcc)
    
        Optional Arguments:
        -------------------
        nx, ny, nz: number of repeated unit cells in the a,b,c directions
        default-  3, 3, 3
"""

def gen_data_hcp(a, c_a, na=3, nb=3, nc=3):
    C=[]; B=[]; A=[]

    x0,y0,z0 = 0.0, 0.0, 0.0
    b1= [0.5*na*a, -1*sqrt(3)/2*nb*a, 0]
    b2= [0.5*na*a, sqrt(3)/2*nb*a, 0]
    b3= [0, 0, c_a*nc*a]

    x = 1./na
    x_shift = x*1./3.
    y = 1./nb
    y_shift = y*2./3.
    z = 1./(2.*nc)    
    for k in range(0,2*nc,2):    
        for j in range(nb):
            for i in range(na):
                A.append(x0+i*x); B.append(y0+j*y); C.append(z0+k*z)
                A.append(x0+i*x+x_shift); B.append(y0+j*y+y_shift); C.append(z0+(k+1)*z)
                
    return A,B,C,b1,b2,b3


def gen_data_fcc(a0, na=3, nb=3, nc=2):
    C=[]; B=[]; A=[]
    
    x0,y0,z0 = 0.0, 0.0, 0.0
    a= a0*sqrt(2)/2
    
    b1= [0.5*na*a, -1*sqrt(3)/2*nb*a, 0]
    b2= [0.5*na*a, sqrt(3)/2*nb*a, 0]
    b3= [0, 0, sqrt(3)*a0*nc]
    x = 1./na    
    x_shift = x*1./3.
    y = 1./nb
    y_shift = y*2./3.
    z = 1./(3.*nc)
    for k in range(0,3*nc,3):
        for j in range(nb):
            for i in range(na):
                A.append(x0+i*x); B.append(y0+j*y); C.append(z0+k*z)
                A.append(x0+i*x+x_shift); B.append(y0+j*y+y_shift); C.append(z0+(k+1)*z)
                A.append(x0+i*x+2*x_shift); B.append(y0+j*y+2*y_shift); C.append(z0+(k+2)*z)
    return A,B,C,b1,b2,b3

def gen_poscar(xa,ya,za,bx,by,bz,filename):
  fout = open(filename,"w")
  fout.write("Title Goes Here\n")
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(bx[0],bx[1],bx[2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(by[0],by[1],by[2]))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(bz[0],bz[1],bz[2]))
  fout.write("%d\n"%len(xa))
#  fout.write("Selective Dynamics\n")
  fout.write("Direct\n")
  for i in range(len(xa)):
    fout.write("%22.16f %22.16f %22.16f\n"%(xa[i],ya[i],za[i]))
#    fout.write("%22.16f %22.16f %22.16f F F T\n"%(xa[i],ya[i],za[i]))
  fout.close()
  return None


#-------------------------------Main program---------------------------------------------------
#       example input
#./hcp_fcc_interp_gen.py a_hcp c_a a_fcc (opt)na (opt)nb (opt)nc
#   if python in front of other arguments make a copy of sys.argv and delete python entry to generalize the code to silly users
inputs=sys.argv[:]
if inputs[0] == 'python':
    del(inputs[0])

if len(inputs) == 4:
    a_hcp = float(inputs[1])
    c_a_ratio = float(inputs[2])
    a_fcc = float(inputs[3])
    na = 3
    nb = 3
    nc = 3
    Ah, Bh, Ch, b1h, b2h, b3h = gen_data_hcp(a_hcp,c_a_ratio)
    gen_poscar(Ah, Bh, Ch, b1h, b2h, b3h, "HCP%d%d%d.vasp"%(na,nb,nc))
    Af, Bf, Cf, b1f, b2f, b3f = gen_data_fcc(a_fcc)
    gen_poscar(Af, Bf, Cf, b1f, b2f, b3f, "FCC%d%d%d.vasp"%(na,nb,int(nc*2/3)))

elif len(inputs) == 7:
    a_latt = float(inputs[1])
    c_a_ratio = float(inputs[2])
    out_style = float(inputs[3])
    na = int(inputs[4])
    na = int(inputs[5])
    na = int(inputs[6])
    Ah, Bh, Ch, b1h, b2h, b3h = gen_data_hcp(a_latt,c_a_ratio,na,nb,nc)
    gen_poscar(Ah, Bh, Ch, b1h, b2h, b3h, "HCP%d%d%d.vasp"%(na,nb,nc))
    Af, Bf, Cf, b1f, b2f, b3f = gen_data_fcc(a_fcc,na,nb,int(nc*3/2))
    gen_poscar(Af, Bf, Cf, b1f, b2f, b3f, "FCC%d%d%d.vasp"%(na,nb,int(nc*2/3)))

else:
    print( "Error: wrong number of arguments" )
    print(usage)
