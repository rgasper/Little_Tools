#!/usr/bin/env python
# Purpose: Generate HCP (001) surface
# Author: Raymond Gasper
# Adapted From FCC surface generator by Sungho Kim and Laalitha Liyanage

import os
import sys
import math

usage="""
        Usage: ./hcp_altlyrs_gen.py a c_a nx ny nz 

        Mandatory arguments
        -------------------
        a - equilibrium lattice constant
        c_a - c/a ratio

        Optional arguments
        -------------------
        nx,ny,nz - periodicity of supercell; DEFAULT (3,3,3)

        MUST provide either none or all of the optional arguments
"""


#Default setting
#--------------------------------------------------------------------
vacuum = 0.0
adatom = 0

#--------------------------Surface (001)-----------------------------

def gen_data_for_001_hcp(a, c_a, nx=3, ny=3, nz=3):
    la=[]; xa = [] ; ya = [] ; za = []
    ax = a
    ay = a * math.sqrt(3)
    az = a * c_a
    x0 = 0.0
    x2 = a * 1/2
    y2 = a * math.sqrt(3)/2
    y3 = a * math.sqrt(3)/3
    y4 = a * math.sqrt(3)*5/6
    bx, by, bz = ax*nx, ay*ny, az*nz
    for i in range(nx):
        for j in range(ny):
            layer = 0
            for k in range(nz):
                la.append(1); xa.append(x0+i*ax); ya.append(x0+j*ay); za.append((layer/2.0)*az)
                la.append(1); xa.append(x2+i*ax); ya.append(y2+j*ay); za.append((layer/2.0)*az); layer += 1
                la.append(2); xa.append(x0+i*ax); ya.append(y3+j*ay); za.append((layer/2.0)*az)
                la.append(2); xa.append(x2+i*ax); ya.append(y4+j*ay); za.append((layer/2.0)*az); layer += 1
    return la,xa,ya,za,bx,by,bz
    return None
#----------------------------lammps generation------------------------------------------------
def gen_lammps(la,xa,ya,za,bx,by,bz):
	fout = open("config.lammps","w")
	fout.write("Title Goes Here\n")
	fout.write("\n")
	fout.write("%d atoms\n"%len(xa))
	fout.write("\n")
	fout.write("2 atom types\n")
	fout.write("\n")
	fout.write(" %22.16f %22.16f xlo xhi\n"%(0,bx))
	fout.write(" %22.16f %22.16f ylo yhi\n"%(0,by))
	fout.write(" %22.16f %22.16f zlo zhi\n"%(0,bz))
	fout.write("\n")
	fout.write("Atoms \n")
	fout.write("\n")
	for i in range(len(xa)):
		index=i+1
		fout.write("%d %d %22.16f %22.16f %22.16f\n"%(index,la[i],xa[i],ya[i],za[i]))
	fout.close()
	return None
	
#-------------------------------Main program---------------------------------------------------
inputs=sys.argv[:]
if inputs[0] == 'python':
    del(inputs[0])

if len(inputs) == 3:
    a_latt = float(inputs[1])
    c_a_ratio = float(inputs[2])
    la,xa,ya,za,bx,by,bz = gen_data_for_001_hcp(a_latt,c_a_ratio)
    gen_lammps(la,xa,ya,za,bx,by,bz)

elif len(inputs) == 6:
    a_latt = float(inputs[1])
    c_a_ratio = float(inputs[2])
    nx = int(inputs[3])
    ny = int(inputs[4])
    nz = int(inputs[5])
    la,xa,ya,za,bx,by,bz = gen_data_for_001_hcp(a_latt,c_a_ratio,nx,ny,nz)
    gen_lammps(la,xa,ya,za,bx,by,bz)

else:
    print( "Error: wrong number of arguments" )
    print(usage)


