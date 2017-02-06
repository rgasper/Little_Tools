#!/usr/bin/env python
# Purpose: Generate HCP (001) surface
# Author: Raymond Gasper
# Adapted From FCC surface generator by Sungho Kim and Laalitha Liyanage

import os
import sys
import math

usage="""
        Usage: ./hcp_surf_gen.py a c_a out_style vacuum nx ny nz adatom

        Mandatory arguments
        -------------------
        a - equilibrium lattice constant
        c_a - c/a ratio
        out_style - vasp or lammps

        Optional arguments
        -------------------
        vacuum - length of vacuum; DEFAULT = 0.0 angstroms
        nx,ny,nz - periodicity of supercell; DEFAULT (2,2,4)
        adatom - 1/0 (True/False); DEFAULT = 0 (False)

        MUST provide either none or all of the optional arguments
        Will generate a bulk structure if you ask for zero angstroms vacuum
"""


#Default setting
#--------------------------------------------------------------------
vacuum = 0.0
adatom = 0

#--------------------------Surface (001)-----------------------------

def gen_data_for_001_hcp(a, c_a, nx=2, ny=2, nz=4):
    xa = [] ; ya = [] ; za = []
    ax = a
    ay = a * math.sqrt(3)
    az = a * c_a
    x0 = 0.0
    x2 = a * 1/2
    y2 = a * math.sqrt(3)/2
    y3 = a * math.sqrt(3)/3
    y4 = a * math.sqrt(3)*5/6
    bx, by, bz = ax*nx, ay*ny, az*nz+vacuum
    for i in range(nx):
        for j in range(ny):
            layer = 0
            for k in range(nz):
                xa.append(x0+i*ax); ya.append(x0+j*ay); za.append((layer/2.0)*az)
                xa.append(x2+i*ax); ya.append(y2+j*ay); za.append((layer/2.0)*az); layer += 1
                xa.append(x0+i*ax); ya.append(y3+j*ay); za.append((layer/2.0)*az)
                xa.append(x2+i*ax); ya.append(y4+j*ay); za.append((layer/2.0)*az); layer += 1
    if adatom != 0:
        xa.append(bx/2.); ya.append(by/2.); za.append(x0+nz*az)
    return xa,ya,za,bx,by,bz
#----------------------------VASP POSCAR generation------------------------------------------------
def gen_poscar(xa, ya, za,bx,by,bz):
    fout = open("config.vasp", "w")
    fout.write("Title Goes Here\n")
    fout.write("1.0\n")
    fout.write(" %22.16f    %22.16f    %22.16f\n"%(bx,0,0))
    fout.write(" %22.16f    %22.16f    %22.16f\n"%(0,by,0))
    fout.write(" %22.16f    %22.16f    %22.16f\n"%(0,0,bz))
    fout.write("Ru\n")
    fout.write("%d\n"%len(xa))
#   fout.write("Selective Dynamics\n")
    fout.write("Cartesian\n")
    for i in range(len(xa)):
        fout.write("%22.16f %22.16f %22.16f\n"%(xa[i],ya[i],za[i]))
#       fout.write("%22.16f %22.16f %22.16f F F T\n"%(xa[i],ya[i],za[i]))
    fout.close()
    return None
#----------------------------lammps generation------------------------------------------------
def gen_lammps(xa,ya,za,bx,by,bz):
	fout = open("config.lammps","w")
	fout.write("Title Goes Here\n")
	fout.write("\n")
	fout.write("%d atoms\n"%len(xa))
	fout.write("\n")
	fout.write("1 atom types\n")
	fout.write("\n")
	fout.write(" %22.16f %22.16f xlo xhi\n"%(0,bx))
	fout.write(" %22.16f %22.16f ylo yhi\n"%(0,by))
	fout.write(" %22.16f %22.16f zlo zhi\n"%(0,bz))
	fout.write("\n")
	fout.write("Masses\n")
	fout.write("\n")
	fout.write("1 user_input_mass\n")  #### 101.07 is the value for Ruthenium, that's what I'm using for now
	fout.write("\n")
	fout.write("Atoms \n")
	fout.write("\n")
	for i in range(len(xa)):
		index=i+1
		fout.write("%d 1 %22.16f %22.16f %22.16f\n"%(index,xa[i],ya[i],za[i]))
	fout.close()
	return None
	
#-------------------------------Main program---------------------------------------------------
inputs=sys.argv[:]
if inputs[0] == 'python':
    del(inputs[0])

if len(inputs) == 4:
    a_latt = float(inputs[1])
    c_a_ratio = float(inputs[2])
    out_style = inputs[3]
    xa,ya,za,bx,by,bz = gen_data_for_001_hcp(a_latt,c_a_ratio)
    if out_style == 'lammps':
        gen_lammps(xa,ya,za,bx,by,bz)
    elif out_style == 'vasp':
        gen_poscar(xa,ya,za,bx,by,bz)

elif len(inputs) == 9:
    a_latt = float(inputs[1])
    c_a_ratio = float(inputs[2])
    out_style = inputs[3]
    vacuum = float(inputs[4])
    nx = int(inputs[5])
    ny = int(inputs[6])
    nz = int(inputs[7])
    adatom = int(inputs[8])
    xa,ya,za,bx,by,bz = gen_data_for_001_hcp(a_latt,c_a_ratio,nx,ny,nz)
    if out_style == 'lammps':
        gen_lammps(xa,ya,za,bx,by,bz)
    elif out_style == 'vasp':
        gen_poscar(xa,ya,za,bx,by,bz)

else:
    print( "Error: wrong number of arguments" )
    print(usage)


