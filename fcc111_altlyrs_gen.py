#!/usr/bin/env python
# Purpose: Generate FCC (111) bulk with layers of alternating elements.
# Author: Raymond Gasper
# Adapted from FCC surface generator by Sungho Kim and Laalitha Liyanage

import os
import sys
import math

usage="""
        Usage: ./fcc111_altlyr_gen.py a nx ny nz

        Mandatory arguments
        -------------------
        a - equilibrium lattice constant

        Optional arguments
        -------------------
        nx,ny,nz - periodicity of supercell; DEFAULT (3,3,1)

        MUST provide either none or all of the optional arguments
	note: minimum size (nz=1; default) is 6 layers tall
"""


#--------------------------Surface (111)-----------------------------

def gen_data_for_111_fcc(a,nx=3,ny=3,nz=1):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  la=[]; xa=[]; ya=[]; za=[]
  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4 * a
  y2 = math.sqrt(6)/4 * a
  y3 = math.sqrt(6)/6 * a
  y4 = math.sqrt(6)*5/12 * a
  y5 = math.sqrt(6)*2/6 * a
  y6 = math.sqrt(6)/12 * a
  bx,by,bz = ax*nx, ay*ny, az*2*nz
  for i in range(nx):
    for j in range(ny):
      layer = 0
      for k in range(nz):
        la.append(1); xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)
        la.append(1); xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1
        la.append(2); xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)
        la.append(2); xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1
        la.append(1); xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)
        la.append(1); xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1
        la.append(2); xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)
        la.append(2); xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1
        la.append(1); xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)
        la.append(1); xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1
        la.append(2); xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)
        la.append(2); xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1
  return la,xa,ya,za,bx,by,bz
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
		fout.write("%d %d %22.16f %22.16f %22.16f\n"%(index, la[i], xa[i],ya[i],za[i]))
	fout.close()
	return None
	
#-------------------------------Main program---------------------------------------------------

inputs=sys.argv[:]
if inputs[0] == 'python':
    del(inputs[0])

if len(inputs) == 2:
    a_latt = float(inputs[1])
    la,xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt)
    gen_lammps(la,xa,ya,za,bx,by,bz)
elif len(inputs) == 5:
    a_latt = float(inputs[1])
    nx = int(inputs[2])
    ny = int(inputs[3])
    nz = int(inputs[4])
    la,xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt,nx,ny,nz)
    gen_lammps(la,xa,ya,za,bx,by,bz)
else:
    print( "Error: wrong number of arguments" )
    print(usage)
