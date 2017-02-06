#!/usr/bin/env python
# Purpose: Generate FCC (111), (110), and (100) surface
# Author: Raymond Gasper
# Adapted from FCC surface generator by Sungho Kim and Laalitha Liyanage

import os
import sys
import math

usage="""
        Usage: ./fcc_surf_gen.py a surf out_style vacuum nx ny nz adatom

        Mandatory arguments
        -------------------
        a - equilibrium lattice constant
        surf - Type of surface 100, 110 or 111
        out_style - vasp or lammps

        Optional arguments
        -------------------
        vacuum - length of vacuum; DEFAULT = 0.0 angstroms
        nx,ny,nz - periodicity of supercell; DEFAULT (2,2,4)
        adatom - 1/0 (True/False); DEFAULT = 0 (False)

        MUST provide either none or all of the optional arguments
        will generate a bulk structure if you leave vacuum at zero angstroms
"""


#Default setting
#--------------------------------------------------------------------
vacuum = 0.0
adatom = 0

#--------------------------Surface (100)-----------------------------

def gen_data_for_100_fcc(a,nx=2,ny=2,nz=4):
  """ Generate datafile of FCC structure with lattice constant a """
  xa=[]; ya=[]; za=[]
  x0 = 0.0
  bx,by,bz = a*nx,a*ny,a*nz+vacuum
  x,y,z = bx,by,bz

  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(0 + i*a); ya.append(  0 + j*a); za.append(  0 + k*a)
        xa.append(  0 + i*a); ya.append(a/2 + j*a); za.append(a/2 + k*a)
        xa.append(a/2 + i*a); ya.append(  0 + j*a); za.append(a/2 + k*a)
        xa.append(a/2 + i*a); ya.append(a/2 + j*a); za.append(  0 + k*a)
  if adatom != 0:
#        xa.append(x0); ya.append(x0); za.append(x0+nz*a)
        xa.append(bx/2.); ya.append(by/2.); za.append(x0+nz*a)
  return xa,ya,za,bx,by,bz

#--------------------------Surface (110)-----------------------------

def gen_data_for_110_fcc(a,nx=2,ny=2,nz=4):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  xa=[]; ya=[]; za=[]
  ax = a*math.sqrt(2)/2
  ay = a*math.sqrt(6)/2
  az = a*math.sqrt(3)
  x0 = 0.0
  x2 = math.sqrt(2)/4. * a
  y2 = math.sqrt(6)/4. * a
  y3 = math.sqrt(6)/6. * a
  y4 = math.sqrt(6)*5./12. * a
  y5 = math.sqrt(6)*2./6. * a
  y6 = math.sqrt(6)/12 * a
  z3 = math.sqrt(3)/3. * a
  z5 = math.sqrt(3)*2./3. * a
  bx,by,bz = ax*nx + vacuum, ay*ny, az*nz
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(x0+k*az)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(x0+k*az)
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(z3+k*az)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(z3+k*az)
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(z5+k*az)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(z5+k*az)
  if adatom != 0:
        xa.append(x0+nx*ax); ya.append(by/2.); za.append(bz/2.)
  return xa,ya,za,bx,by,bz

#--------------------------Surface (111)-----------------------------

def gen_data_for_111_fcc(a,nx=2,ny=2,nz=4):
  """ Generate datafile of FCC surface: 110:x, 112:y, 111:z """
  xa=[]; ya=[]; za=[]
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
  bx,by,bz = ax*nx, ay*ny, az*nz+vacuum
  for i in range(nx):
    for j in range(ny):
      layer = 0
      for k in range(nz):
        xa.append(x0+i*ax); ya.append(x0+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y2+j*ay); za.append(layer/3.0*az); layer += 1
        xa.append(x0+i*ax); ya.append(y3+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y4+j*ay); za.append(layer/3.0*az); layer += 1
        xa.append(x0+i*ax); ya.append(y5+j*ay); za.append(layer/3.0*az)
        xa.append(x2+i*ax); ya.append(y6+j*ay); za.append(layer/3.0*az); layer += 1
  if adatom != 0:
        xa.append(bx/2.); ya.append(by/2.); za.append(x0+nz*az)
  return xa,ya,za,bx,by,bz
#----------------------------POSCAR generation------------------------------------------------
def gen_poscar(xa,ya,za,bx,by,bz):
  fout = open("config.vasp","w")
  fout.write("Title Goes Here\n")
  fout.write("1.0\n")
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(bx,0,0))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(0,by,0))
  fout.write(" %22.16f  %22.16f  %22.16f\n"%(0,0,bz))
  fout.write("%d\n"%len(xa))
#  fout.write("Selective Dynamics\n")
  fout.write("Cart\n")
  for i in range(len(xa)):
    fout.write("%22.16f %22.16f %22.16f\n"%(xa[i],ya[i],za[i]))
#    fout.write("%22.16f %22.16f %22.16f F F T\n"%(xa[i],ya[i],za[i]))
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
	fout.write("1 mass\n")
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
    surf = inputs[2]
    out_style = inputs[3]
    if surf == '100' :
        xa,ya,za,bx,by,bz = gen_data_for_100_fcc(a_latt)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
    elif surf == '110':
        xa,ya,za,bx,by,bz = gen_data_for_110_fcc(a_latt)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
    elif surf == '111':
        xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
        
elif len(inputs) == 9:
    a_latt = float(inputs[1])
    surf = inputs[2]
    out_style = inputs[3]
    vacuum = float(inputs[4])
    nx = int(inputs[5])
    ny = int(inputs[6])
    nz = int(inputs[7])
    adatom = int(inputs[8])
    
    if surf == '100' :
        xa,ya,za,bx,by,bz = gen_data_for_100_fcc(a_latt,nx,ny,nz)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
    
    elif surf == '110':
        xa,ya,za,bx,by,bz = gen_data_for_110_fcc(a_latt,nx,ny,nz)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
						
    elif surf == '111':	
        xa,ya,za,bx,by,bz = gen_data_for_111_fcc(a_latt,nx,ny,nz)
        if out_style == 'lammps':
            gen_lammps(xa,ya,za,bx,by,bz)
        elif out_style == 'vasp':
            gen_poscar(xa,ya,za,bx,by,bz)
else:
    print( "Error: wrong number of arguments" )
    print(usage)
