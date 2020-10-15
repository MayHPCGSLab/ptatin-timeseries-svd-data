
import h5py as h5
import numpy as np
import PetscBinaryIO as pio # This file has been copied from the petsc source tree

def fetch_temperature(fname):
  fin = h5.File(fname,'r')

  dset = fin["Block_0_t000000/Node/temperature"]
  field = np.array(dset,'double')
  print('  *** loaded \"temperature\" ***')
  fin.close()
  return field

def fetch_coordinates(fname):
  fin = h5.File(fname,'r')
  
  dset = fin["Block_0_t000000/Geometry/Points"]
  field = np.array(dset,'double')
  print('  *** loaded \"coordinates\" ***')
  fin.close()
  return field

def write_as_petsc_vec(field,fname):
  io = pio.PetscBinaryIO()
  vec = field.view(pio.Vec)
  io.writeBinaryFile(fname, [vec,])


#filename = "step000200_energy.h5"
#field = fetch_temperature(filename)
#print("temperature [h5]",field.shape)
#f2 = field.flatten()
#print("temperature [vec]",f2.shape)

#stackfile = "svec.pbvec"
#io = pio.PetscBinaryIO()
#vec = f2.view(pio.Vec)
#io.writeBinaryFile(stackfile, [vec,])

#filename = "step000200_energy.h5"
#_field = fetch_temperature(filename)
#field = _field.flatten()
#_field = []
#write_as_petsc_vec(field,"svec.pbvec")


steps = [ 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400 ]
for s in steps:
  filename = "step" + ('%1.6d' % s)  + "_energy.h5"
  _field = fetch_temperature(filename)
  field = _field.flatten()
  _field = []
  ofilename = "step" + ('%1.6d' % s)  + "_energy.pbvec"
  write_as_petsc_vec(field,ofilename)

  _field = fetch_coordinates(filename)
  field = _field.flatten()
  _field = []
  ofilename = "step" + ('%1.6d' % s)  + "_coor.pbvec"
  write_as_petsc_vec(field,ofilename)

