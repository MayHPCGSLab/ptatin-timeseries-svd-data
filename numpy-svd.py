
import os, sys
import numpy as np
import time

try:
  import PetscBinaryIO as pio # This file has been copied from the petsc source tree
except:
  # Hard wire path before importing PetscBinaryIO
  cwd = os.getcwd()
  ppath = os.getenv("PYTHONPATH") # None
  if ppath is None:
    os.environ["PYTHONPATH"] = os.path.join(cwd,"python")
  else:
    os.environ["PYTHONPATH"] += os.path.join(cwd,"python")
  sys.path.append(os.path.join(cwd,"python"))

  import PetscBinaryIO as pio # This file has been copied from the petsc source tree


def load_petsc_vec(fname):
  io = pio.PetscBinaryIO() # Instantiate a petsc binary loader
  with open(fname) as fp:
    objecttype = io.readObjectType(fp)
    v = io.readVec(fp)
  return v

steps = [ 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400 ]
#steps = [ 200, 300, 400 ]

n = len(steps)
print('n',n)

cnt = 0
for s in steps:
  filename = "data/step" + ('%1.6d' % s)  + "_energy.pbvec"
  field = load_petsc_vec(filename)
  if cnt == 0:
    m = len(field)
    print('m',m)
    U = np.zeros((m,n))
  U[:,cnt] = field
  cnt += 1

Umean = U.mean(1)
for i in range(n):
  U[:,i] = U[:,i] - Umean

t0 = time.clock()
u, s, vh = np.linalg.svd(U, full_matrices=False)
t1 = time.clock()
dt = float(t1-t0)
print('np.linalg.svd: ' + ('%1.4e' % dt) + ' (sec)')

print('singular values:\n',s)

# Build basis vectors
phi_normalized = np.zeros((m,n))
for i in range(n):
  phi = u[:,i] * s[i]
  norm = np.linalg.norm(phi)
  phi_normalized[:,i] = phi / norm
