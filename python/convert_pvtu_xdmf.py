# trace generated using paraview version 5.5.2

#### import the simple module from the paraview
from paraview.simple import *

def convert_p2x(step):
  #### disable automatic camera reset on 'Show'
  paraview.simple._DisableFirstRenderCameraReset()
  
  #step100/step000100_energy.pvts
  location = 'step' + str(step)
  in_fname = location + '/' + 'step' + ('%.6d' % step) + '_energy.pvts'
  out_fname = 'step' + ('%.6d' % step) + '_energy.xmf'
  print('<<== Converting',in_fname)

  # create a new 'XML Partitioned Structured Grid Reader'
  energypvts = XMLPartitionedStructuredGridReader(FileName=[in_fname])
  energypvts.CellArrayStatus = ['diffusivity_qp_avg', 'heatsource_qp_avg']
  energypvts.PointArrayStatus = ['temperature']
  
  # Properties modified on step000325_energypvts
  energypvts.CellArrayStatus = []

  # save data
  SaveData(out_fname, proxy=energypvts)



range = [ 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400 ]
for s in range:
  convert_p2x(s)


