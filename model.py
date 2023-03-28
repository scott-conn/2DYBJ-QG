#This code solves the 2D coupled YBJ-QG system
#To run supply the following to the function run_sim:
#Ls = vector containing length of x and y dimensions
#ns = vector containing number of grid points in x and y dimensions
#p = vector of parameters - 1) hyperviscosity coefficient for q
#                           2) hyperviscosity coefficient for phi
#                           3) Coriolis parameter
#                           4) Stratification
#                           5) NIW vertical wavenumber
#                           6) Eddy scale
#                           7) Eddy velocity scale
#                           8) Number of eddy timescales to simulate
#                           9) Number of timesteps after which to write results
#                           10) Fraction of eddy timescale to take as timestep
#init_file which contains an initial wave and streamfunction field as well as the forcing timeseries.

import numpy as np
import dedalus.public as de
from dedalus.extras import flow_tools
import scipy.fft as spfft
from matplotlib import pyplot as plt
import time
import logging
import h5py
root = logging.root
for h in root.handlers:
    h.setLevel("INFO")  
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
logger = logging.getLogger(__name__)


def run_sim(Ls,ns,p,init_file):
  #create periodic dimensions
  Lx, Ly = Ls                                                      
  nx, ny = ns                                                        
  Tx, Ty = (Lx/nx,Ly/ny)
  x_basis = de.Fourier('x', nx, interval=(-Lx/2, Lx/2), dealias=3/2)        
  y_basis = de.Fourier('y', ny, interval=(-Ly/2, Ly/2), dealias=3/2)    

  #add x and y to a domain and then define an IVP with variables q, ψ and φ
  domain = de.Domain([x_basis, y_basis], grid_dtype=np.complex128)
  problem = de.IVP(domain, variables=['q','psi','phi'])
  
  #problem parameters
  kappa = p[0]
  nu = p[1]
  f0 = p[2]
  N0 = p[3]
  m0 = p[4]
  k_eddy = p[5]
  U_e = p[6]
  Nt = p[7]
  Nw = p[8]
  time_step = p[9]
  t_eddy = 1/U_e/k_eddy
  lambda_0 = N0/f0/m0
  hf = h5py.File(init_file, 'r')
  F = hf.get('F')[:]

  #add parameters to problem
  problem.parameters['kappa'] = kappa
  problem.parameters['nu'] = nu
  problem.parameters['f0'] = f0
  problem.parameters['lambda'] = lambda_0
  problem.parameters['eta'] = f0*lambda_0**2
  problem.parameters['F'] = F[0]
  
  #define shorthands for various differential operators, magnitude of a complex number and the wave-induced PV
  problem.substitutions["mag2(f)"] = "f * conj(f)"
  problem.substitutions["J(f,g)"] = "dx(f)*dy(g)-dy(f)*dx(g)"
  problem.substitutions["L(f)"] = "d(f,x=2) + d(f,y=2)"
  problem.substitutions["HD(f)"] = "L(L(f))"
  problem.substitutions["qw"] = "L(mag2(phi))/4/f0 + 1j*J(conj(phi),phi)/2/f0"
  
  #add model equations (note that for the k=l=0 mode we have a degeneracy which we remove by setting the gauge on ψ)
  problem.add_equation("dt(q) + kappa*HD(q) = -J(psi,q)")
  problem.add_equation("dt(phi) - 1j*eta*L(phi)/2 + nu*HD(phi) = -J(psi,phi) - 1j*phi*L(psi)/2 + F ")
  problem.add_equation("q - L(psi) = qw", condition="(nx != 0) or (ny != 0)")
  problem.add_equation("psi = 0.0 ", condition="(nx == 0) and (ny == 0)")
  
  #build solver
  solver = problem.build_solver('RK222')
  
  #specify simulation time
  solver.stop_sim_time = Nt*t_eddy
  solver.stop_wall_time = np.inf
  solver.stop_iteration = np.inf
  
  #read ψ and φ from the input file
  psi = solver.state['psi']
  psi.set_scales(1)
  slices = domain.dist.grid_layout.slices(scales=1)                                           
  psi_data = hf.get('psi')                                                                   
  psi['g'] = psi_data[slices]
  
  phi = solver.state['phi']
  phi.set_scales(1)
  slices = domain.dist.grid_layout.slices(scales=1)                                           
  phi_data = hf.get('phi')                                                                   
  phi['g'] = phi_data[slices]
  hf.close()

  #calculate qw and L(phi) to get q
  qw1 = ((np.conj(phi)*phi).evaluate().differentiate(x=2,y=2)/4/f0).evaluate()
  qw2 = (1j*((np.conj(phi).evaluate()).differentiate(x=1)*phi.differentiate(y=1) - 
             (np.conj(phi).evaluate().differentiate(y=1)*phi.differentiate(x=1))/2/f0)).evaluate()
  qw = (qw1 + qw2).evaluate()
  L = (psi.differentiate(x=2) + psi.differentiate(y=2)).evaluate()
  qw.set_scales(1)
  L.set_scales(1)
  q = solver.state['q']
  q.set_scales(1)
  q['g'] = L['g'] + qw['g']  
  
  #tell dedalus what to store in the output
  state_file = "state" #output file name
  analysis = solver.evaluator.add_file_handler(state_file, iter=Nw)
  analysis.add_system(solver.state, layout='g')
  analysis.add_task('L(psi)',name='zeta')
  analysis.add_task('qw',name='qw')

  #solve IVP
  dt = time_step*t_eddy
  logger.info('Starting loop')
  start_time = time.time()
  while solver.ok:
      solver.step(dt)
      #update forcing
      problem.namespace['F'].value = F[solver.iteration]
      if solver.iteration % 100 == 0:
          # Update plot of scalar field
          display.clear_output()
          logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))

  end_time = time.time()
  hf.close()

  display.clear_output()
  # Print statistics
  logger.info('Run time: %f' %(end_time-start_time))
  logger.info('Iterations: %i' %solver.iteration)

  #cleanup output file
  from dedalus.tools import post
  post.merge_process_files(state_file, cleanup=True)
