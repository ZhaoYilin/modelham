from modelham import *
from modelham.sites import Site
from modelham.tn.mps import MPS
from modelham.tn.mpo import MPO
import numpy as np

site = SpinOneHalfSite()
### initial state |+-+-+-+-+->
Psi = MPS(site,24,4,'MPS')
### Hamiltonian MPO5
H = MPO(site,24,'Heisenberg',{'J':1.0,'Jz':1.0,'h':0.0})
Optimizer = Variational(Psi,H,50,100)
Optimizer.dmrg()
