from modelham import *

d = 1.0 # unit cell length

def linear(d):
    #define a 1D linear lattice cell with vectors a1 a2 and a3
    cell = Cell(1,[d,0.,0.],[0.,0.,0.],[0.,0.,0.])
    #add a site labeled 'A' at positon [0.,0.,0.]
    cell.add_site(LatticeSite([0.,0.,0.],'A'))
    #add a bond from first site in cell [0,0,0] to first site in cell [1,0,0]
    cell.add_bond(LatticeBond([0,0,0],0,[1,0,0],0))
    #buld a lattice of shape [4,1,1] with cell we defined
    lattice = Lattice(cell,[4,1,1])
    return lattice

lattice = linear(d)
#plot the lattice we constructed
lattice.plot() 
