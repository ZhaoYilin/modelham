import numpy as np
import copy

class MPO(object):
    """
    """
    def __init__(self,site,N,model='Heisenberg',P={'J':1.0,'Jz':1.0,'h':0.5}):
        """
        N: numer of site
        P: dic of all parameters for the model hamilton
           {'J':1.0,'Jz':1.0,'h':1.0}
        """
        self.site = site
        self.N=N 
        self.P=P
        if model=='Ising':
            self.vertices = copy.deepcopy(self.Ising()) 
        elif model=='XXX':
            self.vertices = copy.deepcopy(self.XXX())
        elif model=='XXZ' or model=='Heisenberg':
            self.vertices = copy.deepcopy(self.XXZ())
        elif model=='XYZ':
            self.vertices = copy.deepcopy(self.XYZ())
        elif model=='Hubbard':
            pass
        elif model=='t-J':
            pass

    def Ising(self):
        """
        Ising Hamiltonian
            H = JSzSz - hSx
        P: dic of all parameters for the model hamilton
           {'J':1.0,'h':1.0}
        """
        self.D = 3
        P = self.P
        N = self.N
        site = self.site
        Id = site.operators["Id"]
        Sx = site.operators["Sx"]
        Sz = site.operators["Sz"]

        vertices = []

        row = np.zeros((1,2,2,3))
        column = np.zeros((3,2,2,1))
        matrix = np.zeros((3,2,2,3))

        #Left-hand edge is 1x5 matrix
        #[[I,   J*Sz,  -h*Sx]]
        row[0,:, :,0] = Id[:, :]
        row[0,:, :,1] = P['J']*Sz[:, :]
        row[0,:, :,2] = -P['h']*Sx[:, :]

        #Right-hand edge is 5x1 matrix
        #[[-hSx],
        #[Sz],
        #[I]]
        column[0,:, :,0] = -P['h']*Sx[:, :]
        column[1,:, :,0] = Sz[:, :]
        column[2,:, :,0] = Id[:, :]

        #Bulk Hamiltonian MPO
        #[I,    J*Sz,   -h*Sx]
        #[Z,    Z,      Sz]
        #[Z,    Z,      I]
        matrix[0, :, :, :] = row[0, :]
        matrix[:, :, :, 2] = column[:, 0]

        vertices.append(row[:, :])
        for _ in range(1, N-1):
            vertices.append(matrix[:, :, :, :])
        vertices.append(column[:, :])

        return vertices

    def XXX(self):
        pass

    def XXZ(self):
        """
        XXZ Hamiltonian
            H = JSzSz - hSx
        Parameters:
            P: parameters of the model Hamiltonians for each term
               J,Jz,h

        """
        self.D =5
        P = self.P
        N = self.N
        site = self.site
        Id = site.operators["Id"]
        Sp = site.operators["Sp"]
        Sm = site.operators["Sm"]
        Sz = site.operators["Sz"]

        vertices = []
        
        row = np.zeros((1, 5, 2, 2))
        column = np.zeros((5, 1, 2, 2))
        matrix = np.zeros((5, 5, 2, 2))
        #Left-hand edge is 1x5 matrix
        #[[I,   0.5J*Sp,  0.5J*Sm, Jz*Sz,  h*Sz]]
        row[0, 0, :, :] = Id[:, :]
        row[0, 1, :, :] = P['J'] / 2. * Sp[:, :]
        row[0, 2, :, :] = P['J'] / 2. * Sm[:, :]
        row[0, 3, :, :] = P['Jz'] * Sz[:, :]
        row[0, 4, :, :] = - P['h'] * Sz[:, :]

        #Right-hand edge is 5x1 matrix
        #[[-h*Sz],
        #[Sm],
        #[Sp],
        #[Sz],
        #[I]]
        column[0, 0, :, :] = - P['h'] * Sz[:, :]
        column[1, 0, :, :] = Sm[:, :]
        column[2, 0, :, :] = Sp[:, :]
        column[3, 0, :, :] = Sz[:, :]
        column[4, 0, :, :] = Id[:, :]

        #Bulk Hamiltonian MPO
        #[I,    J*Sp,   J*Sm,   Jz*Sz,  -h*Sz]
        #[Z,    Z,      Z,      Z,      Sz]
        #[Z,    Z,      Z,      Z,      Sm]
        #[Z,    Z,      Z,      Z,      Sp]
        #[Z,    Z,      Z,      Z,      I]
        matrix[0, :, :, :] = row[0, :]
        matrix[:, 4, :, :] = column[:, 0]

        #The whole MPO
        vertices.append(row[:, :])
        for _ in range(1, N-1):
            vertices.append(matrix[:, :, :, :])
        vertices.append(column[:, :])

        return vertices


    def XYZ(self):
        pass

    def Hubbard(self):
        pass

    def tJ(self):
        pass

    def one_vertice_exp(self,i,t):
        pass

    def two_vertices_exp(self,i,j,t):
        pass

    def exp(self,t):
        pass


