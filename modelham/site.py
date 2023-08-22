import numpy as np

class Site(object):
    """A class for general single site
    
    Use this class to create a single site object. The site comes with identity
    operator for a given dimension. To build specific site, additional operators
    need be add with add_operator method.
    
    """
    def __init__(self, dim, pg=0):
        """Creates an empty site of dimension dim.

        Parameters
        ----------

        dim : int 
            Size of the Hilbert space for single site. The dimension must be at 
            least 1. A site of dim = 1 is trival which represents the vaccum

        operators : a dictionary of string and numpy array (with ndim = 2).
            Operators for the site.

        """
        super(Site, self).__init__()
        self.dim = dim
        self.pg = pg
        self.states = {}
        self.operators = { "Id" : np.eye(self.dim, self.dim) }
    
    def add_operator(self, operator_name, matrix=None):
        """Adds an operator to the site with zero matrix.
    
        Parameters
        ----------
        operator_name : string
            The operator name.
	
        """
        if matrix==None:
            self.operators[str(operator_name)] = np.zeros((self.dim, self.dim))
        else:
            self.operators[str(operator_name)] = matrix

    def add_state(self, state_name, vector=None):
        """Adds an state to the site with zero list.
    
        Parameters
        ----------
        operator_name : string
            The operator name.
	
        """
        if vector==None:
            self.states[str(state_name)] = np.zeros(self.dim)
        else:
            self.states[str(state_name)] = vector

class SpinHalfSite(Site):
    """A site for spin 1/2 models.
    
    Use this class for spin one-half sites. The Hilbert space is ordered
    such as the first state is the spin down, and the second state is the
    spin up. 

    Notes
    -----
    Postcondition : The site has already built-in the spin operators for
                    Sx, Sy, Sz, Sp, Sm.

    """
    def __init__(self, pg=0):
        """Creates the spin one-half site.

        Notes
        -----
        Postcond : the dimension is set to 2

        """
        super(SpinHalfSite, self).__init__(2, pg)
        # add the operators
        self.add_operator("Sx")
        self.add_operator("Sy")
        self.add_operator("Sz")
        self.add_operator("Sp")
        self.add_operator("Sm")
        # for clarity
        Sx = self.operators["Sx"]
        Sy = self.operators["Sy"]
        Sy = Sy.astype(complex)
        Sz = self.operators["Sz"]
        Sp = self.operators["Sp"]
        Sm = self.operators["Sm"]
        # set the matrix elements different from zero to the right values
        Sx[0, 1] = 0.5
        Sx[1, 0] = 0.5
        Sy[0, 1] = 1j*(-0.5)
        Sy[1, 0] = 1j*0.5
        Sz[0, 0] = -0.5
        Sz[1, 1] = 0.5
        Sp[1, 0] = 1.0
        Sm[0, 1] = 1.0
        # add the states
        self.add_state("up")
        self.add_state("down")
        # for clarity
        up = self.states["up"]
        down = self.states["down"]
        # set the list elements different from zero to the right values
        up[1] = 1.0
        down[0] = 1.0

class ElectronSpinOrbital(Site):
    """A site for electronic models
    
    You use this site for models where the single sites are electron
    sites. The Hilbert space is ordered such as:

    - the first state, labelled 0,  is the empty site,
    - the second, labelled 1, is spin down, 
    - the third, labelled 2, is spin up, and 
    
    Notes
    -----
    Postcond: The site has already built-in the spin operators for: 

    - c_up : destroys an spin up electron,
    - c_up_dag, creates an spin up electron,
    - c_down, destroys an spin down electron,
    - c_down_dag, creates an spin down electron,
    - s_z, component z of spin,
    - s_p, raises the component z of spin,
    - s_m, lowers the component z of spin,
    - n_up, number of electrons with spin up,
    - n_down, number of electrons with spin down,
    - n, number of electrons, i.e. n_up+n_down, and
    - u, number of double occupancies, i.e. n_up*n_down.

    """
    def __init__(self, pg=0):
        super(ElectronSpinOrbital, self).__init__(4,pg)
        # add the operators
        self.add_operator("c_up")
        self.add_operator("c_up_dag")
        self.add_operator("c_down")
        self.add_operator("c_down_dag")
        self.add_operator("s_z")
        self.add_operator("n_up")
        self.add_operator("n_down")
        self.add_operator("n")
        self.add_operator("u")
        # for clarity
        c_up = self.operators["c_up"]
        c_up_dag = self.operators["c_up_dag"]
        c_down = self.operators["c_down"]
        c_down_dag = self.operators["c_down_dag"]
        s_z = self.operators["s_z"]
        n_up = self.operators["n_up"]
        n_down = self.operators["n_down"]
        n = self.operators["n"]
        u = self.operators["u"]
        # set the matrix elements different from zero to the right values
        c_up[0,2] = 1.0
        c_up[1,3] = 1.0
        c_up_dag[2,0] = 1.0
        c_up_dag[3,1] = 1.0
        c_down[0,1] = 1.0
        c_down[2,3] = 1.0
        c_down_dag[1,0] = 1.0
        c_down_dag[3,2] = 1.0
        s_z[1,1] = -1.0
        s_z[2,2] = 1.0
        n_up[2,2] = 1.0
        n_up[3,3] = 1.0
        n_down[1,1] = 1.0
        n_down[3,3] = 1.0
        n[1,1] = 1.0
        n[2,2] = 1.0
        n[3,3] = 2.0
        u[3,3] = 1.0
        # add the states
        self.add_state("empty")
        self.add_state("spin_down")
        self.add_state("spin_up")
        self.add_state("double")
        # for clarity
        state_empty = self.states["empty"]
        state_down = self.states["spin_down"]
        state_up = self.states["spin_up"]
        state_double = self.states["double"]
        # set the list elements different from zero to the right values
        state_empty[0] = 1.0
        state_down[1] = 1.0
        state_up[2] = 1.0
        state_double[3] = 1.0

class ElectronSpatialOrbital(Site):
    """A site for electronic models
    
    You use this site for models where the single sites are electron
    sites. The Hilbert space is ordered such as:

    - the first state, labelled 0,  is the empty site,
    - the second, labelled 1, is spin down, 
    - the third, labelled 2, is spin up, and 
    - the fourth, labelled 3, is double occupancy.
    
    Notes
    -----
    Postcond: The site has already built-in the spin operators for: 

    - c_up : destroys an spin up electron,
    - c_up_dag, creates an spin up electron,
    - c_down, destroys an spin down electron,
    - c_down_dag, creates an spin down electron,
    - s_z, component z of spin,
    - s_p, raises the component z of spin,
    - s_m, lowers the component z of spin,
    - n_up, number of electrons with spin up,
    - n_down, number of electrons with spin down,
    - n, number of electrons, i.e. n_up+n_down, and
    - u, number of double occupancies, i.e. n_up*n_down.

    """
    def __init__(self, pg=0):
        super(ElectronSpatialOrbital, self).__init__(4,pg)
        # add the operators
        self.add_operator("c_up")
        self.add_operator("c_up_dag")
        self.add_operator("c_down")
        self.add_operator("c_down_dag")
        self.add_operator("s_z")
        self.add_operator("n_up")
        self.add_operator("n_down")
        self.add_operator("n")
        self.add_operator("u")
        # for clarity
        c_up = self.operators["c_up"]
        c_up_dag = self.operators["c_up_dag"]
        c_down = self.operators["c_down"]
        c_down_dag = self.operators["c_down_dag"]
        s_z = self.operators["s_z"]
        n_up = self.operators["n_up"]
        n_down = self.operators["n_down"]
        n = self.operators["n"]
        u = self.operators["u"]
        # set the matrix elements different from zero to the right values
        c_up[0,2] = 1.0
        c_up[1,3] = 1.0
        c_up_dag[2,0] = 1.0
        c_up_dag[3,1] = 1.0
        c_down[0,1] = 1.0
        c_down[2,3] = 1.0
        c_down_dag[1,0] = 1.0
        c_down_dag[3,2] = 1.0
        s_z[1,1] = -1.0
        s_z[2,2] = 1.0
        n_up[2,2] = 1.0
        n_up[3,3] = 1.0
        n_down[1,1] = 1.0
        n_down[3,3] = 1.0
        n[1,1] = 1.0
        n[2,2] = 1.0
        n[3,3] = 2.0
        u[3,3] = 1.0
        # add the states
        self.add_state("empty")
        self.add_state("spin_down")
        self.add_state("spin_up")
        self.add_state("double")
        # for clarity
        state_empty = self.states["empty"]
        state_down = self.states["spin_down"]
        state_up = self.states["spin_up"]
        state_double = self.states["double"]
        # set the list elements different from zero to the right values
        state_empty[0] = 1.0
        state_down[1] = 1.0
        state_up[2] = 1.0
        state_double[3] = 1.0
