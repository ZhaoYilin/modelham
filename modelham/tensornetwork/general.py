import numpy as np
from collections import Counter

from modelham.site import SpinHalfSite, ElectronSpinOrbital, ElectronSpatialOrbital
from modelham.symmetry import QuantumNumber

class Vertex(np.ndarray):
    """ TNS(U) := (V,E,H)
    V: set of vertices V= {U_i : i = 1,...,K}. K is the number of vertices.
    E: set of virtual indices dimension E = {m1,m2,...m_}. 
    H: set of physical indices H = {alpha_1,alpha_2,...,alpha_N}
    A is a 3-index tensor, A[s,i,j]

        s
        |
    i - A - j
        
    [s] acts on the local Hilbert space
    [i,j] act on the virtual bonds

    Parameters
    ----------
    np : _type_
        _description_
    """
    def __new__(cls, tensor, edges=None, half_edges=None):
        obj = np.asarray(tensor).view(cls)
        if edges and half_edges is not None and not isinstance((edges,half_edges), list):
            edges = list(edges)
            half_edges = list(half_edges)
        obj.edges = edges
        obj.half_edges = half_edges
        return obj
        
    def __init__(self):
        return None

    @classmethod
    def random(cls, indices):
        if [isinstance(indice, int) for indice in indices]:
            pass
        elif [isinstance(indice, Edge) for indice in indices]:
            indices = [indice.dimension for indice in indices]
        else:
            raise TypeError("Indices must be a list of int or Inidce.")
        tensor = super.random.rand(indices)
        obj = cls(tensor)
        return obj
    
    def selection(self, indices):
        """ In accordance to the selection rule mi−1 = si + mi one
        obtains the block structure of MPS matrices M[i],si as
        follows:

        Parameters
        ----------
        virtual_indices : _type_
            _description_
        physcial_indices : _type_
            _description_
        """        
        
        pass
 
class Edge(Counter):
    """ The symmetry in tensor networks is treated so that each of the indices is 
    decomposed into a quantum number (spin projection in this case) and its degeneracy index.

    MPS(U) := (V,E,H)
    V: set of vertices V= {U_i : i = 1,...,K}. K is the number of vertices.
    E: set of virtual indices dimension E = {m1,m2,...m_}. 
    H: set of physical indices H = {alpha_1,alpha_2,...,alpha_N}

    Attributes
    ----------
    self : Counter
        Dict of quantum label and number of indices.
    """    
    def __init__(self, *args, **kwargs):
        """_summary_

        Parameters
        ----------
        args : List[QuantumNumber]
            List of quantum numbers.
        """        
        super(Edge, self).__init__(*args, **kwargs)

    def __and__(self, other):
        """Intersection."""
        return self.__class__(super().__and__(other))

    def __or__(self, other):
        """Union."""
        return self.__class__(super().__or__(other))

    def __add__(self, other):
        """Sum of number of states."""
        return self.__class__(super().__add__(other))

    def __neg__(self):
        """Negation of keys."""
        return self.__class__({-k: v for k, v in self.items()})

    def __mul__(self, other):
        """Tensor product.
        """
        quanta = self.__class__()
        for ka, va in self.items():
            for kb, vb in other.items():
                quanta[ka + kb] += va * vb
        return quanta

    def __repr__(self):
        return " ".join(["%r = %d" % (k, v) for k, v in sorted(self.items(), key=lambda x: x[0])])

    @property
    def cardinality(self):
        """Total number of independent quantum number."""
        return len(self)
    
    @property
    def dimension(self):
        """Total number of quantum number."""
        return sum(self.values())

    @classmethod
    def virtual(cls, source, target, *sites):
        """_summary_

        Parameters
        ----------
        source : Indice
            _description_
        target : Indice
            _description_

        Returns
        -------
        _type_
            _description_
        """
        N = len(sites)
        physical_indices = [HalfEdge.physical(site) for site in sites]
        source_to_target = [source]
        target_to_source = [target]
        for i in range(N):
            source_to_target.append(source_to_target[-1] * physical_indices[i])
            target_to_source.append(target_to_source[-1] * -physical_indices[i])

        # The intersection gives the needed degeneracy sets
        intersection = [li & ri for li,ri in zip(source_to_target, target_to_source[::-1])]
        for item in intersection[1:-2]:
            item = cls(item)
        return intersection
    
    def filter(self, other):
        return self.__class__({k: min(other[k], v)
                         for k, v in self.items() if k in other})

    def truncate(self, bond_dim, ref=None):
        n_total = self.counting
        if n_total > bond_dim:
            for k, v in self.items():
                self[k] = int(np.ceil(v * bond_dim / n_total + 0.1))
                if ref is not None:
                    self[k] = min(self[k], ref[k])

    def truncate_no_keep(self, bond_dim, ref=None):
        n_total = self.counting
        if n_total > bond_dim:
            for k, v in self.items():
                self[k] = int(np.round(v * bond_dim / n_total + 0.1))
                if ref is not None:
                    self[k] = min(self[k], ref[k])
            for k in list(self.keys()):
                if self[k] == 0:
                    del self[k]

    def keep_maximal(self):
        maxk = max(self)
        return self.__class__({maxk: self[maxk]})


class HalfEdge(Edge):
    """ The symmetry in tensor networks is treated so that each of the indices is 
    decomposed into a quantum number (spin projection in this case) and its degeneracy index.

    MPS(U) := (V,E,H)
    V: set of vertices V= {U_i : i = 1,...,K}. K is the number of vertices.
    E: set of virtual indices dimension E = {m1,m2,...m_}. 
    H: set of physical indices H = {alpha_1,alpha_2,...,alpha_N}

    Attributes
    ----------
    self : Counter
        Dict of quantum label and number of indices.
    """    
    def __init__(self, *args, **kwargs):
        """_summary_

        Parameters
        ----------
        args : List[QuantumNumber]
            List of quantum numbers.
        """        
        super(HalfEdge, self).__init__(*args, **kwargs)

    @classmethod
    def vacuum(cls):
        obj = cls()
        obj[QuantumNumber(0, 0, 0)]=1
        return obj

    @classmethod
    def target(cls, n, twos, pg):
        obj = cls()
        obj[QuantumNumber(n, twos, pg)]=1
        return obj

    @classmethod
    def physical(cls, site):
        """_summary_

        Parameters
        ----------
        site : Site
            Instance of Site.
        """ 
        obj = cls()
        if isinstance(site, SpinHalfSite):
            obj[QuantumNumber(1, -1, site.pg)] = 1
            obj[QuantumNumber(1, 1, site.pg)] = 1
        elif isinstance(site, ElectronSpatialOrbital):
            obj[QuantumNumber(0, 0, site.pg)] = 1
            obj[QuantumNumber(1, -1, site.pg)] = 1
            obj[QuantumNumber(1, 1, site.pg)] = 1
            obj[QuantumNumber(2, 0, site.pg)] = 1
        else:
            raise TypeError('Site must be instance of either SpinSite or ElectronSite.')
        return obj