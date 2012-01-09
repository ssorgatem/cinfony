## Copyright (c) 2008-2011, Noel O'Boyle; 2012 Adrià Cereto-Massagué
## All rights reserved.
##
##  This file is part of Cinfony.
##  The contents are covered by the terms of the BSD license
##  which is included in the file LICENSE_BSD.txt.

"""
template - A Cinfony module template

Global variables:
  informats - a dictionary of supported input formats
  outformats - a dictionary of supported output formats
  descs - a list of supported descriptors
  fps - a list of supported fingerprint types
  forcefields - a list of supported forcefields
"""

import os.path

fps = []
"""A list of supported fingerprint types"""
descs = []
"""A list of supported descriptors"""

informats = []
"""A dictionary of supported input formats"""
outformats = []
"""A dictionary of supported output formats"""
forcefields = []
"""A list of supported forcefields"""


def readfile(format, filename):
    """Iterate over the molecules in a file.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       filename

    You can access the first molecule in a file using the next() method
    of the iterator:
        mol = readfile("smi", "myfile.smi").next()

    You can make a list of the molecules in a file using:
        mols = list(readfile("smi", "myfile.smi"))

    You can iterate over the molecules in a file as shown in the
    following code snippet:
    >>> atomtotal = 0
    >>> for mol in readfile("sdf", "head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    raise NotImplementedError

def readstring(format, string):
    """Read in a molecule from a string.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       string

    Example:
    >>> input = "C1=CC=CS1"
    >>> mymol = readstring("smi", input)
    >>> len(mymol.atoms)
    5
    """
    raise NotImplementedError


class Outputfile(object):
    """Represent a file to which *output* is to be sent.

    Although it's possible to write a single molecule to a file by
    calling the write() method of a molecule, if multiple molecules
    are to be written to the same file you should use the Outputfile
    class.

    Required parameters:
       format - see the outformats variable for a list of available
                output formats
       filename

    Optional parameters:
       overwrite -- if the output file already exists, should it
                   be overwritten? (default is False)

    Methods:
       write(molecule)
       close()
    """
    def __init__(self, format, filename, overwrite=False):
        self.format = format
        self.filename = filename
        if not overwrite and os.path.isfile(self.filename):
            raise IOError("%s already exists. Use 'overwrite=True' to overwrite it." % self.filename)
        self.total = 0 # The total number of molecules written to the file

    def write(self, molecule):
        """Write a molecule to the output file.

        Required parameters:
           molecule
        """
        if not self.filename:
            raise IOError("Outputfile instance is closed.")

        raise NotImplementedError

        self.total += 1

    def close(self):
        """Close the Outputfile to further writing."""
        self.filename = None
        raise NotImplementedError

class Molecule(object):
    """Represent a Pybel Molecule.

    Required parameter:
       Mol -- a backend molecule or any type of cinfony Molecule

    Attributes:
       atoms, charge, conformers, data, dim, energy, exactmass, formula,
       molwt, spin, sssr, title, unitcell.
    (refer to the Open Babel library documentation for more info).

    Methods:
       addh(), calcfp(), calcdesc(), draw(), localopt(), make3D(), removeh(),
       write()

    The underlying backend molecule can be accessed using the attribute:
       Mol
    """
    _cinfony = True

    def __init__(self, Mol):

        if hasattr(Mol, "_cinfony"):
            a, b = Mol._exchange
            if a == 0:
                molecule = readstring("smi", b)
            else:
                molecule = readstring("mol", b)
            Mol = molecule.Mol

        self.Mol = Mol

    @property
    def atoms(self): raise NotImplementedError
    @property
    def charge(self): raise NotImplementedError
    @property
    def conformers(self): raise NotImplementedError
    @property
    def data(self): raise NotImplementedError
    @property
    def dim(self): raise NotImplementedError
    @property
    def energy(self): raise NotImplementedError
    @property
    def exactmass(self): raise NotImplementedError
    @property
    def formula(self): raise NotImplementedError
    @property
    def molwt(self): raise NotImplementedError
    @property
    def spin(self): raise NotImplementedError
    @property
    def sssr(self): return self.OBMol.GetSSSR()
    def _gettitle(self): raise NotImplementedError
    def _settitle(self, val): raise NotImplementedError
    title = property(_gettitle, _settitle)
    @property
    def unitcell(self):
        raise NotImplementedError
    @property
    def _exchange(self):
        raise NotImplementedError
#        if self.Mol has coordinate data:
#            return (1, self.write("mol"))
#        else: #If it has no coordinate data:
#            return (0, self.write("can").split()[0])

    def __iter__(self):
        """Iterate over the Atoms of the Molecule.

        This allows constructions such as the following:
           for atom in mymol:
               print atom
        """
        return iter(self.atoms)

    def calcdesc(self, descnames=[]):
        """Calculate descriptor values.

        Optional parameter:
           descnames -- a list of names of descriptors

        If descnames is not specified, all available descriptors are
        calculated. See the descs variable for a list of available
        descriptors.
        """
        if not descnames:
            descnames = descs
        ans = {}
        for descname in descnames:
            raise NotImplementedError
        return ans

    def calcfp(self, fptype="FP2"):
        """Calculate a molecular fingerprint.

        Optional parameters:
           fptype -- the fingerprint type (default is "FP2"). See the
                     fps variable for a list of of available fingerprint
                     types.
        """
        raise NotImplementedError
        return Fingerprint(fp)

    def write(self, format="smi", filename=None, overwrite=False):
        """Write the molecule to a file or return a string.

        Optional parameters:
           format -- see the informats variable for a list of available
                     output formats (default is "smi")
           filename -- default is None
           overwite -- if the output file already exists, should it
                       be overwritten? (default is False)

        If a filename is specified, the result is written to a file.
        Otherwise, a string is returned containing the result.

        To write multiple molecules to the same file you should use
        the Outputfile class.
        """
        raise NotImplementedError

    def localopt(self, forcefield='', steps=500):
        """Locally optimize the coordinates.

        Optional parameters:
           forcefield -- default is "". See the forcefields variable
                         for a list of available forcefields.
           steps -- default is 500

        If the molecule does not have any coordinates, make3D() is
        called before the optimization. Note that the molecule needs
        to have explicit hydrogens. If not, call addh().
        """
        forcefield = forcefield.lower()
        if self.dim != 3:
            self.make3D(forcefield)
        raise NotImplementedError

    def make3D(self, forcefield = "", steps = 50):
        """Generate 3D coordinates.

        Optional parameters:
           forcefield -- default is "". See the forcefields variable
                         for a list of available forcefields.
           steps -- default is 50

        Once coordinates are generated, hydrogens are added and a quick
        local optimization is carried out with 50 steps and the
        MMFF94 forcefield. Call localopt() if you want
        to improve the coordinates further.
        """
        forcefield = forcefield.lower()
        self.addh()
        self.localopt(forcefield, steps)

    def addh(self):
        """Add hydrogens."""
        raise NotImplementedError

    def removeh(self):
        """Remove hydrogens."""
        raise NotImplementedError

    def __str__(self):
        return self.write()

    def draw(self, show=True, filename=None, update=False, usecoords=False):
        """Create a 2D depiction of the molecule.

        Optional parameters:
          show -- display on screen (default is True)
          filename -- write to file (default is None)
          update -- update the coordinates of the atoms to those
                    determined by the structure diagram generator
                    (default is False)
          usecoords -- don't calculate 2D coordinates, just use
                       the current coordinates (default is False)

        """
        raise NotImplementedError

class Fingerprint(object):
    """A Molecular Fingerprint.

    Required parameters:
       fingerprint -- a string of 0's and 1's representing a binary fingerprint

    Attributes:
       fp -- the underlying fingerprint object
       bits -- a list of bits set in the Fingerprint

    Methods:
       The "|" operator can be used to calculate the Tanimoto coeff. For example,
       given two Fingerprints 'a', and 'b', the Tanimoto coefficient is given by:
          tanimoto = a | b
    """
    def __init__(self, fingerprint):
        self.fp = fingerprint
    def __or__(self, other):
        mybits = set(self.bits)
        otherbits = set(other.bits)
        return len(mybits&otherbits) / float(len(mybits|otherbits))
    @property
    def bits(self):
        raise NotImplementedError
    def __str__(self):
        return ", ".join([str(x) for x in self.fp])

def _findbits(fp, bitsperint):
    """Find which bits are set in a list/vector.

    This function is used by the Fingerprint class.

    >>> _findbits([13, 71], 8)
    [1, 3, 4, 9, 10, 11, 15]
    """
    ans = []
    start = 1
    for x in fp:
        i = start
        while x > 0:
            if x % 2:
                ans.append(i)
            x >>= 1
            i += 1
        start += bitsperint
    return ans

class Atom(object):
    """Represent a backend atom.

    Required parameter:
       Atom -- a backend Atom

    Attributes:
       atomicmass, atomicnum, cidx, coords, coordidx, exactmass,
       formalcharge, heavyvalence, heterovalence, hyb, idx,
       implicitvalence, isotope, partialcharge, spin, type,
       valence, vector.

    (refer to the backend documentation for more info).

    The original Atom can be accessed using the attribute:
       Atom
    """

    def __init__(self, Atom):
        self.Atom = Atom

    @property
    def coords(self):
        raise NotImplementedError
    @property
    def atomicmass(self): raise NotImplementedError
    @property
    def atomicnum(self): raise NotImplementedError
    @property
    def cidx(self): raise NotImplementedError
    @property
    def coordidx(self): raise NotImplementedError
    @property
    def exactmass(self): raise NotImplementedError
    @property
    def formalcharge(self): raise NotImplementedError
    @property
    def heavyvalence(self): raise NotImplementedError
    @property
    def heterovalence(self): raise NotImplementedError
    @property
    def hyb(self): raise NotImplementedError
    @property
    def idx(self): raise NotImplementedError
    @property
    def implicitvalence(self): raise NotImplementedError
    @property
    def isotope(self): raise NotImplementedError
    @property
    def partialcharge(self): raise NotImplementedError
    @property
    def spin(self): raise NotImplementedError
    @property
    def type(self): raise NotImplementedError
    @property
    def valence(self): raise NotImplementedError
    @property
    def vector(self): raise NotImplementedError

    def __str__(self):
        c = self.coords
        return "Atom: %d (%.2f %.2f %.2f)" % (self.atomicnum, c[0], c[1], c[2])


class Smarts(object):
    """A Smarts Pattern Matcher

    Required parameters:
       smartspattern

    Methods:
       findall(molecule)

    Example:
    >>> mol = readstring("smi","CCN(CC)CC") # triethylamine
    >>> smarts = Smarts("[#6][#6]") # Matches an ethyl group
    >>> print smarts.findall(mol)
    [(1, 2), (4, 5), (6, 7)]

    The numbers returned are the indices (starting from 1) of the atoms
    that match the SMARTS pattern. In this case, there are three matches
    for each of the three ethyl groups in the molecule.
    """
    def __init__(self,smartspattern):
        """Initialise with a SMARTS pattern."""
        raise NotImplementedError
    def findall(self,molecule):
        """Find all matches of the SMARTS pattern to a particular molecule.

        Required parameters:
           molecule
        """
        raise NotImplementedError
    def match(self, molecule):
        """Does a SMARTS pattern match a particular molecule?

        Required parameters:
           molecule
        """
        raise NotImplementedError

class MoleculeData(object):
    """Store molecule data in a dictionary-type object

    Required parameters:
      mol -- a backend Mol

    Methods and accessor methods are like those of a dictionary except
    that the data is retrieved on-the-fly from the underlying Mol.

    Example:
    >>> mol = readfile("sdf", 'head.sdf').next()
    >>> data = mol.data
    >>> print data
    {'Comment': 'CORINA 2.61 0041  25.10.2001', 'NSC': '1'}
    >>> print len(data), data.keys(), data.has_key("NSC")
    2 ['Comment', 'NSC'] True
    >>> print data['Comment']
    CORINA 2.61 0041  25.10.2001
    >>> data['Comment'] = 'This is a new comment'
    >>> for k,v in data.iteritems():
    ...    print k, "-->", v
    Comment --> This is a new comment
    NSC --> 1
    >>> del data['NSC']
    >>> print len(data), data.keys(), data.has_key("NSC")
    1 ['Comment'] False
    """
    def __init__(self, mol):
        self._mol = mol
    def _data(self):
        raise NotImplementedError
    def _testforkey(self, key):
        if not key in self:
            raise KeyError("'%s'" % key)
    def keys(self):
        raise NotImplementedError
    def values(self):
        raise NotImplementedError
    def items(self):
        return zip(self.keys(), self.values())
    def __iter__(self):
        return iter(self.keys())
    def iteritems(self):
        return iter(self.items())
    def __len__(self):
        return len(self._data())
    def __contains__(self, key):
        raise NotImplementedError
    def __delitem__(self, key):
        self._testforkey(key)
        raise NotImplementedError
    def clear(self):
        for key in self:
            del self[key]
    def has_key(self, key):
        return key in self
    def update(self, dictionary):
        for k, v in dictionary.iteritems():
            self[k] = v
    def __getitem__(self, key):
        self._testforkey(key)
        raise NotImplementedError
    def __setitem__(self, key, value):
        raise NotImplementedError
    def __repr__(self):
        return dict(self.iteritems()).__repr__()

if __name__=="__main__": #pragma: no cover
    import doctest
    doctest.testmod(verbose=True)
