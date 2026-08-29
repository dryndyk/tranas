"""This module defines TraNaS classes for Scanning Probe Microscopy (SPM).

Author: Dmitry A. Ryndyk <dmitry.ryndyk@tu-dresden.de>
"""


class SPM:
    """ TraNaS class for SPM.
    
    """    
    
    
    def __init__(self):
        
        self.atoms = None
        
    def shift(self, x,y,z):
        """Shifts the tip coordinates."""
        
        atoms = self.atoms
        for n in range(self.tip1-1,self.tip2):
            atoms.positions[n] = self.atoms.positions[n]+[x,y,z]
        
        return atoms
    
    def point(self, x,y,z):
        """Calculates STM for given tip coordinates."""
        
        # Change parameters of "atoms" object.   
        atoms = self.shift(x,y,z)
        
        # Make CP2K input file from "inp" and "atoms" objects.
        precalc = self.precalc
        precalc.atoms = atoms
        inp_file = precalc._generate_input()
        inp_fname = precalc.label + '.inp'
        file = open(inp_fname,'w')
        file.write(inp_file)
        file.close()
        from ase.io import write
        write(precalc.label+'.xyz', atoms)

        # Start CP2K calculation. 
        from subprocess import call, run
        #run(['/home/tranas/cp2k/install_psmp/bin/launch','mpirun', '-n', '8', 'cp2k', inp_fname])
        #run("/home/tranas/cp2k/install_psmp/bin/launch mpirun -n 8 cp2k "+inp_fname, shell=True)
        run("/home/tranas/cp2k/install_psmp/bin/launch mpirun -n 8 cp2k "+inp_fname, shell=True)

