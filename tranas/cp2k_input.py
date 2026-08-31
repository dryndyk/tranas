"""TraNaS python library
Author: Dmitry A. Ryndyk <dmitry.ryndyk@tu-dresden.de>

Input precalculator for CP2K (https://www.cp2k.org/)
Based on ASE-Calculator for CP2K
Author: Ole Schuett <ole.schuett@mat.ethz.ch>
"""

from warnings import warn
from ase.units import Rydberg
from ase.calculators.cp2k import CP2K, InputSection, parse_input


class CP2Kinput(CP2K):
    """TraNaS precalculator for CP2K.

    Arguments:

    auto_write: bool
        Flag to enable the auto-write mode. If enabled the
        ``write()`` routine is called after every
        calculation, which mimics the behavior of the
        ``FileIOCalculator``. Default is ``False``.
    basis_set: str
        Name of the basis set to be use.
        The default is ``DZVP-MOLOPT-SR-GTH``.
    basis_set_file: str
        Filename of the basis set file.
        Default is ``BASIS_MOLOPT``.
        Set the environment variable $CP2K_DATA_DIR
        to enabled automatic file discovered.
    charge: float
        The total charge of the system.  Default is ``0``.
    command: str
        The command used to launch the CP2K-shell.
        If ``command`` is not passed as an argument to the
        constructor, the class-variable ``CP2K.command``,
        and then the environment variable
        ``$ASE_CP2K_COMMAND`` are checked.
        Eventually, ``cp2k.psmp -s`` is used as default.
    cutoff: float
        The cutoff of the finest grid level.  Default is ``400 * Rydberg``.
    debug: bool
        Flag to enable debug mode. This will print all
        communication between the CP2K-shell and the
        CP2K-calculator. Default is ``False``.
    force_eval_method: str
        The method CP2K uses to evaluate energies and forces.
        The default is ``Quickstep``, which is CP2K's
        module for electronic structure methods like DFT.
    inp: str
        CP2K input template. If present, the calculator will
        augment the template, e.g. with coordinates, and use
        it to launch CP2K. Hence, this generic mechanism
        gives access to all features of CP2K.
        Note, that most keywords accept ``None`` to disable the generation
        of the corresponding input section.

        This input template is important for advanced CP2K
        inputs, but is also needed for e.g. controlling the Brillouin
        zone integration. The example below illustrates some common
        options::

            inp = '''&FORCE_EVAL
               &DFT
                 &KPOINTS
                   SCHEME MONKHORST-PACK 12 12 8
                 &END KPOINTS
                 &SCF
                   ADDED_MOS 10
                   &SMEAR
                     METHOD FERMI_DIRAC
                     ELECTRONIC_TEMPERATURE [K] 500.0
                   &END SMEAR
                 &END SCF
               &END DFT
             &END FORCE_EVAL
            '''

    max_scf: int
        Maximum number of SCF iteration to be performed for
        one optimization. Default is ``50``.
    multiplicity: int, default=None
        Select the multiplicity of the system
        (two times the total spin plus one).
        If None, multiplicity is not explicitly given in the input file.
    poisson_solver: str
        The poisson solver to be used. Currently, the only supported
        values are ``auto`` and ``None``. Default is ``auto``.
    potential_file: str
        Filename of the pseudo-potential file.
        Default is ``POTENTIAL``.
        Set the environment variable $CP2K_DATA_DIR
        to enabled automatic file discovered.
    pseudo_potential: str
        Name of the pseudo-potential to be use.
        Default is ``auto``. This tries to infer the
        potential from the employed XC-functional,
        otherwise it falls back to ``GTH-PBE``.
    stress_tensor: bool
        Indicates whether the analytic stress-tensor should be calculated.
        Default is ``True``.
    uks: bool
        Requests an unrestricted Kohn-Sham calculations.
        This is need for spin-polarized systems, ie. with an
        odd number of electrons. Default is ``False``.
    xc: str
        Name of exchange and correlation functional.
        Accepts all functions supported by CP2K itself or libxc.
        Default is ``LDA``.
    xtb: str
         Name of the xTB method. Currently, the supported values are 'GFN1' and 'GFN2'(default).
    print_level: str
        PRINT_LEVEL of global output.
        Possible options are:
        DEBUG Everything is written out, useful for debugging purposes only
        HIGH Lots of output
        LOW Little output
        MEDIUM Quite some output
        SILENT Almost no output
        Default is 'LOW'
    set_pos_file: bool
        Send updated positions to the CP2K shell via file instead of
        via stdin, which can bypass limitations for sending large
        structures via stdin for CP2K built with some MPI libraries.
        Requires CP2K 2024.2
    """

    default_parameters = dict(
        # Parameters that are not explicitly defined here. Can be defined by class arguments.
        # The input template or CP2K default values (if exist) are used if not overwritten by class arguments.
        project=None,
        # print_level=None,
        force_eval_method=None,
        stress_tensor=None,
        charge=None,
        multiplicity=None,
        uks=None,
        # xtb=None
        # xc=None,
        # basis_set_file=None,
        # basis_set=None,
        # potential_file=None,
        # pseudo_potential=None,
        poisson_solver=None,
        cutoff=None,
        # max_scf=None,

        # Parameters that are explicitly defined here. Can be overwritten by class arguments.
        # To use the input template or CP2K default values (if exist), use "None" class argument.
        # project='project',
        print_level='low',
        # force_eval_method="Quickstep",
        # stress_tensor=False,
        # charge=0,
        # multiplicity=1,
        # uks=False,
        xtb='GFN2',
        xc='pbe',
        basis_set_file='BASIS_MOLOPT',
        basis_set='DZVP-MOLOPT-SR-GTH',
        potential_file='POTENTIAL',
        pseudo_potential='auto',
        # poisson_solver='auto',
        # cutoff=400 * Rydberg,
        max_scf=1000,

        # The input template.
        inp="",

        # Additional parameters
        set_pos_file=False,
        auto_write=False,
    )

    def _generate_input(self):
        """Generates a CP2K input file"""
        p = self.parameters
        root = parse_input(p.inp)

        if p.xtb:
            p.xc=None
            p.basis_set_file=None
            p.basis_set=None
            p.potential_file=None
            p.pseudo_potential=None

        root.add_keyword('GLOBAL', 'RUN_TYPE NEGF')

        if p.project:
            root.add_keyword('GLOBAL', 'PROJECT ' + p.project)

        if p.print_level:
            root.add_keyword('GLOBAL', 'PRINT_LEVEL ' + p.print_level)

        if p.force_eval_method:
            root.add_keyword('FORCE_EVAL', 'METHOD ' + p.force_eval_method)

        if p.stress_tensor:
            root.add_keyword('FORCE_EVAL', 'STRESS_TENSOR ANALYTICAL')
            root.add_keyword('FORCE_EVAL/PRINT/STRESS_TENSOR',
                             '_SECTION_PARAMETERS_ ON')

        if p.basis_set_file:
            root.add_keyword('FORCE_EVAL/DFT',
                             'BASIS_SET_FILE_NAME ' + p.basis_set_file)

        if p.potential_file:
            root.add_keyword('FORCE_EVAL/DFT',
                             'POTENTIAL_FILE_NAME ' + p.potential_file)

        if p.cutoff:
            root.add_keyword('FORCE_EVAL/DFT/MGRID',
                             'CUTOFF [eV] %.18e' % p.cutoff)

        if p.max_scf:
            root.add_keyword('FORCE_EVAL/DFT/SCF', 'MAX_SCF %d' % p.max_scf)
            root.add_keyword('FORCE_EVAL/DFT/LS_SCF', 'MAX_SCF %d' % p.max_scf)

        if p.xtb:
            root.add_keyword('FORCE_EVAL/DFT/QS', 'METHOD xTB')
            #root.add_keyword('FORCE_EVAL/DFT/QS/XTB', 'CHECK_ATOMIC_CHARGES F')
            if p.xtb == 'GFN2':
                root.add_keyword('FORCE_EVAL/DFT/QS/XTB', 'GFN_TYPE TBLITE')
                root.add_keyword('FORCE_EVAL/DFT/QS/XTB/TBLITE', 'METHOD GFN2')

        if p.xc:
            legacy_libxc = ""
            for functional in p.xc.split():
                functional = functional.replace("LDA", "PADE")  # resolve alias
                xc_sec = root.get_subsection('FORCE_EVAL/DFT/XC/XC_FUNCTIONAL')
                # libxc input section changed over time
                if functional.startswith("XC_") and self._shell.version < 3.0:
                    legacy_libxc += " " + functional  # handled later
                elif functional.startswith("XC_") and self._shell.version < 5.0:
                    s = InputSection(name='LIBXC')
                    s.keywords.append('FUNCTIONAL ' + functional)
                    xc_sec.subsections.append(s)
                elif functional.startswith("XC_"):
                    s = InputSection(name=functional[3:])
                    xc_sec.subsections.append(s)
                else:
                    s = InputSection(name=functional.upper())
                    xc_sec.subsections.append(s)
            if legacy_libxc:
                root.add_keyword('FORCE_EVAL/DFT/XC/XC_FUNCTIONAL/LIBXC',
                                 'FUNCTIONAL ' + legacy_libxc)

        if p.uks:
            root.add_keyword('FORCE_EVAL/DFT', 'UNRESTRICTED_KOHN_SHAM ON')

        if p.multiplicity:
            root.add_keyword('FORCE_EVAL/DFT',
                             'MULTIPLICITY %d' % p.multiplicity)

        if p.charge and p.charge != 0:
            root.add_keyword('FORCE_EVAL/DFT', 'CHARGE %d' % p.charge)

        # add Poisson solver if needed
        if p.poisson_solver == 'auto' and not any(self.atoms.get_pbc()):
            root.add_keyword('FORCE_EVAL/DFT/POISSON', 'PERIODIC NONE')
            root.add_keyword('FORCE_EVAL/DFT/POISSON', 'PSOLVER  MT')

        # write coords
        syms = self.atoms.get_chemical_symbols()
        atoms = self.atoms.get_positions()
        for elm, pos in zip(syms, atoms):
            line = f'{elm} {pos[0]:.18e} {pos[1]:.18e} {pos[2]:.18e}'
            root.add_keyword('FORCE_EVAL/SUBSYS/COORD', line, unique=False)

        # write cell
        pbc = ''.join([a for a, b in zip('XYZ', self.atoms.get_pbc()) if b])
        if len(pbc) == 0:
            pbc = 'NONE'
        root.add_keyword('FORCE_EVAL/SUBSYS/CELL', 'PERIODIC ' + pbc)
        c = self.atoms.get_cell()
        for i, a in enumerate('ABC'):
            line = f'{a} {c[i, 0]:.18e} {c[i, 1]:.18e} {c[i, 2]:.18e}'
            root.add_keyword('FORCE_EVAL/SUBSYS/CELL', line)

        # determine pseudo-potential
        potential = p.pseudo_potential
        if p.pseudo_potential == 'auto':
            if p.xc and p.xc.upper() in ('LDA', 'PADE', 'BP', 'BLYP', 'PBE',):
                potential = 'GTH-' + p.xc.upper()
            else:
                msg = 'No matching pseudo potential found, using GTH-PBE'
                warn(msg, RuntimeWarning)
                potential = 'GTH-PBE'  # fall back

        # write atomic kinds
        subsys = root.get_subsection('FORCE_EVAL/SUBSYS').subsections
        kinds = {s.params: s for s in subsys if s.name == "KIND"}
        for elem in set(self.atoms.get_chemical_symbols()):
            if elem not in kinds.keys():
                s = InputSection(name='KIND', params=elem)
                subsys.append(s)
                kinds[elem] = s
            if p.basis_set:
                kinds[elem].keywords.append('BASIS_SET ' + p.basis_set)
            if potential:
                kinds[elem].keywords.append('POTENTIAL ' + potential)

        output_lines = ['!!! Generated by ASE !!!'] + root.write()
        return '\n'.join(output_lines)
