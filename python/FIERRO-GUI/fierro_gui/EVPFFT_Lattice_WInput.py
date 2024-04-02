import os
import tempfile

# ===============================================
# ======= EVPFFT/LATTICE WRITE INPUT FILE =======
# ===============================================

def EVPFFT_Lattice_WInput(self, BC_index):

    # Plastic Input File
    plastic_parameters = open(self.PLASTIC_PARAMETERS,"w")
    no_plasticity_input = 'SLIP SYSTEMS FOR CUBIC CRYSTAL\n' \
                          'CUBIC             icryst\n' \
                          '   1.   1.   1.   crystal axis (cdim(i))\n' \
                          '   0              nmodesx (total # of modes listed in the file)\n' \
                          '   0              nmodes (# of modes to be used in the calculation)'
    plastic_parameters.write(no_plasticity_input)
    plastic_parameters.close()

    # Elastic Input File
    for i in range(self.TMaterials.rowCount()):
        if self.TMaterials.item(i,2).text() == 'Isotropic' or 'Transversely Isotropic' in self.TMaterials.item(i,2).text() or self.TMaterials.item(i,2).text() == 'Orthotropic':
            if i == 0:
                print("generating EP0")
                elastic_parameters = open(self.ELASTIC_PARAMETERS_0,"w")
            else:
                print("generating EP1")
                elastic_parameters = open(self.ELASTIC_PARAMETERS_1,"w")
            iso = '0\n'
            elastic_parameters.write(iso)
            stiffness = f'  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,5).text()}  0  0  0     Cu (MPa)\n' \
                        f'  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,9).text()}  {self.TMaterials.item(i,10).text()}  0  0  0\n' \
                        f'  {self.TMaterials.item(i,5).text()}  {self.TMaterials.item(i,10).text()}  {self.TMaterials.item(i,14).text()}  0  0  0\n' \
                        f'  0  0  0   {self.TMaterials.item(i,18).text()}  0  0\n' \
                        f'  0  0  0  0  {self.TMaterials.item(i,21).text()}  0\n' \
                        f'  0  0  0  0  0  {self.TMaterials.item(i,23).text()}'
            elastic_parameters.write(stiffness)
            elastic_parameters.close()
        elif self.TMaterials.item(i,2).text() == 'Anisotropic':
            if i == 0:
                elastic_parameters = open(self.ELASTIC_PARAMETERS_0,"w")
            else:
                elastic_parameters = open(self.ELASTIC_PARAMETERS_1,"w")
            iso = '0\n'
            elastic_parameters.write(iso)
            stiffness = f'  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,5).text()}  {self.TMaterials.item(i,6).text()}  {self.TMaterials.item(i,7).text()}  {self.TMaterials.item(i,8).text()}     Cu (MPa)\n' \
                        f'  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,9).text()}  {self.TMaterials.item(i,10).text()}  {self.TMaterials.item(i,11).text()}  {self.TMaterials.item(i,12).text()}  {self.TMaterials.item(i,13).text()}\n' \
                        f'  {self.TMaterials.item(i,5).text()}  {self.TMaterials.item(i,10).text()}  {self.TMaterials.item(i,14).text()}  {self.TMaterials.item(i,15).text()}  {self.TMaterials.item(i,16).text()}  {self.TMaterials.item(i,17).text()}\n'  \
                        f'  {self.TMaterials.item(i,6).text()}  {self.TMaterials.item(i,11).text()}  {self.TMaterials.item(i,15).text()}  {self.TMaterials.item(i,18).text()}  {self.TMaterials.item(i,19).text()}  {self.TMaterials.item(i,20).text()}\n'  \
                        f'  {self.TMaterials.item(i,7).text()}  {self.TMaterials.item(i,12).text()}  {self.TMaterials.item(i,16).text()}  {self.TMaterials.item(i,19).text()}  {self.TMaterials.item(i,21).text()}  {self.TMaterials.item(i,22).text()}\n'  \
                        f'  {self.TMaterials.item(i,8).text()}  {self.TMaterials.item(i,13).text()}  {self.TMaterials.item(i,17).text()}  {self.TMaterials.item(i,20).text()}  {self.TMaterials.item(i,22).text()}  {self.TMaterials.item(i,23).text()}\n'
            elastic_parameters.write(stiffness)
            elastic_parameters.close()
    
    # EVPFFT input parameters file
    evpfft_lattice_input = open(self.EVPFFT_INPUT,"w")
    modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
    evpfft_lattice_input.write(modes)
    dimensions = f'{self.TParts.item(0,7).text()} {self.TParts.item(0,8).text()} {self.TParts.item(0,9).text()}               x-dim, y-dim, z-dim\n'
    evpfft_lattice_input.write(dimensions)
    dx = float(self.TParts.item(0,4).text())/float(self.TParts.item(0,7).text())
    dy = float(self.TParts.item(0,5).text())/float(self.TParts.item(0,8).text())
    dz = float(self.TParts.item(0,6).text())/float(self.TParts.item(0,9).text())
    nph_delt = '2                      number of phases (nph)\n' + f'{dx:.4f} {dy:.4f} {dz:.4f}             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
    evpfft_lattice_input.write(nph_delt)
    vtkfile = f'{self.VTK_OUTPUT}\n'
    evpfft_lattice_input.write(vtkfile)
    for i in range(2):
        if not self.TMaterials.item(i,2) or self.TMaterials.item(i,2).text() == 'Ideal Gas':
            if not self.TMaterials.item(i,2) and i == 1 or self.TMaterials.item(i,1).text() == 'Void':
                phase1 = '*INFORMATION ABOUT PHASE #1\n' + \
                         '1                          igas(iph)\n' + \
                         '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + \
                         'dummy\n' + \
                         'dummy\n'
            else:
                phase2 = '*INFORMATION ABOUT PHASE #2\n' + \
                         '1                          igas(iph)\n' + \
                         '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + \
                         'dummy\n' + \
                         'dummy\n'
        else:
            if i == 0:
                efile = f'{self.ELASTIC_PARAMETERS_0}'
            else:
                efile = f'{self.ELASTIC_PARAMETERS_1}'
                
            if self.TMaterials.item(i,1).text() == 'Void':
                phase1 = '*INFORMATION ABOUT PHASE #1\n' + \
                         '0                          igas(iph)\n' + \
                         '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' +  \
                         f'{self.PLASTIC_PARAMETERS}\n' + \
                         efile + '\n'
            else:
                phase2 = '*INFORMATION ABOUT PHASE #2\n' + \
                         '0                          igas(iph)\n' + \
                         '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' +  \
                         f'{self.PLASTIC_PARAMETERS}\n' + \
                         efile + '\n'
    evpfft_lattice_input.write(phase1)
    evpfft_lattice_input.write(phase2)
    if self.TBCs.item(BC_index,1).text() == "x-direction":
        if "Tension" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       0       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       0                     |\n' + \
                              '                                          |\n' + \
                              '    1.0     0.      0.          udot      |    vel.grad\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    0       0       0           iscau     |    flag for Cauchy\n' + \
                              '            1       0                     |\n' + \
                              '                    1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        elif "Compression" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       0       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       0                     |\n' + \
                              '                                          |\n' + \
                              '    -1.0    0.      0.          udot      |    vel.grad\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    0       0       0           iscau     |    flag for Cauchy\n' + \
                              '            1       0                     |\n' + \
                              '                    1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        else:
            print("INVALID BOUNDARY CONDITION")
    elif self.TBCs.item(BC_index,1).text() == "y-direction":
        if "Tension" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    0       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       1       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          udot      |    vel.grad\n' + \
                              '    0.      1.0     0.                    |\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    1       0       0           iscau     |    flag for Cauchy\n' + \
                              '            0       0                     |\n' + \
                              '                    1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        elif "Compression" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    0       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       1       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          udot      |    vel.grad\n' + \
                              '    0.      -1.0    0.                    |\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    1       0       0           iscau     |    flag for Cauchy\n' + \
                              '            0       0                     |\n' + \
                              '                    1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        else:
            print("INVALID BOUNDARY CONDITION")
    elif self.TBCs.item(BC_index,1).text() == "z-direction":
        if "Tension" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    0       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       0       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          udot      |    vel.grad\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '    0.      0.      1.0                   |\n' + \
                              '                                          |\n' + \
                              '    1       0       0           iscau     |    flag for Cauchy\n' + \
                              '            1       0                     |\n' + \
                              '                    0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        elif "Compression" in self.TBCs.item(BC_index,0).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    0       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       0       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          udot      |    vel.grad\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '    0.      0.      -1.0                  |\n' + \
                              '                                          |\n' + \
                              '    1       0       0           iscau     |    flag for Cauchy\n' + \
                              '            1       0                     |\n' + \
                              '                    0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        else:
            print("INVALID BOUNDARY CONDITION")
    elif "Shear" in self.TBCs.item(BC_index,0).text():
        if "xy-direction" in self.TBCs.item(BC_index,1).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       1       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      1.0     0.          udot      |    vel.grad\n' + \
                              '    1.0     0.      0.                    |\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    0       0       0           iscau     |    flag for Cauchy\n' + \
                              '            0       0                     |\n' + \
                              '                    0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        elif "xz-direction" in self.TBCs.item(BC_index,1).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       1       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      1.0          udot     |    vel.grad\n' + \
                              '    0.      0.      0.                    |\n' + \
                              '    1.0     0.      0.                    |\n' + \
                              '                                          |\n' + \
                              '    0       0       0           iscau     |    flag for Cauchy\n' + \
                              '            0       0                     |\n' + \
                              '                    0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
        elif "yz-direction" in self.TBCs.item(BC_index,1).text():
            test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                              '* boundary conditions\n' + \
                              '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                              '    1       1       1                     |    (0:unknown-1:known)\n' + \
                              '    1       1       1                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          udot      |    vel.grad\n' + \
                              '    0.      0.      1.0                   |\n' + \
                              '    0.      1.0     0.                    |\n' + \
                              '                                          |\n' + \
                              '    0       0       0           iscau     |    flag for Cauchy\n' + \
                              '            0       0                     |\n' + \
                              '                    0                     |\n' + \
                              '                                          |\n' + \
                              '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                              '            0.      0.                    |\n' + \
                              '                    0.                    @\n'
    else:
        print("INVALID BOUNDARY CONDITION")
    evpfft_lattice_input.write(test_conditions)
    other = '* other\n' + \
            '0.0001         eqincr (if ictrl>=0) or tdot (if ictrl=-1)\n' + \
            '-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)\n'
    evpfft_lattice_input.write(other)
    run_conditions = '*INFORMATION ABOUT RUN CONDITIONS\n' + \
                     self.INNumberOfSteps.text() + '              nsteps\n' + \
                     '0.00001         err\n' + \
                     '50              itmax\n' + \
                     '0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?\n' + \
                     '0               ISAVE write grain states in STRESS.OUT (1) or not (0)?\n' + \
                     '1               IUPDATE update tex & RVE dim (1) or not (0)?\n' + \
                     '0               IUPHARD\n' + \
                     '1               IWTEX\n' + \
                     '1 10            IWFIELDS,IWSTEP\n' + \
                     '0               ITHERMO (if ithermo=1, next line is filethermo)\n' + \
                     'dummy\n'
    evpfft_lattice_input.write(run_conditions)
    evpfft_lattice_input.close()
