import os
import tempfile

# ===============================================
# ======= EVPFFT/LATTICE WRITE INPUT FILE =======
# ===============================================

def Homogenization_WInput(self, BC_index):

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
    materials_used = []
    for j in range(self.TMaterialAssignment.rowCount()):
        for i in range(self.TMaterials.rowCount()):
            if self.TMaterialAssignment.item(j,1).text() == self.TMaterials.item(i,0).text():
                materials_used.append(i)
                if 'Isotropic' in self.TMaterials.item(i,1).text() or 'Transversely Isotropic' in self.TMaterials.item(i,1).text() or 'Orthotropic' in self.TMaterials.item(i,1).text():
                    if j == 0:
                        elastic_parameters = open(self.ELASTIC_PARAMETERS_0,"w")
                    else:
                        elastic_parameters = open(self.ELASTIC_PARAMETERS_1,"w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = f'  {self.TMaterials.item(i,2).text()}  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,4).text()}  0  0  0     Cu (MPa)\n' \
                                f'  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,8).text()}  {self.TMaterials.item(i,9).text()}  0  0  0\n' \
                                f'  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,9).text()}  {self.TMaterials.item(i,13).text()}  0  0  0\n' \
                                f'  0  0  0   {self.TMaterials.item(i,17).text()}  0  0\n' \
                                f'  0  0  0  0  {self.TMaterials.item(i,20).text()}  0\n' \
                                f'  0  0  0  0  0  {self.TMaterials.item(i,22).text()}'
                    elastic_parameters.write(stiffness)
                    elastic_parameters.close()
                elif self.TMaterials.item(i,1).text() == 'Anisotropic':
                    if i == 0:
                        elastic_parameters = open(self.ELASTIC_PARAMETERS_0,"w")
                    else:
                        elastic_parameters = open(self.ELASTIC_PARAMETERS_1,"w")
                    iso = '0\n'
                    elastic_parameters.write(iso)
                    stiffness = f'  {self.TMaterials.item(i,2).text()}  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,5).text()}  {self.TMaterials.item(i,6).text()}  {self.TMaterials.item(i,7).text()}     Cu (MPa)\n' \
                                f'  {self.TMaterials.item(i,3).text()}  {self.TMaterials.item(i,8).text()}  {self.TMaterials.item(i,9).text()}  {self.TMaterials.item(i,10).text()}  {self.TMaterials.item(i,11).text()}  {self.TMaterials.item(i,12).text()}\n' \
                                f'  {self.TMaterials.item(i,4).text()}  {self.TMaterials.item(i,9).text()}  {self.TMaterials.item(i,13).text()}  {self.TMaterials.item(i,14).text()}  {self.TMaterials.item(i,15).text()}  {self.TMaterials.item(i,16).text()}\n'  \
                                f'  {self.TMaterials.item(i,5).text()}  {self.TMaterials.item(i,10).text()}  {self.TMaterials.item(i,14).text()}  {self.TMaterials.item(i,17).text()}  {self.TMaterials.item(i,18).text()}  {self.TMaterials.item(i,19).text()}\n'  \
                                f'  {self.TMaterials.item(i,6).text()}  {self.TMaterials.item(i,11).text()}  {self.TMaterials.item(i,15).text()}  {self.TMaterials.item(i,18).text()}  {self.TMaterials.item(i,20).text()}  {self.TMaterials.item(i,21).text()}\n'  \
                                f'  {self.TMaterials.item(i,7).text()}  {self.TMaterials.item(i,12).text()}  {self.TMaterials.item(i,16).text()}  {self.TMaterials.item(i,19).text()}  {self.TMaterials.item(i,21).text()}  {self.TMaterials.item(i,22).text()}\n'
                    elastic_parameters.write(stiffness)
                    elastic_parameters.close()
    # Select number of material phases based on input file type (.vtk = 2, .txt = 1)
    vtkfile = self.TParts.item(0,10).text() + '\n'
    file_type = os.path.splitext(vtkfile)[1]
    
    # EVPFFT input parameters file
    evpfft_lattice_input = open(self.EVPFFT_INPUT,"w")
    if ".vtk" in file_type:
        modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
        phases = '2                      number of phases (nph)\n'
        nphases = 2
    elif ".txt" in file_type:
        modes = '1 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
        phases = '1                      number of phases (nph)\n'
        nphases = 1
    else:
        print("WARNING: bad file type for homogenization solver")

    evpfft_lattice_input.write(modes)
    Nx = int(self.TParts.item(0,7).text())
    Ny = int(self.TParts.item(0,8).text())
    Nz = int(self.TParts.item(0,9).text())
    dimensions = f'{Nx} {Ny} {Nz}               x-dim, y-dim, z-dim\n'
    evpfft_lattice_input.write(dimensions)
    dx = float(self.TParts.item(0,4).text())/Nx
    dy = float(self.TParts.item(0,5).text())/Ny
    dz = float(self.TParts.item(0,6).text())/Nz
    nph_delt = phases + f'{dx:.4f} {dy:.4f} {dz:.4f}             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
    evpfft_lattice_input.write(nph_delt)
    evpfft_lattice_input.write(vtkfile)
    for i in range(nphases):
        if i < len(materials_used):
            j = materials_used[i]
            if i == 0:
                efile = f'{self.ELASTIC_PARAMETERS_0}'
                phase1 = '*INFORMATION ABOUT PHASE #1\n' + \
                         '0                          igas(iph)\n' + \
                         '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' +  \
                         f'{self.PLASTIC_PARAMETERS}\n' + \
                         efile + '\n'
            else:
                efile = f'{self.ELASTIC_PARAMETERS_1}'
            if self.TMaterialAssignment.item(i,0).text() == 'global':
                if self.TMaterials.item(j,1).text() == 'Ideal Gas':
                    phase1 = '*INFORMATION ABOUT PHASE #1\n' + \
                             '1                          igas(iph)\n' + \
                             '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + \
                             'dummy\n' + \
                             'dummy\n'
                else:
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
        else:
            phase1 = '*INFORMATION ABOUT PHASE #1\n' + \
                     '1                          igas(iph)\n' + \
                     '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' + \
                     'dummy\n' + \
                     'dummy\n'
    if nphases == 2:
        evpfft_lattice_input.write(phase1)
        evpfft_lattice_input.write(phase2)
    elif nphases == 1:
        evpfft_lattice_input.write(phase1)
    else:
        print("ERROR: Number of phases is greater than 2 or less than 1")
    # Tension x-direction
    if BC_index == 0:
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
    # Tension y-direction
    if BC_index == 1:
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
    # Tension z-direction
    if BC_index == 2:
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
    # Shear xy-direction
    if BC_index == 3:
        test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                          '* boundary conditions\n' + \
                          '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                          '    1       1       1                     |    (0:unknown-1:known)\n' + \
                          '    1       1       1                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.5     0.          udot      |    vel.grad\n' + \
                          '    0.5     0.      0.                    |\n' + \
                          '    0.      0.      0.                    |\n' + \
                          '                                          |\n' + \
                          '    0       0       0           iscau     |    flag for Cauchy\n' + \
                          '            0       0                     |\n' + \
                          '                    0                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                          '            0.      0.                    |\n' + \
                          '                    0.                    @\n'
    # Shear xz-direction
    if BC_index == 4:
        test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                          '* boundary conditions\n' + \
                          '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                          '    1       1       1                     |    (0:unknown-1:known)\n' + \
                          '    1       1       1                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.      0.5          udot     |    vel.grad\n' + \
                          '    0.      0.      0.                    |\n' + \
                          '    0.5     0.      0.                    |\n' + \
                          '                                          |\n' + \
                          '    0       0       0           iscau     |    flag for Cauchy\n' + \
                          '            0       0                     |\n' + \
                          '                    0                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                          '            0.      0.                    |\n' + \
                          '                    0.                    @\n'
    # Shear yz-direction
    if BC_index == 5:
        test_conditions = '*INFORMATION ABOUT TEST CONDITIONS\n' + \
                          '* boundary conditions\n' + \
                          '    1       1       1           iudot     |    flag for vel.grad.\n' + \
                          '    1       1       1                     |    (0:unknown-1:known)\n' + \
                          '    1       1       1                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.      0.          udot      |    vel.grad\n' + \
                          '    0.      0.      0.5                   |\n' + \
                          '    0.      0.5     0.                    |\n' + \
                          '                                          |\n' + \
                          '    0       0       0           iscau     |    flag for Cauchy\n' + \
                          '            0       0                     |\n' + \
                          '                    0                     |\n' + \
                          '                                          |\n' + \
                          '    0.      0.      0.          scauchy   |    Cauchy stress\n' + \
                          '            0.      0.                    |\n' + \
                          '                    0.                    @\n'
    evpfft_lattice_input.write(test_conditions)
    other = '* other\n' + \
            '0.0001         eqincr (if ictrl>=0) or tdot (if ictrl=-1)\n' + \
            '-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)\n'
    evpfft_lattice_input.write(other)
    run_conditions = f'*INFORMATION ABOUT RUN CONDITIONS\n' + \
                     f'{self.INNumberOfSteps.text()}             nsteps\n' + \
                     f'{self.INErrorTolerance.text()}         err\n' + \
                     f'{self.INMaxIterations.text()}              itmax\n' + \
                     f'0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?\n' + \
                     f'0               ISAVE write grain states in STRESS.OUT (1) or not (0)?\n' + \
                     f'1               IUPDATE update tex & RVE dim (1) or not (0)?\n' + \
                     f'0               IUPHARD\n' + \
                     f'1               IWTEX\n' + \
                     f'1 {self.INNumberOfSteps.text()}           IWFIELDS,IWSTEP\n' + \
                     f'0               ITHERMO (if ithermo=1, next line is filethermo)\n' + \
                     f'dummy\n'
    evpfft_lattice_input.write(run_conditions)
    evpfft_lattice_input.close()
