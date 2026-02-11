import os
import tempfile

def Bulk_Forming_WInput(self):
    # Plastic Input File
    plastic_parameters = open(self.BULK_FORMING_PLASTIC_PARAMETERS,"w")
    if self.TMaterials_2.item(0,23) is None or not self.TMaterials_2.item(0,23).text().strip():
        no_plasticity_input = 'SLIP SYSTEMS FOR CUBIC CRYSTAL\n' \
                      'CUBIC             icryst\n' \
                      '   1.   1.   1.   crystal axis (cdim(i))\n' \
                      '   0              nmodesx (total # of modes listed in the file)\n' \
                      '   0              nmodes (# of modes to be used in the calculation)'
        plastic_parameters.write(no_plasticity_input)
        modes = '2 0 0 0               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
    else:
        header = 'SLIP SYSTEMS FOR CUBIC CRYSTAL\n' \
                 'CUBIC             icryst\n'
        plastic_parameters.write(header)
        crystal_axis = f'   {self.TMaterials_2.item(0,23).text()}   {self.TMaterials_2.item(0,24).text()}   {self.TMaterials_2.item(0,25).text()}   crystal axis (cdim(i))\n'
        plastic_parameters.write(crystal_axis)
        system_names = self.TMaterials_2.item(0,26).text()
        num_systems = system_names.count(',') + 1
        nmodes = f'   {num_systems}              nmodesx (total # of modes listed in the file)\n' \
                 f'   {num_systems}              nmodes (# of modes to be used in the calculation)\n'
        plastic_parameters.write(nmodes)
        modei = '   1'
        for i in range(num_systems-1):
            modei = modei + f'  {i+2}'
        modei = modei + '              mode(i) (label of the modes to be used)\n'
        plastic_parameters.write(modei)
        slip_system_names = system_names.split(',')
        nsmx_max = 0
        for i in range(num_systems):
            slip_system_name = slip_system_names[i].split('.', 1)[1].strip()
            dict = slip_system_names[i].split('.', 1)[0].strip()
            dvar = f'T{dict}'
            slip_table = getattr(self, dvar)
            nsmx = slip_table.rowCount()
            if nsmx > nsmx_max:
                nsmx_max = nsmx
            slip_parameters = f'   {slip_system_name} SLIP\n' \
                              f'  {i+1}   {nsmx}   {self.TMaterials_2.item(0,27).text().split(",")[i]}   {self.TMaterials_2.item(0,28).text().split(",")[i]}   0.0   0           modex,nsmx,nrsx,gamd0x,twshx,isectwx\n' \
                              f'  {self.TMaterials_2.item(0,29).text().split(",")[i]}   {self.TMaterials_2.item(0,30).text().split(",")[i]}   {self.TMaterials_2.item(0,31).text().split(",")[i]}   {self.TMaterials_2.item(0,32).text().split(",")[i]}   {self.TMaterials_2.item(0,33).text().split(",")[i]}         tau0xf,tau0xb,tau1x,thet0,thet1\n' \
                              f'  {self.TMaterials_2.item(0,34).text().split(",")[i]}'
            for ii in range(num_systems):
                slip_parameters = slip_parameters + f'   {self.TMaterials_2.item(0,35).text().split(",")[i]}'
            slip_parameters = slip_parameters + '                            hselfx,hlatex\n'
            plastic_parameters.write(slip_parameters)
            for j in range(nsmx):
                slip_plane = slip_table.item(j,0).text().split(',')
                slip_system = '   '
                for k in range(len(slip_plane)):
                    slip_system = slip_system + f'{slip_plane[k]}   '
                slip_direction = slip_table.item(j,1).text().split(',')
                for k in range(len(slip_direction)):
                    slip_system = slip_system + f'   {slip_direction[k]}'
                if j == 0:
                    slip_system = slip_system + '        SLIP (n-b)\n'
                else:
                    slip_system = slip_system + '\n'
                plastic_parameters.write(slip_system)
        modes = f'1 1 1 {nsmx_max}               NPHMX, NMODMX, NTWMMX, NSYSMX\n'
    plastic_parameters.close()

    # Elastic Input File
    if 'Isotropic' in self.TMaterials_2.item(0,1).text() or 'Transversely Isotropic' in self.TMaterials_2.item(0,1).text() or 'Orthotropic' in self.TMaterials_2.item(0,1).text():
        elastic_parameters = open(self.BULK_FORMING_ELASTIC_PARAMETERS,"w")
        iso = '0\n'
        elastic_parameters.write(iso)
        stiffness = f'  {self.TMaterials_2.item(0,2).text()}  {self.TMaterials_2.item(0,3).text()}  {self.TMaterials_2.item(0,4).text()}  0  0  0     Cu (MPa)\n' \
                    f'  {self.TMaterials_2.item(0,3).text()}  {self.TMaterials_2.item(0,8).text()}  {self.TMaterials_2.item(0,9).text()}  0  0  0\n' \
                    f'  {self.TMaterials_2.item(0,4).text()}  {self.TMaterials_2.item(0,9).text()}  {self.TMaterials_2.item(0,13).text()}  0  0  0\n' \
                    f'  0  0  0   {self.TMaterials_2.item(0,17).text()}  0  0\n' \
                    f'  0  0  0  0  {self.TMaterials_2.item(0,20).text()}  0\n' \
                    f'  0  0  0  0  0  {self.TMaterials_2.item(0,22).text()}'
        elastic_parameters.write(stiffness)
        elastic_parameters.close()
    elif self.TMaterials_2.item(0,1).text() == 'Anisotropic':
        elastic_parameters = open(self.BULK_FORMING_ELASTIC_PARAMETERS,"w")
        iso = '0\n'
        elastic_parameters.write(iso)
        stiffness = f'  {self.TMaterials_2.item(0,2).text()}  {self.TMaterials_2.item(0,3).text()}  {self.TMaterials_2.item(0,4).text()}  {self.TMaterials_2.item(0,5).text()}  {self.TMaterials_2.item(0,6).text()}  {self.TMaterials_2.item(0,7).text()}     Cu (MPa)\n' \
                    f'  {self.TMaterials_2.item(0,3).text()}  {self.TMaterials_2.item(0,8).text()}  {self.TMaterials_2.item(0,9).text()}  {self.TMaterials_2.item(0,10).text()}  {self.TMaterials_2.item(0,11).text()}  {self.TMaterials_2.item(0,12).text()}\n' \
                    f'  {self.TMaterials_2.item(0,4).text()}  {self.TMaterials_2.item(0,9).text()}  {self.TMaterials_2.item(0,13).text()}  {self.TMaterials_2.item(0,14).text()}  {self.TMaterials_2.item(0,15).text()}  {self.TMaterials_2.item(0,16).text()}\n'  \
                    f'  {self.TMaterials_2.item(0,5).text()}  {self.TMaterials_2.item(0,10).text()}  {self.TMaterials_2.item(0,14).text()}  {self.TMaterials_2.item(0,17).text()}  {self.TMaterials_2.item(0,18).text()}  {self.TMaterials_2.item(0,19).text()}\n'  \
                    f'  {self.TMaterials_2.item(0,6).text()}  {self.TMaterials_2.item(0,11).text()}  {self.TMaterials_2.item(0,15).text()}  {self.TMaterials_2.item(0,18).text()}  {self.TMaterials_2.item(0,20).text()}  {self.TMaterials_2.item(0,21).text()}\n'  \
                    f'  {self.TMaterials_2.item(0,7).text()}  {self.TMaterials_2.item(0,12).text()}  {self.TMaterials_2.item(0,16).text()}  {self.TMaterials_2.item(0,19).text()}  {self.TMaterials_2.item(0,21).text()}  {self.TMaterials_2.item(0,22).text()}\n'
        elastic_parameters.write(stiffness)
        elastic_parameters.close()

    # Bulk Forming input parameters file
    bulk_forming_input = open(self.BULK_FORMING_INPUT,"w")
    bulk_forming_input.write(modes)
    Nx = int(self.TParts.item(0,7).text())
    Ny = int(self.TParts.item(0,8).text())
    Nz = int(self.TParts.item(0,9).text())
    dimensions = f'{Nx} {Ny} {Nz}               x-dim, y-dim, z-dim\n'
    bulk_forming_input.write(dimensions)
    nph = '1                      number of phases (nph)\n'
    bulk_forming_input.write(nph)
    dx = float(self.TParts.item(0,4).text())/Nx
    dy = float(self.TParts.item(0,5).text())/Ny
    dz = float(self.TParts.item(0,6).text())/Nz
    delt = f'{dx:.4f} {dy:.4f} {dz:.4f}             RVE dimensions (delt)\n' + '* name and path of microstructure file (filetext)\n'
    bulk_forming_input.write(delt)
    partfile = self.TParts.item(0,10).text() + '\n'
    bulk_forming_input.write(partfile)
    phase1 = '*INFORMATION ABOUT PHASE #1\n' \
             '0                          igas(iph)\n' \
             '* name and path of single crystal files (filecryspl, filecrysel) (dummy if igas(iph)=1)\n' \
             f'{self.BULK_FORMING_PLASTIC_PARAMETERS}\n' \
             f'{self.BULK_FORMING_ELASTIC_PARAMETERS}\n' \
             '*INFORMATION ABOUT TEST CONDITIONS\n'
    bulk_forming_input.write(phase1)
    boundary_conditions_title = '* boundary conditions\n'
    bulk_forming_input.write(boundary_conditions_title)
    iudot11 = 0 if self.TVgrad.item(0,0) is None or self.TVgrad.item(0,0).text() == "" else 1
    iudot12 = 0 if self.TVgrad.item(0,1) is None or self.TVgrad.item(0,1).text() == "" else 1
    iudot13 = 0 if self.TVgrad.item(0,2) is None or self.TVgrad.item(0,2).text() == "" else 1
    iudot21 = 0 if self.TVgrad.item(1,0) is None or self.TVgrad.item(1,0).text() == "" else 1
    iudot22 = 0 if self.TVgrad.item(1,1) is None or self.TVgrad.item(1,1).text() == "" else 1
    iudot23 = 0 if self.TVgrad.item(1,2) is None or self.TVgrad.item(1,2).text() == "" else 1
    iudot31 = 0 if self.TVgrad.item(2,0) is None or self.TVgrad.item(2,0).text() == "" else 1
    iudot32 = 0 if self.TVgrad.item(2,1) is None or self.TVgrad.item(2,1).text() == "" else 1
    iudot33 = 0 if self.TVgrad.item(2,2) is None or self.TVgrad.item(2,2).text() == "" else 1
    iudot = f'    {iudot11}       {iudot12}       {iudot13}           iudot     |    flag for vel.grad.\n' \
            f'    {iudot21}       {iudot22}       {iudot23}                     |    (0:unknown-1:known)\n' \
            f'    {iudot31}       {iudot32}       {iudot33}                     |\n' \
            f'                                          |\n'
    bulk_forming_input.write(iudot)
    udot11 = 0. if self.TVgradi.item(0,0) is None or self.TVgradi.item(0,0).text() == "" else self.TVgradi.item(0,0).text()
    udot12 = 0. if self.TVgradi.item(0,1) is None or self.TVgradi.item(0,1).text() == "" else self.TVgradi.item(0,1).text()
    udot13 = 0. if self.TVgradi.item(0,2) is None or self.TVgradi.item(0,2).text() == "" else self.TVgradi.item(0,2).text()
    udot21 = 0. if self.TVgradi.item(1,0) is None or self.TVgradi.item(1,0).text() == "" else self.TVgradi.item(1,0).text()
    udot22 = 0. if self.TVgradi.item(1,1) is None or self.TVgradi.item(1,1).text() == "" else self.TVgradi.item(1,1).text()
    udot23 = 0. if self.TVgradi.item(1,2) is None or self.TVgradi.item(1,2).text() == "" else self.TVgradi.item(1,2).text()
    udot31 = 0. if self.TVgradi.item(2,0) is None or self.TVgradi.item(2,0).text() == "" else self.TVgradi.item(2,0).text()
    udot32 = 0. if self.TVgradi.item(2,1) is None or self.TVgradi.item(2,1).text() == "" else self.TVgradi.item(2,1).text()
    udot33 = 0. if self.TVgradi.item(2,2) is None or self.TVgradi.item(2,2).text() == "" else self.TVgradi.item(2,2).text()
    udot = f'    {udot11}       {udot12}       {udot13}          udot   |    vel.grad\n' \
           f'    {udot21}       {udot22}       {udot23}                 |\n' \
           f'    {udot31}       {udot32}       {udot33}                 |\n' \
           f'                                          |\n'
    bulk_forming_input.write(udot)
    iscau11 = 0 if self.TCstress.item(0,0) is None or self.TCstress.item(0,0).text() == "" else 1
    iscau12 = 0 if self.TCstress.item(0,1) is None or self.TCstress.item(0,1).text() == "" else 1
    iscau13 = 0 if self.TCstress.item(0,2) is None or self.TCstress.item(0,2).text() == "" else 1
    iscau22 = 0 if self.TCstress.item(1,1) is None or self.TCstress.item(1,1).text() == "" else 1
    iscau23 = 0 if self.TCstress.item(1,2) is None or self.TCstress.item(1,2).text() == "" else 1
    iscau33 = 0 if self.TCstress.item(2,2) is None or self.TCstress.item(2,2).text() == "" else 1
    iscau = f'    {iscau11}       {iscau12}       {iscau13}           iscau     |    flag for Cauchy\n' \
            f'            {iscau22}       {iscau23}                     |\n' \
            f'                    {iscau33}                     |\n' \
            f'                                          |\n'
    bulk_forming_input.write(iscau)
    scauchy11 = 0. if self.TCstress.item(0,0) is None or self.TCstress.item(0,0).text() == "" else self.TCstress.item(0,0).text()
    scauchy12 = 0. if self.TCstress.item(0,1) is None or self.TCstress.item(0,1).text() == "" else self.TCstress.item(0,1).text()
    scauchy13 = 0. if self.TCstress.item(0,2) is None or self.TCstress.item(0,2).text() == "" else self.TCstress.item(0,2).text()
    scauchy22 = 0. if self.TCstress.item(1,1) is None or self.TCstress.item(1,1).text() == "" else self.TCstress.item(1,1).text()
    scauchy23 = 0. if self.TCstress.item(1,2) is None or self.TCstress.item(1,2).text() == "" else self.TCstress.item(1,2).text()
    scauchy33 = 0. if self.TCstress.item(2,2) is None or self.TCstress.item(2,2).text() == "" else self.TCstress.item(2,2).text()
    scauchy = f'    {scauchy11}       {scauchy12}       {scauchy13}          scauchy     |    Cauchy stress\n' \
              f'              {scauchy22}       {scauchy23}                     |\n' \
              f'                        {scauchy33}                     @\n'
    bulk_forming_input.write(scauchy)
    other = '* other\n' + \
            f'{self.INBFdt.text()}         eqincr (if ictrl>=0) or tdot (if ictrl=-1)\n' + \
            '-1              ictrl (1-6: strain comp, 0: VM eq, -1: tdot)\n'
    bulk_forming_input.write(other)
    run_conditions = f'*INFORMATION ABOUT RUN CONDITIONS\n' + \
                     f'{self.INBFloadsteps.text()}             nsteps\n' + \
                     f'{self.INBFerrortol.text()}         err\n' + \
                     f'{self.INBFmaxiter.text()}              itmax\n' + \
                     f'0               IRECOVER read grain states from STRESS.IN  (1) or not (0)?\n' + \
                     f'0               ISAVE write grain states in STRESS.OUT (1) or not (0)?\n' + \
                     f'1               IUPDATE update tex & RVE dim (1) or not (0)?\n' + \
                     f'1               IUPHARD\n' + \
                     f'1               IWTEX\n' + \
                     f'1 {self.INBFoutputsteps.text()}           IWFIELDS,IWSTEP\n' + \
                     f'0               ITHERMO (if ithermo=1, next line is filethermo)\n' + \
                     f'dummy\n'
    bulk_forming_input.write(run_conditions)
    bulk_forming_input.close()

