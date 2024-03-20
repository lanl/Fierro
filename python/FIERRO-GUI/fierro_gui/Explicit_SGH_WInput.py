# =============================================
# ======= EXPLICIT/SGH WRITE INPUT FILE =======
# =============================================

def Explicit_SGH_WInput(self):
#    InputFile = open("Explicit_SGH_Input.yaml","w")
    InputFile = open(self.EXPLICIT_SGH_INPUT,"w")
    
    dimensions = 'num_dims: 3\n'
    InputFile.write(dimensions)
    
    dynamic_opts = 'dynamic_options:' + '\n' \
                   '  time_final: 1.0' + '\n' \
                   '  dt_min: 1.e-8' + '\n' \
                   '  dt_max: 1.e-2' + '\n' \
                   '  dt_start: 1.e-5' + '\n' \
                   '  cycle_stop: 2000000' + '\n\n'
    InputFile.write(dynamic_opts)
   
    mesh_gen_opts = f'input_options:' + '\n' \
                    f'  mesh_file_format: vtk' + '\n' \
                    f'  mesh_file_name: {self.GLOBAL_MESH_OUTPUT}' + '\n' \
                    f'  element_type: hex8'+ '\n' \
                    f'  zero_index_base: true' + '\n\n'
    InputFile.write(mesh_gen_opts)
    
    ouput_opts = 'output_options:' + '\n' \
                 '  timer_output_level: thorough' + '\n' \
                 '  output_file_format: vtk' + '\n' \
                 '  graphics_step: 0.25' + '\n\n'
    InputFile.write(ouput_opts)
    
    materials_label = 'materials:' + '\n'
    InputFile.write(materials_label)
    for i in range(self.TMaterialsSGH.rowCount()):
        materials = f'  - id: {i}' + '\n' \
                    f'    eos_model: {self.TMaterialsSGH.item(i,1).text()}' + '\n' \
                    f'    q1: {self.TMaterialsSGH.item(i,2).text()}' + '\n' \
                    f'    q2: {self.TMaterialsSGH.item(i,3).text()}' + '\n' \
                    f'    q1ex: {self.TMaterialsSGH.item(i,4).text()}' + '\n' \
                    f'    q2ex: {self.TMaterialsSGH.item(i,5).text()}' + '\n' \
                     '    global_vars:' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,6).text()}' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,7).text()}' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,8).text()}' + '\n' \
                    f'    strength_model: {self.TMaterialsSGH.item(i,9).text()}' + '\n\n'
        InputFile.write(materials)
    
    regions_label = 'regions:' + '\n'
    InputFile.write(regions_label)
    for i in range(self.Tassignmat.rowCount()):
        if self.Tassignmat.item(i,0).text() == "global":
            type = f'  - volume:' + '\n' \
                   f'      type: global' + '\n'
        else:
            vtk_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.Tassignmat.item(i,0).text()) + '.vtk'
            type = f'  - volume:' + '\n' \
                   f'      type: vtk' + '\n' \
                   f'      stl_file_path: {vtk_path}' + '\n'
        InputFile.write(type)
        
        for j in range(self.TMaterialsSGH.rowCount()):
            if self.Tassignmat.item(j,1).text() == self.TMaterialsSGH.item(j,0).text():
                material_id = f'    material_id: {j}' + '\n'
                InputFile.write(material_id)
                break
        
        reg_info = f'    den: {self.Tassignmat.item(i,2).text()}' + '\n' \
                   f'    sie: {self.Tassignmat.item(i,3).text()}' + '\n\n' \
                   f'    velocity: cartesian' + '\n' \
                   f'    u: {self.Tassignmat.item(i,4).text()}' + '\n' \
                   f'    v: {self.Tassignmat.item(i,5).text()}' + '\n' \
                   f'    w: {self.Tassignmat.item(i,6).text()}' + '\n\n'
        InputFile.write(reg_info)
        
    fea_parameters_label = 'fea_module_parameters:' + '\n' \
                           '  - type: SGH' + '\n' \
                           '    material_id: 0' + '\n' \
                           '    boundary_conditions:' + '\n'
    InputFile.write(fea_parameters_label)
    for i in range(self.TBoundaryConditions.rowCount()):
        surface = f'      - surface:' + '\n' \
                  f'          type: {self.TBoundaryConditions.item(i,0).text()}' + '\n' \
                  f'          plane_position: {self.TBoundaryConditions.item(i,1).text()}' + '\n' \
                  f'        type: {self.TBoundaryConditions.item(i,2).text()}' + '\n\n'
        InputFile.write(surface)

    InputFile.close()
