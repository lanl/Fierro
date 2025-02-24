# =============================================
# ======= EXPLICIT/SGH WRITE INPUT FILE =======
# =============================================

def SGH_WInput(self):
#    InputFile = open("Explicit_SGH_Input.yaml","w")
    InputFile = open(self.EXPLICIT_SGH_INPUT,"w")
    
    dimensions = 'num_dims: 3\n'
    InputFile.write(dimensions)
    
    dynamic_opts = f'dynamic_options:' + '\n' \
                   f'  time_final: {self.INTime.text()}' + '\n' \
                   f'  dt_min: {self.INMindt.text()}' + '\n' \
                   f'  dt_max: {self.INMaxdt.text()}' + '\n' \
                   f'  dt_start: {self.INInitialdt.text()}' + '\n' \
                   f'  cycle_stop: {self.INmaxcycles.text()}' + '\n\n'
    InputFile.write(dynamic_opts)
   
    mesh_gen_opts = f'input_options:' + '\n' \
                    f'  mesh_file_format: vtk' + '\n' \
                    f'  mesh_file_name: {self.GLOBAL_MESH_OUTPUT}' + '\n' \
                    f'  element_type: hex8'+ '\n' \
                    f'  zero_index_base: true' + '\n\n'
    InputFile.write(mesh_gen_opts)
    
    ouput_opts = f'output_options:' + '\n' \
                 f'  timer_output_level: thorough' + '\n' \
                 f'  output_file_format: vtk' + '\n' \
                 f'  output_file_location: {self.directory}/outputs/' + '\n' \
                 f'  graphics_step: {self.INGraphicsOutput.text()}' + '\n\n'
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
                     '    eos_global_vars:' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,6).text()}' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,7).text()}' + '\n' \
                    f'      - {self.TMaterialsSGH.item(i,8).text()}' + '\n\n'
#                    f'    strength_model: {self.TMaterialsSGH.item(i,9).text()}' + '\n\n'
        InputFile.write(materials)
    
    regions_label = 'regions:' + '\n'
    InputFile.write(regions_label)
    for i in range(self.Tassignmat.rowCount()):
        if self.Tassignmat.item(i,0).text() == "global":
            type = f'  - volume:' + '\n' \
                   f'      type: global' + '\n'
        else:
            for k in range(self.TParts.rowCount()):
                if self.Tassignmat.item(i,0).text() == self.TParts.item(k,0).text():
                    vtk_path = self.voxelizer_dir + '/VTK_Geometry_' + str(self.Tassignmat.item(i,0).text()) + '.vtk'
                    type = f'  - volume:' + '\n' \
                           f'      type: vtk' + '\n' \
                           f'      vtk_file_path: {vtk_path}' + '\n'
                    break
            for kk in range(self.TBasicGeometries.rowCount()):
                if self.Tassignmat.item(i,0).text() == self.TBasicGeometries.item(kk,0).text():
                    if self.TBasicGeometries.item(kk,1).text() == 'box':
                        type = f'  - volume:' + '\n' \
                               f'      type: {self.TBasicGeometries.item(kk,1).text()}' + '\n' \
                               f'      x1: {self.TBasicGeometries.item(kk,2).text()}' + '\n' \
                               f'      x2: {self.TBasicGeometries.item(kk,3).text()}' + '\n' \
                               f'      y1: {self.TBasicGeometries.item(kk,4).text()}' + '\n' \
                               f'      y2: {self.TBasicGeometries.item(kk,5).text()}' + '\n' \
                               f'      z1: {self.TBasicGeometries.item(kk,6).text()}' + '\n' \
                               f'      z2: {self.TBasicGeometries.item(kk,7).text()}' + '\n'
                        break
                    if self.TBasicGeometries.item(kk,1).text() == 'sphere':
                        type = f'  - volume:' + '\n' \
                               f'      type: {self.TBasicGeometries.item(kk,1).text()}' + '\n' \
                               f'      radius1: {self.TBasicGeometries.item(kk,8).text()}' + '\n' \
                               f'      radius2: {self.TBasicGeometries.item(kk,9).text()}' + '\n' \
                               f'      orig_x: {self.TBasicGeometries.item(kk,10).text()}' + '\n' \
                               f'      orig_y: {self.TBasicGeometries.item(kk,11).text()}' + '\n' \
                               f'      orig_z: {self.TBasicGeometries.item(kk,12).text()}' + '\n'
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
