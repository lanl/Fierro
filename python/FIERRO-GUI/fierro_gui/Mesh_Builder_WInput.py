import fierro_mesh_builder
import paraview.simple as pvsimple
from PySide6.QtCore import (QProcess)

# =============================================
# ======= MESH BUILDER WRITE INPUT FILE =======
# =============================================

#GLOBAL_MESH = 'global_mesh.yaml'

def Mesh_Builder_WInput(self, GLOBAL_MESH, global_mesh_dir):
#    def global_mesh_click():
    global_mesh_input = open(GLOBAL_MESH,"w")
    output = f'output:\n' + \
             f'  name: mesh\n' + \
             f'  file_location: {global_mesh_dir}\n' + \
             f'  file_type: VTK\n\n'
    global_mesh_input.write(output)
    if str(self.INElementType.currentText()) == "Linear":
        order = '  p_order: 1'
    elif str(self.INElementType.currentText()) == "Quadratic":
        order = '  p_order: 2'
    else:
        order = '  p_order: 3'
    if str(self.INCoordinateSystem.currentText()) == "Rectangular":
        type = '  type: Box\n'
        if str(self.INDimension.currentText()) == "2D":
            length = '  length: [' + self.INLengthXR2D.text() + ', ' + self.INLengthYR2D.text() + ']\n'
            elements = '  num_elems: [' + self.INElementsXR2D.text() + ', ' + self.INElementsYR2D.text() + ']\n'
            origin = '  origin: [' + self.INOriginXR2D.text() + ', ' + self.INOriginYR2D.text() + ']\n'
        else:
            length = '  length: [' + self.INLengthXR3D.text() + ', ' + self.INLengthYR3D.text() + ', ' + self.INLengthZR3D.text() + ']\n'
            elements = '  num_elems: [' + self.INElementsXR3D.text() + ', ' + self.INElementsYR3D.text() + ', ' + self.INElementsZR3D.text() + ']\n'
            origin = '  origin: [' + self.INOriginXR3D.text() + ', ' + self.INOriginYR3D.text() + ', ' + self.INOriginZR3D.text() + ']\n'
        input = 'input:\n' + type + length + elements + origin + order
    else:
        type = '  type: Cylinder\n'
        if str(self.INDimension.currentText()) == "2D":
            length = '  length: [' + self.INLengthOutRadC2D.text() + ', ' + self.INLengthThetaC2D.text() + ']\n'
            elements = '  num_elems: [' + self.INElementsRadialC2D.text() + ', ' + self.INElementsArcC2D.text() + ']\n'
            origin = '  origin: [' + self.INOriginXC2D.text() + ', ' + self.INOriginYC2D.text() + ']\n'
            inner_radius = '\n  inner_radius: ' + self.INInnerRadiusC2D.text()
        else:
            length = '  length: [' + self.INLengthOutRadC3D.text() + ', ' + self.INLengthThetaC3D.text() + ', ' + self.INLengthZC3D.text() + ']\n'
            elements = '  num_elems: [' + self.INElementsRadC3D.text() + ', ' + self.INElementsArcC3D.text() + ', ' + self.INElementsZC3D.text() + ']\n'
            origin = '  origin: [' + self.INOriginXC3D.text() + ', ' + self.INOriginYC3D.text() + ', ' + self.INOriginZC3D.text() + ']\n'
            inner_radius = '\n  inner_radius: ' + self.INInnerRadiusC3D.text()
        input = 'input:\n' + type + length + elements + origin + order + inner_radius
    global_mesh_input.write(input)
    global_mesh_input.close()
    
    # Run Mesh Builder Tool
    import subprocess
    executable_path = "/Users/shankins/Documents/FY24/Github/XcodeFierro/Fierro/build-fierro-serial/bin/fierro-mesh-builder"
    arguments = [self.GLOBAL_MESH]
    command = [executable_path] + arguments
    process = subprocess.Popen(command)
    process.wait()
    self.mesh_builder = QProcess()
#    self.mesh_builder.start("fierro-mesh-builder","/Users/shankins/Documents/FY24/Github/XcodeFierro/Fierro/python/FIERRO-GUI/global_mesh.yaml")
#            self.p.start("evpfft",["-f", EVPFFT_INPUT, "-m", "2"])

#    fierro_mesh_builder.build_mesh_from_file(GLOBAL_MESH)
    
    # View Global Mesh in Paraview Window
    mesh_dir = global_mesh_dir + '/mesh.vtk'
    self.mesh = pvsimple.LegacyVTKReader(FileNames = mesh_dir)
    pvsimple.SetDisplayProperties(Representation = "Wireframe")
    pvsimple.Show(self.mesh, self.render_view)
    pvsimple.ResetCamera(view=None)
#    self.BGenerateGlobalMesh.clicked.connect(global_mesh_click)
