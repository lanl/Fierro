# Examples

## Calculating an average

The folder titled "average" contains an example showing how to calculate an average on the mesh. To run the example, follow the instructions on the main page to build the code, and then go to

```
/build/examples/average/test/
```
where "build" is the name of the folder that the code was built in.  To run the code, the user will type
```
./Average my_mesh.geo
```
where "my_mesh.geo" is the supplied mesh.  The user must supply a mesh when executing the code, a range of meshes are provided in the `meshes/` folder in the repository, see
```
examples/average/meshes
```
It is important to realize that the boundary planes specified in `input.cpp`
need to correspond to the physical extents of the mesh being run. Since each mesh has
it's own extent, `input.cpp` must be modified accordingly. For example, using
`mesh_12x1x1.geo`, the values of planes 3,4,5 need to be changed to `1.2, .12,` and
`.12`, respectively.

The output of the "Average" program, in EnSight format, can be found in the "ensight" directory. To visualize open "enisght/Testing.case" with Paraview.


### Index naming conventions
The global index spaces for the mesh (but local to a rank) are denoted with a _gid_.  The index spaces for the local mesh entities, relative to a _gid_, are denoted with a _lid_.  The index spaces in a reference element, which comes from the elements library and are not in **SWAGE**, are denoted with a _rid_.  A local refernce index, relative to a _rid_, is denoted with a _rlid_.


### Connectivity functions naming conventions
The general form of all connectivity structures, in terms of an arbitrary mesh entity, is,
```
// entity is generic for illustrative purposes and can be cell, element, node, etc., 
// likewise, something is generic for illustrative purposes and is some mesh entity 
for (int something_lid = 0; something_lid < mesh.num_something_in_entity(entity_gid); something_lid++){
int something_gid = mesh.something_in_entity(entity_gid, something_lid); 
// ...
}
```
See the **SWAGE** documentation for additional details and many examples.


## figures

The folder title "figures" contains diagrams and plots. 



