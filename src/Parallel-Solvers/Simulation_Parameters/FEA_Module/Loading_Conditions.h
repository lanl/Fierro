#pragma once

#include "yaml-serializable.h"
#include "Simulation_Parameters/FEA_Module/Boundary_Conditions.h"
#include "Simulation_Parameters/SerializableGeometry.h"

SERIALIZABLE_ENUM(LOADING_CONDITION_TYPE, surface_traction, surface_heat_flux, body_force, surface_force)
SERIALIZABLE_ENUM(LOADING_SPECIFICATION, normal, coordinated)

struct loading_t {
    LOADING_CONDITION_TYPE condition_type = LOADING_CONDITION_TYPE::body_force;
    Surface surface;
    volume_t volume;
    double x, y, z;

    KOKKOS_FUNCTION
    bool contains(const double* elem_coords) {
      return volume.contains(elem_coords);
    }
    
    KOKKOS_FUNCTION
    size_t planar_surface_index() {
      switch (surface.type) {
        case BOUNDARY_TYPE::x_plane:
          return 0;
        case BOUNDARY_TYPE::y_plane:
          return 1;
        case BOUNDARY_TYPE::z_plane:
          return 2;
        default:
          // Not sure what to do about this error case, since we can't
          // throw on the device.
          return -1;
      }
    }
};

struct Loading_Condition 
  : virtual loading_t, 
    Yaml::TypeDiscriminated<Loading_Condition, LOADING_CONDITION_TYPE>, 
    Yaml::DerivedFields {
  void derive() {
    loading_t::condition_type = type;
  }

  virtual ~Loading_Condition() { }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Loading_Condition, type)
IMPL_YAML_SERIALIZABLE_FOR(Loading_Condition, type)

struct Surface_Loading : virtual Loading_Condition { };
YAML_ADD_REQUIRED_FIELDS_FOR(Surface_Loading, surface)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Surface_Loading, Loading_Condition, surface)

struct Surface_Traction_Condition 
  : Surface_Loading, Loading_Condition::Register<Surface_Traction_Condition, LOADING_CONDITION_TYPE::surface_traction> {
  
  double component_x;
  double component_y;
  double component_z;
  
  void derive() {
    loading_t::x = component_x;
    loading_t::y = component_y;
    loading_t::z = component_z;
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Surface_Traction_Condition, component_x, component_y, component_z)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Surface_Traction_Condition, Surface_Loading, component_x, component_y, component_z)

struct Surface_Flux_Condition 
  : Surface_Loading, Loading_Condition::Register<Surface_Flux_Condition, LOADING_CONDITION_TYPE::surface_heat_flux> {
  
  void derive() { }

  double flux_value;
  LOADING_SPECIFICATION specification;
};
YAML_ADD_REQUIRED_FIELDS_FOR(Surface_Flux_Condition, flux_value, specification)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Surface_Flux_Condition, Surface_Loading, flux_value, specification)


struct Volume_Loading : virtual Loading_Condition {
  std::shared_ptr<Volume> volume;
  void derive() {
    loading_t::volume = (volume_t)*volume;
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Volume_Loading, volume)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Volume_Loading, Loading_Condition, volume)

struct Body_Force_Condition 
  : Volume_Loading, 
    Loading_Condition::Register<Body_Force_Condition, LOADING_CONDITION_TYPE::body_force> {

  double component_x;
  double component_y;
  double component_z;
  
  void derive() {
    loading_t::x = component_x;
    loading_t::y = component_y;
    loading_t::z = component_z;
  }
};
YAML_ADD_REQUIRED_FIELDS_FOR(Body_Force_Condition, component_x, component_y, component_z)
IMPL_YAML_SERIALIZABLE_WITH_BASE(Body_Force_Condition, Volume_Loading, component_x, component_y, component_z)
