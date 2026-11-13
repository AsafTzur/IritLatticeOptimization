#ifndef INP_FEM_PARAMS_H
#define INP_FEM_PARAMS_H

#include <stddef.h>
#include <stdint.h>

#include "inc_irit/triv_lib.h"

#define INP_FEM_NAME_LENGTH 64
#define INP_FEM_UNIT_LENGTH 32

typedef enum INP_FEM_SIMULATION_TYPE {
    INP_FEM_SIMULATION_STRUCTURAL = 1,
    INP_FEM_SIMULATION_STEADY_STATE_THERMAL = 2
} INP_FEM_SIMULATION_TYPE;

typedef enum INP_FEM_UNIT_SYSTEM {
    INP_FEM_UNITS_SI = 1
} INP_FEM_UNIT_SYSTEM;

typedef enum INP_FEM_REGION_TYPE {
    /* Export the nodes lying on the selected trivariate boundaries. */
    INP_FEM_REGION_NODE_SET = 1,
    /* Export the element faces lying on the selected boundaries. */
    INP_FEM_REGION_SURFACE_SET = 2,
    /* Export the solid elements of the selected object; Boundaries is ignored. */
    INP_FEM_REGION_VOLUME_SET = 3
} INP_FEM_REGION_TYPE;

/*
 * A named pre-mesh selection stored in the ITD file.  The INP exporter turns
 * it into a node set, element-face surface, or element set after meshing.
 *
 * ObjectName is empty when the region belongs to the object carrying the
 * INP_FEM_PARAMS attribute.  It names a child object when the attribute is
 * attached to a containing list.  Boundaries uses IRIT's existing values,
 * for example:
 *
 *   IRIT_TRIV_U_MIN_BNDRY
 *   IRIT_TRIV_W_MAX_BNDRY
 *   IRIT_TRIV_U_MIN_BNDRY | IRIT_TRIV_W_MAX_BNDRY
 */
typedef struct INP_FEM_REGION {
    char Name[INP_FEM_NAME_LENGTH];
    char ObjectName[INP_FEM_NAME_LENGTH];
    INP_FEM_REGION_TYPE Type;
    IritTrivTVBndryType Boundaries;
} INP_FEM_REGION;

typedef struct INP_FEM_STRUCTURAL_MATERIAL {
    double YoungsModulusPa;
    double PoissonRatio;
    double DensityKgPerM3;
} INP_FEM_STRUCTURAL_MATERIAL;

enum INP_FEM_DOF_MASK {
    INP_FEM_DOF_X  = 1u << 0,
    INP_FEM_DOF_Y  = 1u << 1,
    INP_FEM_DOF_Z  = 1u << 2,
    INP_FEM_DOF_RX = 1u << 3,
    INP_FEM_DOF_RY = 1u << 4,
    INP_FEM_DOF_RZ = 1u << 5
};

typedef struct INP_FEM_SUPPORT {
    /* Must reference an INP_FEM_REGION_NODE_SET region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    uint32_t ConstrainedDofs;
    /* Prescribed SI displacement/rotation for X, Y, Z, RX, RY, RZ. */
    double Values[6];
} INP_FEM_SUPPORT;

typedef enum INP_FEM_STRUCTURAL_LOAD_TYPE {
    INP_FEM_LOAD_NODAL_FORCE = 1,
    INP_FEM_LOAD_SURFACE_TRACTION = 2,
    INP_FEM_LOAD_PRESSURE = 3,
    INP_FEM_LOAD_GRAVITY = 4
} INP_FEM_STRUCTURAL_LOAD_TYPE;

typedef struct INP_FEM_STRUCTURAL_LOAD {
    /* Nodal force uses a node region, traction/pressure use a surface region,
     * and gravity uses a volume region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    INP_FEM_STRUCTURAL_LOAD_TYPE Type;
    double MagnitudeSI;
    /* Unit vector in the global coordinate system; unused for pressure. */
    double Direction[3];
} INP_FEM_STRUCTURAL_LOAD;

typedef struct INP_FEM_STRUCTURAL_PARAMS {
    INP_FEM_STRUCTURAL_MATERIAL Material;
    size_t SupportCount;
    const INP_FEM_SUPPORT *Supports;
    size_t LoadCount;
    const INP_FEM_STRUCTURAL_LOAD *Loads;
} INP_FEM_STRUCTURAL_PARAMS;

typedef struct INP_FEM_THERMAL_MATERIAL {
    double ConductivityWPerMK;
} INP_FEM_THERMAL_MATERIAL;

typedef struct INP_FEM_FIXED_TEMPERATURE {
    /* Must reference an INP_FEM_REGION_NODE_SET region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    double TemperatureK;
} INP_FEM_FIXED_TEMPERATURE;

typedef enum INP_FEM_HEAT_LOAD_TYPE {
    INP_FEM_HEAT_FLUX = 1,
    INP_FEM_HEAT_GENERATION = 2
} INP_FEM_HEAT_LOAD_TYPE;

typedef struct INP_FEM_HEAT_LOAD {
    /* Flux uses a surface region; generation uses a volume region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    INP_FEM_HEAT_LOAD_TYPE Type;
    /* W/m^2 for flux and W/m^3 for volumetric generation. */
    double ValueSI;
} INP_FEM_HEAT_LOAD;

typedef struct INP_FEM_CONVECTION {
    /* Must reference an INP_FEM_REGION_SURFACE_SET region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    double CoefficientWPerM2K;
    double AmbientTemperatureK;
} INP_FEM_CONVECTION;

typedef struct INP_FEM_RADIATION {
    /* Must reference an INP_FEM_REGION_SURFACE_SET region. */
    char RegionName[INP_FEM_NAME_LENGTH];
    double Emissivity;
    double AmbientTemperatureK;
} INP_FEM_RADIATION;

typedef struct INP_FEM_THERMAL_PARAMS {
    INP_FEM_THERMAL_MATERIAL Material;
    size_t FixedTemperatureCount;
    const INP_FEM_FIXED_TEMPERATURE *FixedTemperatures;
    size_t HeatLoadCount;
    const INP_FEM_HEAT_LOAD *HeatLoads;
    size_t ConvectionCount;
    const INP_FEM_CONVECTION *Convections;
    size_t RadiationCount;
    const INP_FEM_RADIATION *Radiations;
} INP_FEM_THERMAL_PARAMS;

typedef enum INP_FEM_RESULT_TYPE {
    INP_FEM_RESULT_DISPLACEMENT = 1,
    INP_FEM_RESULT_STRESS = 2,
    INP_FEM_RESULT_STRAIN = 3,
    INP_FEM_RESULT_REACTION_FORCE = 4,
    INP_FEM_RESULT_STRAIN_ENERGY = 5,
    INP_FEM_RESULT_TEMPERATURE = 6,
    INP_FEM_RESULT_HEAT_FLUX = 7,
    INP_FEM_RESULT_THERMAL_GRADIENT = 8
} INP_FEM_RESULT_TYPE;

typedef enum INP_FEM_RESULT_COMPONENT {
    INP_FEM_COMPONENT_ALL = 0,
    INP_FEM_COMPONENT_X = 1,
    INP_FEM_COMPONENT_Y = 2,
    INP_FEM_COMPONENT_Z = 3,
    INP_FEM_COMPONENT_XY = 4,
    INP_FEM_COMPONENT_YZ = 5,
    INP_FEM_COMPONENT_XZ = 6,
    INP_FEM_COMPONENT_MAGNITUDE = 7,
    INP_FEM_COMPONENT_EQUIVALENT = 8,
    INP_FEM_COMPONENT_PRINCIPAL_MAX = 9,
    INP_FEM_COMPONENT_PRINCIPAL_MIN = 10
} INP_FEM_RESULT_COMPONENT;

typedef enum INP_FEM_RESULT_REDUCTION {
    INP_FEM_REDUCTION_NONE = 0,
    INP_FEM_REDUCTION_MINIMUM = 1,
    INP_FEM_REDUCTION_MAXIMUM = 2,
    INP_FEM_REDUCTION_AVERAGE = 3,
    INP_FEM_REDUCTION_SUM = 4
} INP_FEM_RESULT_REDUCTION;

typedef struct INP_FEM_RESULT_REQUEST {
    INP_FEM_RESULT_TYPE Type;
    /* Empty means the complete model. */
    char RegionName[INP_FEM_NAME_LENGTH];
    INP_FEM_RESULT_COMPONENT Component;
    INP_FEM_RESULT_REDUCTION Reduction;
    /* Explicit text keeps compound units extensible, e.g. "Pa" or "W/m^2". */
    char OutputUnit[INP_FEM_UNIT_LENGTH];
} INP_FEM_RESULT_REQUEST;

typedef struct INP_FEM_PARAMS {
    uint32_t SchemaVersion;
    INP_FEM_SIMULATION_TYPE SimulationType;
    INP_FEM_UNIT_SYSTEM Units;

    size_t RegionCount;
    const INP_FEM_REGION *Regions;

    union {
        INP_FEM_STRUCTURAL_PARAMS Structural;
        INP_FEM_THERMAL_PARAMS Thermal;
    } Analysis;

    size_t ResultRequestCount;
    const INP_FEM_RESULT_REQUEST *ResultRequests;
} INP_FEM_PARAMS;

#endif /* INP_FEM_PARAMS_H */
