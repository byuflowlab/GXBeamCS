module GXBeamCS

import ForwardDiff
using SparseArrays
using StaticArrays
using FLOWMath
using LinearAlgebra
using DelimitedFiles
using RecipesBase
using ImplicitAD
using ReverseDiff

# from common
export Material, MaterialPlane, Layer
export compliance_matrix, mass_matrix, strains_and_stresses, tsai_wu
export plotmesh, plotsoln

# from fem
export Node, MeshElement, FEM
export initialize_cache

# from afmesh
export afmesh

# from clt
export BeamSection, CLT, clt_stiffness_matrix
export compliance_and_operators, strains_stresses_from_B, StrainOperatorSet
export num_strain_locs, node_operator, operator_fields
export strain_loc_index, strain_loc_tuple, strain_loc_indices

# shared definitions and functions
include("common.jl")

# section properties using finite element based approach
include("fem.jl")

# airfoil meshing
include("afmesh.jl")

# section properties using classical laminate theory
include("clt.jl")

# meshing tools
include("meshtools.jl")

end