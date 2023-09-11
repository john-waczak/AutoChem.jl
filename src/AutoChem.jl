module AutoChem

using JSON
using DelimitedFiles
using CSV, DataFrames


# assets folder
#const assets_path = @path normpath(joinpath(@__DIR__, "../assets"))
const assets_path = normpath(joinpath(@__DIR__, "../assets"))
# @assert ispath(assets_path)



abstract type Reaction end
export Reaction


include("bimolecular-reactions.jl")
export BimolecularReaction, parse_bimol_d, read_bimol

include("trimolecular-reactions.jl")
export TrimolecularReaction, parse_trimol_d, read_trimol

include("photolysis-reactions.jl")
export PhotolysisReaction, parse_photolysis_d, read_photolysis
export FittedPhotolysisReaction, read_fitted_photolysis

include("reference_measurements.jl")
export generate_densities

include("initialize.jl")
export generate_init_dict

include("stoich_mats.jl")
export generate_stoich_mat

include("derivatives.jl")
export RxnDerivative, BimolecularDerivativeTerm, TrimolecularDerivativeTerm, PhotolysisDerivativeTerm, get_derivative_terms, get_bimolecular_derivatives, get_trimolecular_derivatives, get_photolysis_derivatives






end
