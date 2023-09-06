module AutoChem

using JSON
using DelimitedFiles



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


end
