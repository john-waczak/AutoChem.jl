module AutoChem

using JSON
using DelimitedFiles

abstract type Reaction end
export Reaction

# Write your package code here.
include("bimolecular-reactions.jl")


export BimolecularReaction, parse_bimol_d

end
