module AutoChem

using JSON
using DelimitedFiles
using RelocatableFolders



# assets folder
const autochem_data_path = @path normpath(joinpath(@__DIR__, "../assets", "autochem", "data"))
@assert ispath(autochem_data_path)



abstract type Reaction end
export Reaction

include("bimolecular-reactions.jl")
export BimolecularReaction, parse_bimol_d, read_bimol

include("trimolecular-reactions.jl")
export TrimolecularReaction, parse_trimol_d, read_trimol

include("photolysis-reactions.jl")
export PhotolysisReaction, parse_photolysis_d, read_photolysis


end
