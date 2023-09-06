var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AutoChem","category":"page"},{"location":"#AutoChem","page":"Home","title":"AutoChem","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AutoChem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AutoChem]","category":"page"},{"location":"#AutoChem.BimolecularReaction-Tuple{Any, Any, Any}","page":"Home","title":"AutoChem.BimolecularReaction","text":"(rxn::BimolecularReaction)(T,P,M)\n\nGiven a reaction rxn of type BimolecularReaction, compute the reaction rate coefficient as a function of\n\nT: temperature in Kelvin\nP: pressure in mbar\nM: the total particle number density\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.TrimolecularReaction-Tuple{Any, Any, Any}","page":"Home","title":"AutoChem.TrimolecularReaction","text":"(rxn::TrimolecularReaction)(T,P,M)\n\nGiven a reaction rxn of type TrimolecularReaction, compute the reaction rate coefficient as a function of\n\nT: temperature in Kelvin\nP: pressure in mbar\nM: the total particle number density\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.parse_bimol_d-Tuple{Any}","page":"Home","title":"AutoChem.parse_bimol_d","text":"parse_bimol_d(path)\n\nParse an AutoChem bimolecular reaction database, returning a vector of BimolecularReaction objects and a second vector containing the indices of any reactions that failed to parse.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.parse_photolysis_d-Tuple{Any}","page":"Home","title":"AutoChem.parse_photolysis_d","text":"parse_photolysis_d(path)\n\nParse an AutoChem photolysis reaction database, returning a vector of PhotolysisReaction objects and a second vector containing the indices of any reactions that failed to parse.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.parse_trimol_d-Tuple{Any}","page":"Home","title":"AutoChem.parse_trimol_d","text":"parse_trimol_d(path)\n\nParse an AutoChem trimolecular reaction database, returning a vector of TrimolecularReaction objects and a second vector containing the indices of any reactions that failed to parse.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.read_bimol-Tuple{Any}","page":"Home","title":"AutoChem.read_bimol","text":"read_bimol(path)\n\nParse a bimolecular reaction database (in JSON format) returning a vector of BimolecularReaction objects.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.read_fitted_photolysis-Tuple{Any}","page":"Home","title":"AutoChem.read_fitted_photolysis","text":"read_fitted_photolysis(path)\n\nParse a photolysis reaction database (in JSON format) returning a vector of FittedPhotolysisReaction objects.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.read_photolysis-Tuple{Any}","page":"Home","title":"AutoChem.read_photolysis","text":"read_photolysis(path)\n\nParse a photolysis reaction database (in JSON format) returning a vector of PhotolysisReaction objects.\n\n\n\n\n\n","category":"method"},{"location":"#AutoChem.read_trimol-Tuple{Any}","page":"Home","title":"AutoChem.read_trimol","text":"read_trimol(path)\n\nParse a trimolecular reaction database (in JSON format) returning a vector of TrimolecularReaction objects.\n\n\n\n\n\n","category":"method"}]
}
