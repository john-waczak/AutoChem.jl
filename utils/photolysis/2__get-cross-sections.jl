using DelimitedFiles
using CSV, DataFrames
# using Plots
using Unitful
using ProgressMeter

# 0. Set up paths

basepath = "../../assets/mpi-mainz-uvviz/raw"
crosssections_path = joinpath(basepath, "cross-sections")
quantumyields_path = joinpath(basepath, "quantum-yields")

@assert ispath(basepath)

if !ispath(crosssections_path)
    mkpath(crosssections_path)
end

if !ispath(quantumyields_path)
    mkpath(quantumyields_path)
end


# 1. Get list of all summary.csv files
summary_list = String[]

for (root, dirs, files) in walkdir(crosssections_path)
    for f ∈ files
        if occursin("summary", f) && endswith(f, ".csv")
            push!(summary_list, joinpath(root, f))
        end
    end
end


# generate list of urls for testing purposes:
url_list = []

@showprogress 1 for i ∈ 1:length(summary_list)
    summary_df = CSV.File(summary_list[i]) |> DataFrame
    push!(url_list, summary_df.download_url)
end

url_list = vcat(url_list...)


function download_data(url)
    tmp = download(url)
    data = readdlm(tmp)

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkynes,polyynes+radicals/Alkynes/C3H4_Benilan(1999)_293K_185-215nm.txt"
        data[888, 2] = 6.92406e-21
        data[889, 1] = 185.0000
        data[889, 2] = 3.97507e-18
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Aromatic%20compounds/Benzene,%20biphenyl/C6H6_Richardson(1982)_999K_231-316nm.txt"
        data[147, 2] = 82.6369
        data = [data; [231.12 3.29e-19 ""]]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Aromatic%20compounds/Nitro%20compounds/m-CH3-o-C6H3(OH)(NO2)_IUPAC(2012)_298K_320-450nm(rec).txt"
        # this one has a weird bit about a "corrected value"
        data[105, 2] = 1.0e-19
        data = data[1:end-1, :]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Aromatic%20compounds/Nitro%20compounds/m-CH3-o-C6H3(OH)(NO2)_Chen(2011)_293K_320-450nm.txt"
        # this one has a weird bit about a "corrected"value
        data[105, 2] = 1.0e-19
        data = data[1:end-1, :]
    elseif all(elem == "" for elem ∈ data[end,:])
        data = data[1:end-1, :]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Carbon-oxides/CO_SunWeissler(1955)_295K_37-131nm.txt"
        data = data[1:end-1, :]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Carbon-oxides/CO2_Jensen(1997)_1523K_250nm.txt"
        data = [
            "wavel." "cross sect." "error limit"
            250 2.9e-21 0.1e-21
        ]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Carbon-oxides/CO2_Jensen(1997)_1818K_250nm.txt"
        data = [
            "wavel." "cross sect." "error limit"
            250 1.1e-20 0.3e-20
        ]
    elseif data[end, 2] == "means"
        data = data[1:end-1, :]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ethers+alkyl%20radicals/Cyclic%20ethers/C5H8O_Doner(2019)_323K_125-240nm.txt"
        data[1671, 2] = 8.436e-18
        data[1671, 3] = 8.9e-20

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Br%20oxides/BrO_JPL-2010(2011)_228K_286.46-381.27nm(band%20peaks,rec).txt"
        data = data[1:end-1, 3:4]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Br%20oxides/BrO_JPL-2010(2011)_298K_286.46-384.87nm(band%20peaks,rec).txt"
        data = data[1:end-4, 3:4]

    elseif (size(data, 2) == 7) && data[1,7] == "%"
        data = data[1:end, 2:3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_GillespieDonovan(1976)_293K_277.2-285.2nm.txt"
        data = data[3:end, 2:3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_Jourdain(1978)_298K_272.9-312.5nm.txt"
        data = data[3:end, 2:3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_MandelmanNicholls(1977)_298K_248.0-307.8nm.txt"
        data = [
            data[2:8, 1:2]
            data[10:end, 2:3]
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_Nicholls(1975)_298K_303.5-323.9nm.txt"
        data = data[3:end, 2:3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_Porter(1950)_298K_263.63-303.45nm(band%20heads).txt"
        data = data[2:end, 2:3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClOOCl_Burkholder(1990)_250K_245nm(max).txt"
        data = [245 6.5e-18 0.8e-18]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/I%20oxides/IO_CoxCoker(1983)_303K_426.9nm.txt"
        data = [426.9 3.1e-17 2.0e-17]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogenated%20N-compounds(inorg)/BrNO_HouelVandenBergh(1977)_300K_215-708nm.txt"
        data = data[2:end-1, 1:2]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogenated%20N-compounds(inorg)/FNO_JPL-2010(2011)_298K_180-350nm(rec).txt"
        data[131, 2] = 1.18e-20
        [
            data
            [180 5.24e-19 ""]
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogenated%20N-compounds(inorg)/INO_IUPAC(2000)_298K_230-460nm(rec).txt"
        data = data[1:end-1, :]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Chloroalkanes/CCl4_PrahladKumar(1995)_300K_186-240nm.txt"
        data[:,2] = ([parse(Float64, data[i,2]*"$(data[i,3])") for i ∈ 1:size(data,1)])
        data = data[:, 1:2]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Fluoroalkanes/C4F10_Ravishankara(1993)_298K_121.6nm.txt"
        data = [121.6 2.0e-18]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Freons-CFC(C,F,Cl)/CCl3CClF2(CFC-112a)_Davis(2016)_270K_192.5-235nm.txt"
        data[5, 3] = data[5, 4]

        data = data[1:end-1, :]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Halons(C,H,F,Cl,Br)/CF2Br2(Halon-1202)_Doucet(1975)_298K_121.9-220.1nm.txt"
        data[104, 2] = 3.75e-18
        data = [
            data
            [121.87 9.98e-17 ""]
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Iodoalkanes(C,H,F,Cl,Br,I)/n-C3H7I_BoschiSalahub(1972)_298K_213.4-299.5nm.txt"
        data[99, 2] = 2.68e-18
        data = [
            data
            [213.40 5.19e-20 ""]
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogens+mixed%20halogens/BrCl_JPL-2010(2011)_298K_200-600nm(rec).txt"
        data = data[1:end-1, :]

    elseif data[end, 1] == "*)"
        data = data[1:end-1, :]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Hydrogen+water/H2O_Aldener(2005)_298K_870nm.txt"
        data = [870 9.2e-27 0.2e-27]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Nitrogen%20oxides/NO2_JPL-2010(2011)_220K_240-662.5nm(rec).txt"
        # this one gives upper and lower λ band at each level
        λ = Float64[]
        σ = data[:,2]

        for i ∈ 1:size(data, 1)
            splits = split(data[i,1], "-")
            push!(λ, (parse(Float64, splits[1]) + parse(Float64, splits[2]))/2)
        end
        data[:,1] = λ

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Nitrogen%20oxides/NO2_JPL-2010(2011)_294K_240-662.5nm(rec).txt"
        # this one gives upper and lower λ band at each level
        λ = Float64[]
        σ = data[:,2]

        for i ∈ 1:size(data, 1)
            splits = split(data[i,1], "-")
            push!(λ, (parse(Float64, splits[1]) + parse(Float64, splits[2]))/2)
        end
        data[:,1] = λ

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Nitrogen%20oxides/NO2_Voigt(2002)_223K_250-312nm(1000mbar).txt"
        data[2777, 2] = -2.811e-19
        data = [
            data[1:277, :]
            data[279:end, :]
        ]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Nitrogen%20oxides/NO2_Voigt(2002)_246K_250-312nm(1000mbar).txt"
        data[15409, 2] = 4.450e-20
        data[36381, 1] = 299.19812

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(N-compounds)/Nitriles/C6H5CN_Rajasekhar(2022)_298K_112.85-295.87nm.txt"
        data = [
            data[1:1204, :]
            data[1206:end, :]
        ]

    elseif typeof(data[end, 1]) <: SubString{String}
        if occursin("*)", data[end, 1])
            data = data[1:end-1, :]
        end

    elseif data[end, 2] == "<"
        data[end, 2] = data[end, 3]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O_GentieuMentall(1970)_298K_123.6,147nm.txt"
        data = [
            123.6 (1.62e-17 + 2.44e-17)/2
            147.0 5.28e-18
        ]
    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O_IUPAC(2013)_223K_251.7-350.0nm(rec,atm).txt"
        # this one is weird... has some fitted "Temperature coefficients Γ"
        λ = data[:,1]
        σ = data[:, 3]
        data = hcat(λ, σ)

    elseif data[1, 2] == "wavel."  && data[1,3] == "range"
        λ = data[:,1]
        σ = data[:, 3]
        data = hcat(λ, σ)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O_JPL-2015(2015)_298,223K_226-370nm(rec,atm).txt"
        # this is a weird one that has two temperature measurements. We should come back to this
        λ = data[:,1]
        σ = data[:,5]
        data = hcat(λ,σ)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O_JPL-2015(2015)_298K_226-375nm(rec,1nm).txt"
        # this is a weird one that has two temperature measurements. We should come back to this
        λ = data[:,1]
        σ = data[:,2]
        data = hcat(λ,σ)
    end

    if size(data, 1) > 1 && typeof(data[2,1]) <: SubString{String}
        if !isempty(data[2,1]) && isnumeric(data[2,1][1]) && occursin("-", data[2,1])
            data = data[2:end, :]

            λ1 = [parse(Float64, split(elem, "-")[1]) for elem ∈ data[:,1]]
            λ2 = [parse(Float64, split(elem, "-")[2]) for elem ∈ data[:,1]]

            data[:,1] = (λ1 .+ λ2) ./ 2
        end
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Bicarbonyls/FC(O)CF2CF2C(O)F_McGillen(2020)_296K_190-264nm.txt"
        data = [
            data[1:4, :]
            data[6:end, :]
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Bicarbonyls/HCOCH=CHCHO_TangZhu(2005)_293K_193-351nm.txt"
        data = [
            193.0	6.88e-18 0.39e-18
            248.0	3.62e-19 0.69e-19
            308.0	6.0e-21 ""
            351.0	6.0e-21 ""
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Carbonyl%20oxides/CH2OO_Ting(2014)_295K_308.4,351.8nm.txt"
        data[4, 1] = 351.8
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Halogenated%20aldehydes/COHF_Meller(1992)_296K_200-266nm(1nm).txt"
        data[:,2] = parse.(Float64, replace.(data[:,2], "-"=>"e-"))
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Ketones,ketenes/CH3C(O)C2H5_IUPAC(2021)_296K_200-355nm(rec).txt"
        data = data[:, 1:2]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Ketones,ketenes/CH3C(O)CH3_GierczakBurkholder(1998,2005)_280K_215-349nm(calc).txt"
        data[135, 2] = 7.29e-24
        data = [
            data
            [215 1.68e-21 ""]
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_Brion(1998)_295K_345.00-830.00nm.txt"
        λ = data[:,1]
        σ = data[:, 2]

        data = hcat(λ, σ)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_JPL-2010(2011)_293-298K_121.6-827.5nm(rec).txt"
        data = hcat(data[:,1], data[:,4])
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Peroxides/Hydrogen%20peroxide%20H2O2/H2O2_Jorand(2000)_300K_190-350nm(sol).txt"
        data = data[1:end-1, :]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/SO2_McGeeBurris(1987)_295K_300-320nm.txt"
        data = [
            data[1:645, :]
            data[647:end, :]
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Organic%20S-compounds/CH3SCH2Cl_Copeland(2014)_298K_290-355nm.txt"
        data = data[2:end, :]
        λ1 = data[:,1]
        λ2 = data[:,3]
        σ = data[:,4]
        ϵ = data[:,5]

        data = hcat((λ1 .+ λ2) ./ 2, σ, ϵ)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Organic%20S-compounds/CH3SNO_IUPAC(2004)_298K_190-545nm(rec).txt"
        data = data[1:end-1, :]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkyl%20radicals/C2H5_Khamaganov(2007)_296K_216.4,220,222nm.txt"
        data = [
            216.4	3.0e-18 NaN
            220	3.0e-18 NaN
            222	2.5e-18	0.6e-18
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkenes,polyenes+radicals/Polyenes/CH2=CHCH=CH2_Hanf(2010)_298K_121.6,193nm.txt"
        data = [
            121.6	4.4e-17	0.2e-17
            193	6.0e-17	0.1e-17
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogen%20oxides/Cl%20oxides/ClO_Rigaud(1977)_298K_272.6-307.9nm.txt"
        data = [
            272.6	3.90e-18 NaN
            272.9	7.22e-18 NaN
	          274.6	3.77e-18 NaN
            275.0	7.41e-18 NaN
	          276.8	3.34e-18 NaN
            277.2	7.20e-18 NaN
	          279.2	3.19e-18 NaN
            279.6	7.03e-18 NaN
	          281.9	3.14e-18 NaN
            282.2	6.98e-18 NaN
	          284.9	3.00e-18 NaN
            285.2	6.21e-18 NaN
	          287.8	2.53e-18 NaN
            288.4	5.14e-18 NaN
	          291.2	2.12e-18 NaN
            291.8	4.16e-18 NaN
	          295.0	1.77e-18 NaN
            295.4	3.81e-18 NaN
	          298.9	1.33e-18 NaN
            299.3	2.61e-18 NaN
	          303.2	9.90e-19 NaN
            303.5	1.84e-18 NaN
	          307.3	6.40e-19 NaN
            307.9	1.03e-18 NaN
        ]
    end


    return data
end


function replace_weird_characters!(data, j)
    for i ∈ 1:size(data,1)
        if typeof(data[i,j]) <: SubString{String}

            val =  replace(data[i,j],
                           "**)" => "",
                           "*)" => "",
                           "\xdf" => "",
                           "ß" => "",
                           "1.80-19" => "1.80e-19",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkali%20compounds/Na%20compounds/NaO3_SelfPlane(2002)_300K_216.3-395.2nm.txt
                           "1.40-24" => "1.40e-24",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkanes/CH4_LucchesiniGozzini(2007)_294K_838.180-846.696nm.txt
                           ")e" => "e",               # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkyl%20radicals/i-C3H7,(CH3)2CH_AdachiBasco(1981)_298K_238nm(max).txt
                           "~"=>"",                   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkynes,polyynes+radicals/Alkyne%20radicals/C2H_Fahr(2003)_298K_235-260nm.txt
                           "3.83-20" => "3.8e-20",    # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Aromatic%20compounds/Aldehydes/C6H5CHO_ZhuCronin(2000)_294K_280,285,308nm.txt
                           "4.20e-20." => "4.20e-20", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Aromatic%20compounds/Heterocyclic/C9H7N_Leach(2018)_298K_116-330nm.txt
                           "4.2-22" => "4.2e-22",     # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogenated%20N-compounds(inorg)/ClONO2_Rowland(1976)_298K_186-460nm.txt
                           "1.0-25" => "1.0e-25",     # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Chloroalkanes/CCl4_PrahladKumar(1995)_300K_186-240nm.txt
                           "3.0-21" => "3.0e-21",     # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Chloroalkanes/CHCl3_Robbins(1976)_295K_174-226nm.txt
                           "<" => "",
                           "1.20-21" => "1.20e-21",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Freons-CFC(C,F,Cl)/CFCl3(CFC-11)_IUPAC(2008)_210K_186-230nm(rec).txt
                           "2.73-21" => "2.73e-21",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Freons-HCFC(C,H,F,Cl)/CF3CHFCl(HCFC-124)_JPL-97(1997)_298K_190-220nm(rec).txt
                           "9.483-26" => "9.483e-26", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Halons(C,H,F,Cl,Br)/CF2Br2(Halon-1202)_JPL-2015(2015)_296K_170-350nm(rec).txt
                           ">"=>"",                   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogeno-alkanes+radicals/Iodoalkanes(C,H,F,Cl,Br,I)/c-C6H11I...Cl_Enami(2005)_263K_440nm.txt
                           "1.4-21" => "1.4e-21",     # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Halogens+mixed%20halogens/Br2_Passchier(1967)_423K_300-750nm.txt
                           "0.004-19" => "0.004e-19", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Nitrogen%20oxides/N2O_RontuCarlon(2010)_243K_184.950-228.802nm.txt
                           "3.7-19" => "3.7e-19",     # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(N-compounds)/Peroxynitrates/CH3C(O)O2NO2_Stephens(1969)_298K_220-270nm.txt
                           "351.8*)" => "351.8",        # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Carbonyl%20oxides/CH2OO_Ting(2014)_295K_308.4,351.8nm.txt
                           "1.591-21" => "1.691e-21",  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Ketones,ketenes/CF3C(O)CF(CF3)2_Yu(2018)_298K_250.09-405.51nm(calc).txt
                           "1.293-20" => "1.29e-20",  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Organics%20(carbonyls)/Ketones,ketenes/CH3C(O)CH3_Hynes(1992)_288K_300-340nm.txt
                           "3.10-46" => "3.10e-46",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Oxygen/O4_ThalmanVolkamer(2013)_233K_343.95-630.1nm(max).txt
                           "1.863-23" => "1.86e-23",   # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_Axson(2011)_296K_365,405,455nm.txt
                           "5.70-22" => "5.70e-22",    # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_BassPaur(1981)_243K_261-345nm.txt
                           "0.024-21" => "0.024e-21",  # https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_ElHelou(2005)_175K_543.667-632.991nm(LASIM).txt
                           "4.52593-21"=>"4.5259e-21",  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_Serdyuchenko(2014)_223K_543.667-1046.766nm.txt
                           "3.283-19" => "3.28e-19",  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Ozone/O3_Yoshino(1993)_195K_185-250nm(absolute).txt
                           "8.50-19" => "8.50e-19",  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Peroxy%20radicals/Alkylperoxy+F+Cl/CF3CHClO2_Mogelberg(1995)_296K_225-290nm.txt
                           "1.74e-21)" => "1.7e-21", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/S2Cl2_Speth(1998)_298K_230-367nm.txt
                           "1.9894-19" => "1.989e-19", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Organic%20S-compounds/OCS_Limao-Vieira(2015)_298K_115-320nm.txt
                           "8-0e-21"=>"8.0e-21", # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Sulfur%20compounds/Organic%20S-compounds/OCS_Limao-Vieira(2015)_298K_115-320nm.txt
                           )
            if data[i,j] != ""
                try
                    data[i,j] = parse(Float64, val)
                catch
                end
            end

            if data[i,j] == "*)"
                data[i,j] = NaN
            end

        end
    end
end


function replace_blank_with_nan!(data, j)
    for i ∈ 1:size(data, 1)
        if data[i,j] == ""
            data[i,j] = NaN
        elseif data[i,j] == "-"
            data[i,j] = NaN
        elseif data[i,j] == "--"
            data[i,j] = NaN
        end
    end
end



function interpret_data(data::Matrix)
    # make sure we have at least 2 columns
    @assert size(data, 2) > 1

    # start data when first column becomes a number
    i_start = findfirst([typeof(elem) <: Number for elem ∈data[:,1]])
    data = data[i_start:end, :]


    # deal with weird lines with zeros as in https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkali%20compounds/Li,%20K,%20Rb,%20Cs%20compounds/CsBr_DavidovitsBrodhead(1967)_1071K_213,256,275nm(max).txt
    idxs = [elem != 0 for elem ∈ data[:,1]]
    data = data[idxs, :]

    # deal with weird error limits as in https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkali%20compounds/Li,%20K,%20Rb,%20Cs%20compounds/CsCl_DavidovitsBrodhead(1967)_1123K_247.5nm(max).txt
    idxs =  [!(data[i,1] < 0 && all([data[i, j] == "" for j∈2:size(data,2)])) for i ∈ 1:size(data, 1)]
    data = data[idxs, :]


    # replace weird accidental characters
    for j ∈ 1:size(data,2)
        try
            replace_weird_characters!(data, j)
        catch
            continue
        end
    end

    if size(data, 2) == 2
        replace_blank_with_nan!(data, 2)  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkanes/CH4_LucchesiniGozzini(2007)_294K_838.180-846.696nm.txt

        λ = Float64.(data[:, 1])
        σ = Float64.(data[:, 2])
        Δσ = Float64.([NaN for i∈1:size(data,1)])
    elseif size(data, 2) == 3 && data[1,3] == ""

        λ = Float64.(data[:, 1])
        σ = Float64.(data[:, 2])
        Δσ = Float64.([NaN for i∈1:size(data,1)])
    else
        replace_blank_with_nan!(data, 2)  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkenes,polyenes+radicals/Alkenes/CH2=CH2_LucchesiniGozzini(2011)_294K_825.303-837.215nm.txt 
        replace_blank_with_nan!(data, 3)  # see https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkyl%20radicals/C2H5_Khamaganov(2007)_296K_216.4,220,222nm.txt
        λ = Float64.(data[:,1])
        σ = Float64.(data[:,2])
        Δσ = Float64.(data[:,3])
    end

    return λ, σ, Δσ
end



# -----------------------------------------------------------------------
# try it out!
# -----------------------------------------------------------------------
ignore_urls = String[
    "",
]

T_vs_σ_urls = String[
    "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/cross_sections/Alkanes+alkyl%20radicals/Alkyl%20radicals/CH3_Macpherson(1985)_298-537K_216.36nm.txt",
]


# function test_ingestion()
#     @showprogress 1 for i ∈ 1:length(url_list)
# #    @showprogress 1 for i ∈ 7338:length(url_list)
#         if !(url_list[i] ∈ ignore_urls) && !(url_list[i] ∈ T_vs_σ_urls)
#             try
#                 data = download_data(url_list[i])
#                 λ,σ,Δσ = interpret_data(data)

#             catch e
#                 println("FAILED")
#                 println(e)
#                 return i, url_list[i]
#             end
#         end
#     end
#     return length(url_list), ""
# end

# idx, failed_url = test_ingestion()
# failed_url

# ""
# ""

# ""

# data = download_data(failed_url)
# interpret_data(data)


# i_start = findfirst([typeof(elem) <: SubString{String} for elem ∈data[:,2]])
# data[([typeof(elem) <: SubString{String} for elem ∈data[:,2]]), :]
# data[i_start+1, :]
# typeof.(data[i_start, :])

# now go through and download all the data
failures = String[]

Threads.nthreads()
summary_df = CSV.File(summary_list[1]) |> DataFrame

p = Progress(length(summary_list))
Threads.@threads for i ∈ 1:length(summary_list)
    summary_df = CSV.File(summary_list[i]) |> DataFrame
    for row ∈ eachrow(summary_df)
        try
            data = download_data(row.download_url)
            λ,σ,Δσ = interpret_data(data)
            out_df = DataFrame(λ=λ, σ=σ, Δσ=Δσ)
            CSV.write(joinpath(basepath, row.fname), out_df)
        catch e
            push!(failures, row.download_url)
        end
    end
    next!(p)
end

# okay let's load in some test data for formaldehyde
# path_to_CH2O = "/home/jwaczak/gitrepos/ActivePure/Photolysis.jl/data/cross-sections/Organics (carbonyls)/Aldehydes(aliphatic)/CH2O"
# summary_df = CSV.File(joinpath(path_to_CH2O, "summary.csv")) |> DataFrame


# res_df = summary_df[(summary_df.T1 .< 250.0) .&& isnan.(summary_df.T2),:]
# cs = cgrad(:inferno, res_df.T1)
# p = plot()
# #plot!(p, df.λ, df.σ) #, linecolor=res_df.T1[1])
# for i ∈ 1:10
#     df = CSV.File(joinpath("../", res_df.fname[i])) |> DataFrame
#     plot!(p, df.λ, df.σ, linecolor=cs[i], c=cs, label="", colobar=true, show=true)
# end

# h2 = scatter([0,0], [0,1], zcolor=[minimum(res_df.T1),maximum(res_df.T1)],
#              xlims=(1,1.1), xshowaxis=false, yshowaxis=false, label="", c=:inferno, colorbar_title="T [K]", grid=false)
# l = @layout [a b{0.01w}]

# ylims!(p, 0, 1e-19)
# xlims!(p, 200, 400)

# plot(p, h2, layout=l)


