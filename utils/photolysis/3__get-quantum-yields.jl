using DelimitedFiles
using CSV, DataFrames
using Unitful
using ProgressMeter

# 0. Set up paths
basepath = "../../assets/mpi-mainz-uvviz/raw"
crosssections_path = joinpath(basepath, "cross-sections")
quantumyields_path = joinpath(basepath, "quantum-yields")

isdir(basepath)
isdir(crosssections_path)
isdir(quantumyields_path)


# 1. Get list of all summary.csv files
summary_list = String[]


for (root, dirs, files) in walkdir(quantumyields_path)
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


function data_and_summary(url)
    return url, replace(url, "uvvis_data" => "uvvis")
end


function download_data(url)
    tmp = download(url)
    data = readdlm(tmp)

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/Cl2O%7BCl+ClO%7D_Nelson(1994)_298K_193,248,308nm.txt"
        # I'm not sure what "minor" and "major" are supposed to be
        data = [
            308.0 1.0 NaN :Cl
            308.0 1.0 NaN :ClO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/Cl2O%7BCl+ClO%7D_SanderFriedl(1989)_220-400K_180-400nm.txt"
        data = [
            180.0 0.25 0.05 :Cl
            180.0 0.25 0.05 :ClO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/Cl2O%7BO+Cl2%7D_SanderFriedl(1989)_220-400K_180-400nm.txt"
        λ = collect(180.0:1:400.0)
        Φ = 0.25 * ones(length(λ))
        ΔΦ = 0.05 * ones(length(λ))
        species = [:O for _ ∈ 1:length(λ)]
        species2 = [:Cl₂ for _ ∈ 1:length(λ)]
        return vcat(λ, λ), vcat(Φ, Φ), vcat(ΔΦ, ΔΦ), vcat(species, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClO+ClO%7D_Moore(1999)_298K_248-308nm.txt"
        data = [
            248.0 0.12 0.07 :ClO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClO+ClO%7D_Plenge(2004)_170K_250-308nm.txt"
        data = [
            250.0 0.02 NaN :ClO
            308.0 0.10 NaN :ClO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClOO+Cl%7D_CoxHayman(1988)_203-233K_254nm.txt"
        data = [
            245.0 1.00 NaN :ClO2
            245.0 1.00 NaN :Cl
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClOO+Cl%7D_Molina(1990)_235K_308nm.txt"
        data = [
            308.0 1.03 0.12 :ClO2
            308.0 1.03 0.12 :Cl
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClOO+Cl%7D_Moore(1999)_298K_248-308nm.txt"
        data = [
            248.0 0.88 0.07 :ClO2
            308.0 0.90 0.1 :ClO2
            248.0 0.88 0.07 :Cl
            308.0 0.90 0.1 :Cl
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/ClOOCl%7BClOO+Cl%7D_Plenge(2004)_170K_250-308nm.txt"
        data = [
            245.0 0.98 NaN :ClO2
            308.0 0.90 NaN :ClO2
            245.0 0.98 NaN :Cl
            308.0 0.90 NaN :Cl
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogen%20oxides/Cl%20oxides/OClO%7BCl+O2%7D_Bishenden(1991)_298K_362nm.txt"
        data = [
            362.0 0.15 0.06 :Cl
            362.0 0.15 0.06 :O2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/BrONO2%7BBr+NO3%7D_Harwood(1998)_298K_248-352.5nm.txt"
        data = [
            248.0 0.28 0.09 :NO3
            305.0 1.01 0.35 :NO3
            352.5 0.92 0.43 :NO3
            248.0 0.5 NaN :Br
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/BrONO2%7BBr+NO3%7D_Soller(2002)_298K_248-355nm.txt"
        data = [
            248.0 0.35 0.08 :Br
            266 0.65 0.14 :Br
            308.0 0.62 0.11 :Br
            355.0 0.77 0.19 :Br
            248.0 0.66 0.15 Symbol("O(3P)")
            266.0 0.18 0.04 Symbol("O(3P)")
            308.0 0.13 0.03 Symbol("O(3P)")
            355.0 0.02 NaN Symbol("O(3P)")
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/BrONO2%7BBrO+NO2%7D_Harwood(1998)_298K_248nm.txt"
        data = [
            248.0 0.5 NaN :BrO
            248.0 0.5 NaN :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/BrONO2%7BBrO+NO2%7D_Soller(2002)_298K_266,355nm.txt"
        data = [
            266.0 0.37 0.12 :BrO
            355.0 0.23 0.08 :BrO
            266.0 0.37 0.12 :NO2
            355.0 0.23 0.08 :NO2
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClNO2%7BCl+NO2%7D_Carter(1999)_298K_235nm.txt"
        data = [
            235.0 0.15 0.05 :Cl
            235.0 0.15 0.05 :NO2
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClNO2%7BCl+NO2%7D_NelsonJohnston(1981)_298K_350nm.txt"
        data = [
            350.0 0.93 0.15 :Cl
            350.0 0.93 0.15 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClNO2%7BCl+NO2x%7D_Carter(1999)_298K_235nm.txt"
        data = [
            235.0 0.85 0.05 :Cl
            235.0 0.85 0.05 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClNO2%7BClNO+O%7D_NelsonJohnston(1981)_298K_350nm.txt"
        data = [
            350.0 0.02 NaN :ClNO
            350.0 0.02 NaN :O
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Burrows(1988)_298K_254nm.txt"
        # not sure how to treat "high" vs "low" ClONO₂
        data = [
            254.0 1.04 0.04 :Cl
            254.0 1.04 0.04 :NO3
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Chang(1979)_298K_270nm.txt"
        data = [
            270.0 1.0 0.2 :Cl
            270.0 1.0 0.2 :NO3
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Goldfarb(1997)_298K_193-308nm.txt"
        data = [
            193.0 0.53 0.10 :Cl
            222.0 0.46 0.10 :Cl
            248.0 0.41 0.13 :Cl
            308.0 0.64 0.20 :Cl
            193.0 0.53 0.10 :NO3
            222.0 0.46 0.10 :NO3
            248.0 0.41 0.13 :NO3
            308.0 0.64 0.20 :NO3
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Margitan(1983)_298K_266-355nm.txt"
        data = [
            266.0 0.9 0.1 :Cl
            355.0 0.9 0.1 :Cl
            266.0 0.9 0.1 :NO3
            355.0 0.9 0.1 :NO3
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_MarinelliJohnston(1979)_298K_249nm.txt"
        data = [
            249.0 0.55 NaN :Cl
            249.0 0.55 NaN :NO3
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Minton(1992)_298K_193,248nm.txt"
        data = [
            193.0 0.64 0.08 :Cl
            248.0 0.54 0.08 :Cl
            193.0 0.64 0.08 :NO3
            248.0 0.54 0.08 :NO3
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Moore(1995)_298K_308nm.txt"
        data = [
            308.0 0.67 0.06 :Cl
            308.0 0.67 0.06 :NO3
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Nelson(1996)_298K_193,248nm.txt"
        data = [
            193.0 0.66 0.08 :Cl
            248.0 0.54 0.08 :Cl
            193.0 0.66 0.08 :Cl
            248.0 0.54 0.08 :Cl
        ]
    end



    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Nickolaisen(1996)_298K_200-350nm.txt"
        λ1 = 200.0:1.0:300.0
        λ2 = 300.0:1.0:350.0

        ϕ1 = 0.39*ones(length(λ1))
        ϕ2 = 0.56*ones(length(λ2))

        Δϕ1 = 0.20*ones(length(λ1))
        Δϕ2 = 0.08*ones(length(λ2))

        λ = vcat(λ1, λ2)
        ϕ = vcat(ϕ1, ϕ2)
        Δϕ = vcat(Δϕ1, Δϕ2)
        species1 = [:Cl for _ ∈ 1:(length(λ))]
        species2 = [:NO₃ for _ ∈ 1:(length(λ))]

        return vcat(λ, λ), vcat(ϕ, ϕ), vcat(Δϕ, Δϕ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Tyndall(1997)_298K_308nm.txt"
        data = [
            308.0 0.80 0.08 :Cl
            308.0 0.80 0.08 :NO3
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BCl+NO3%7D_Yokelson(1997)_298K_193-308nm.txt"
        data = [
            193.0 0.18 0.04 :NO3
            248.0 0.60 0.09 :NO3
            308.0 0.67 0.09 :NO3
            352.5 0.93 0.24 :NO3
            193.0 0.45 0.08 :Cl
            248.0 0.60 0.12 :Cl
            308.0 0.73 0.14 :Cl
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Goldfarb(1997)_298K_193-308nm.txt"
        data = [
            193.0 0.29 0.20 :ClO
            222.0 0.64 0.20 :ClO
            248.0 0.39 0.19 :ClO
            308.0 0.37 0.19 :ClO
            193.0 0.29 0.20 :NO2
            222.0 0.64 0.20 :NO2
            248.0 0.39 0.19 :NO2
            308.0 0.37 0.19 :NO2
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Minton(1992)_298K_193,248nm.txt"
        data = [
            193.0 0.36 0.08 :ClO
            248.0 0.46 0.08 :ClO
            193.0 0.36 0.08 :NO2
            248.0 0.46 0.08 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Moore(1995)_298K_308nm.txt"
        data = [
            308.0 0.33 0.06 :ClO
            308.0 0.33 0.06 :NO2
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Nelson(1996)_298K_193,248nm.txt"
        data = [
            193.0 0.34 0.08 :ClO
            248.0 0.46 0.08 :ClO
            193.0 0.34 0.08 :NO2
            248.0 0.46 0.08 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Nickolaisen(1996)_298K_200-350nm.txt"
        λ1 = 200.0:1.0:300.0
        λ2 = 300.0:1.0:350.0

        ϕ1 = 0.61*ones(length(λ1))
        ϕ2 = 0.44*ones(length(λ2))

        Δϕ1 = 0.20*ones(length(λ1))
        Δϕ2 = 0.08*ones(length(λ2))

        λ = vcat(λ1, λ2)
        ϕ = vcat(ϕ1, ϕ2)
        Δϕ = vcat(Δϕ1, Δϕ2)
        species1 = [:ClO for _ ∈ 1:(length(λ))]
        species2 = [:NO₃ for _ ∈ 1:(length(λ))]

        return vcat(λ, λ), vcat(ϕ, ϕ), vcat(Δϕ, Δϕ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BClO+NO2%7D_Tyndall(1997)_298K_308nm.txt"
        data = [
            308.0 0.28 0.12 :ClO
            308.0 0.28 0.12 :NO2
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BO(3P)+ClONO%7D_Burrows(1988)_298K_254nm.txt"
        data = [
            254.0 0.24 NaN Symbol("O(3P)")
            254.0 0.2 NaN :ClONO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BO(3P)+ClONO%7D_Goldfarb(1997)_298K_193-308nm.txt"
        data = [
            193.0 0.73 0.08 Symbol("O(3P)")
            222.0 0.17 0.08 Symbol("O(3P)")
            248.0 0.10 NaN Symbol("O(3P)")
            308.0 0.05 NaN Symbol("O(3P)")
            193.0 0.73 0.08 :ClONO
            222.0 0.17 0.08 :ClONO
            248.0 0.10 NaN :ClONO
            308.0 0.05 NaN :ClONO
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BO(3P)+ClONO%7D_Margitan(1983)_298K_266-355nm.txt"
        data = [
            266.0 0.1 NaN Symbol("O(3P)")
            355.0 0.1 NaN Symbol("O(3P)")
            266.0 0.1 NaN :ClONO
            355.0 0.1 NaN :ClONO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BO(3P)+ClONO%7D_Tyndall(1997)_298K_308nm.txt"
        data = [
            308.0 0.05 NaN Symbol("O(3P)")
            308.0 0.05 NaN :ClONO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/ClONO2%7BO(3P)+ClONO%7D_Yokelson(1997)_298K_193,248nm.txt"
        data = [
            193.0 0.9 NaN Symbol("O(3P)")
            248.0 0.4 NaN Symbol("O(3P)")
            193.0 0.9 NaN :ClONO
            248.0 0.4 NaN :ClONO
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/IONO2%7BI+NO3%7D_Joseph(2007)_298K_248nm.txt"
        data = [
            248.0 0.21 0.09 :I
            248.0 0.21 0.09 :NO3
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Halogenated%20N-compounds(inorg)/IONO2%7BIO+NO2%7D_Joseph(2007)_298K_248nm.txt"
        data = [
            248.0 0.02 NaN :IO
            248.0 0.02 NaN :IO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO+O(3P)%7D_Barker(1985)_298K_290nm.txt"
        data = [
            290.0 0.1 NaN :NO3
            290.0 0.1 NaN :NO
            290.0 0.1 NaN Symbol("O(3P)")
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO+O(3P)%7D_Ravishankara(1986)_298K_248-289nm.txt"
        data = [
            248.0 0.77 0.10 :NO3
            248.0 0.72 0.09 Symbol("O(3P)")
            266.0 0.38 0.08 Symbol("O(3P)")
            287.0 0.21 0.03 Symbol("O(3P)")
            289.0 0.15 0.03 Symbol("O(3P)")
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO2%7D_Barker(1985)_298K_290nm.txt"
        data = [
            290.0 0.8 0.2 :NO3
            290.0 0.8 0.2 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO2%7D_Burrows(1984)_298K_254nm.txt"
        data = [
            254.0 0.75 NaN :NO3
            254.0 0.75 NaN :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO2%7D_Harwood(1998)_298K_248-352.5nm.txt"
        data = [
            248.0 0.64 0.20 :NO3
            308.0 0.96 0.15 :NO3
            352.5 1.03 0.15 :NO3
            248.0 0.64 0.20 :NO2
            308.0 0.96 0.15 :NO2
            352.5 1.03 0.15 :NO2
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Davenport(1978)_223K_390-420nm.txt"
        data = [
            389.8 0.74 NaN :NO
            395.1 0.79 NaN :NO
            399.4 0.56 NaN :NO
            405.2 0.26 NaN :NO
            410.2 0.041 NaN :NO
            419.7 0.0053 NaN :NO
            389.8 0.74 NaN Symbol("O(3P)")
            395.1 0.79 NaN Symbol("O(3P)")
            399.4 0.56 NaN Symbol("O(3P)")
            405.2 0.26 NaN Symbol("O(3P)")
            410.2 0.041 NaN Symbol("O(3P)")
            419.7 0.0053 NaN Symbol("O(3P)")
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Davenport(1978)_300K_390-421nm.txt"
        data = [
            390.1 0.75 NaN :NO
            395.2 0.82 NaN :NO
            399.8 0.64 NaN :NO
            405.2 0.33 NaN :NO
            410.1 0.12 NaN :NO
            420.8 0.018 NaN :NO
            390.1 0.75 NaN Symbol("O(3P)")
            395.2 0.82 NaN Symbol("O(3P)")
            399.8 0.64 NaN Symbol("O(3P)")
            405.2 0.33 NaN Symbol("O(3P)")
            410.1 0.12 NaN Symbol("O(3P)")
            420.8 0.018 NaN Symbol("O(3P)")
        ]
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_GaedtkeTroe(1975)_296K_285-424nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]

        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Gardner(1987)_298K_285-424nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Harker(1977)_298K_375-420nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_JonesBayes(1973)_300K_299.5-445.8nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Roehl(1994)_248K_390-426nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Roehl(1994)_298K_390-428nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Troe(2000)_248K_300-415nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO2%7BNO+O(3P)%7D_Troe(2000)_298K_300-415nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [Symbol("O(3P)") for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_Johnston(1996)_190K_586-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_Johnston(1996)_230K_586-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_Johnston(1996)_298K_586-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_MagnottaJohnston(1980)_296K_482-632nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_Orlando(1993)_298K_586-633nm.txt"
        λ = data[1:end-1, 1]
        Φ = data[1:end-1, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO+O2%7D_Wayne(1991)_296K_585-639nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_GrahamJohnston(1978)_298K_510-700nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_GrahamJohnston(1978)_329K_510-700nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_Johnston(1996)_190K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_Johnston(1996)_230K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_Johnston(1996)_298K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_MagnottaJohnston(1980)_296K_472-632nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_Orlando(1993)_298K_586-637nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7BNO2+O%7D_Wayne(1991)_296K_580-639nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:O for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_Johnston(1996)_190K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:NO2 for _ ∈ 1:length(λ)]
        species2 = [:total for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_Johnston(1996)_230K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_Johnston(1996)_298K_585-640nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_MagnottaJohnston(1980)_296K_472-632nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_Orlando(1993)_298K_586-637nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/NO3%7Btotal%7D_Wayne(1991)_296K_580-639nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_Clark(1978)_294K_299.1-353.5nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_Moortgat(1983)_220K_270-353.1nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_Moortgat(1983)_300K_253-354nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_MoortgatWarneck(1979)_300K_276.7-355nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_Roeth(2015)_300K_250-360nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_SperlingToby(1973)_373K_313-366nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7B(total)%7D_Turco(1975)_298K_280-400nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species = [:total for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Calvert(1972)_298K_290-360nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Clark(1978)_294K_299.1-339.2nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Gorrotxategi-Carbajo(2008)_294K_303.70-329.51nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_HorowitzCalvert(1978)_298K_289-338nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Lewis(1976)_296K_284-339nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_LewisLee(1978)_296K_275.4,288.2,303.5nm.txt"
        data = [
            275.4 0.68 0.10 :H
            288.2 0.64 0.10 :H
            303.5 0.68 0.05 :H
            275.4 0.68 0.10 :HCO
            288.2 0.64 0.10 :HCO
            303.5 0.68 0.05 :HCO
        ]
    end


    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Marling(1977)_298K_304-337.8nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)
    end

    if url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_McQuiggCalvert(1969)_353pm20K_279.6-369.4nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url =="https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Moortgat(1983)_220K_270.6-353.1nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Moortgat(1983)_298K_253-353nm.txt"
        data = [
            253.0 0.317 0.074 :H
            262.1 0.305 0.055 :H
            270.6 0.385 0.055 :H
            277.4 0.479 0.095 :H
            284.9 0.666 0.071 :H
            295.0 0.796 0.216 :H
            304.2 0.750 0.151 :H
            315.8 0.823 0.114 :H
            317.0 0.666 0.153 :H
            326.8 0.355 0.106 :H
            339.8 0.009 0.002 :H
            354.0 0.004 0.001 :H
            253.0 0.317 0.074 :HCO
            262.1 0.305 0.055 :HCO
            270.6 0.385 0.055 :HCO
            277.4 0.479 0.095 :HCO
            284.9 0.666 0.071 :HCO
            295.0 0.796 0.216 :HCO
            304.2 0.750 0.151 :HCO
            315.8 0.823 0.114 :HCO
            317.0 0.666 0.153 :HCO
            326.8 0.355 0.106 :HCO
            339.8 0.009 0.002 :HCO
            354.0 0.004 0.001 :HCO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_MoortgatWarneck(1979)_298K_276-355nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Pope(2005)_294K_308.86,314.13nm.txt"
        data = [
            308.86 0.71 NaN :H
            314.13 0.69 NaN :H
            308.86 0.71 NaN :HCO
            314.13 0.69 NaN :HCO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Roeth(2015)_300K_250-360nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Smith(2002)_298K_268.75-338.75nm.txt"
        data = [
            268.75 0.41 0.06 :H
            278.75 0.55 0.06 :H
            283.75 0.65 0.07 :H
            288.75 0.72 0.07 :H
            293.75 0.67 0.07 :H
            298.75 0.62 0.07 :H
            301.75 0.70 0.11 :H
            303.75 0.753 NaN :H
            306.25 0.66 0.07 :H
            308.75 0.71 0.08 :H
            311.25 0.68 0.09 :H
            313.75 0.69 0.07 :H
            316.25 0.65 0.07 :H
            321.25 0.64 0.08 :H
            326.25 0.51 0.05 :H
            328.75 0.36 0.04 :H
            331.25 0.46 0.06 :H
            333.75 0.30 0.07 :H
            336.25 0.07 0.01 :H
            338.75 0.04 0.01 :H
            268.75 0.41 0.06 :HCO
            278.75 0.55 0.06 :HCO
            283.75 0.65 0.07 :HCO
            288.75 0.72 0.07 :HCO
            293.75 0.67 0.07 :HCO
            298.75 0.62 0.07 :HCO
            301.75 0.70 0.11 :HCO
            303.75 0.753 NaN :HCO
            306.25 0.66 0.07 :HCO
            308.75 0.71 0.08 :HCO
            311.25 0.68 0.09 :HCO
            313.75 0.69 0.07 :HCO
            316.25 0.65 0.07 :HCO
            321.25 0.64 0.08 :HCO
            326.25 0.51 0.05 :HCO
            328.75 0.36 0.04 :HCO
            331.25 0.46 0.06 :HCO
            333.75 0.30 0.07 :HCO
            336.25 0.07 0.01 :HCO
            338.75 0.04 0.01 :HCO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_SperlingToby(1973)_373K_313,334,366nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Tang(1979)_296K_294-321.5nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_Calvert(1972)_298K_290-360nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H2 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_Clark(1978)_294K_299.1-353.5nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H2 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_McQuiggCalvert(1969)_353pm20K_279.6-369.4nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H2 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_Moortgat(1983)_220K_253-353nm.txt"
        data = [
            270.6 0.388 0.111 :H2
            278.1 0.283 0.059 :H2
            285.6 0.266 0.033 :H2
            294.9 0.272 0.030 :H2
            304.4 0.213 0.027 :H2
            315.3 0.277 0.020 :H2
            327.0 0.717 0.010 :H2
            339.3 0.509 0.036 :H2
            353.1 0.124 0.038 :H2
            270.6 0.388 0.111 :CO
            278.1 0.283 0.059 :CO
            285.6 0.266 0.033 :CO
            294.9 0.272 0.030 :CO
            304.4 0.213 0.027 :CO
            315.3 0.277 0.020 :CO
            327.0 0.717 0.010 :CO
            339.3 0.509 0.036 :CO
            353.1 0.124 0.038 :CO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_Moortgat(1983)_298K_253-353nm.txt"
        data = [
            253.0 0.483 0.071 :H2
            262.0 0.495 0.061 :H2
            270.6 0.477 0.062 :H2
            277.4 0.350 0.067 :H2
            284.9 0.308 0.021 :H2
            295.0 0.287 0.043 :H2
            304.2 0.216 0.036 :H2
            315.8 0.290 0.030 :H2
            317.0 0.414 0.058 :H2
            326.8 0.660 0.109 :H2
            339.8 0.639 0.154 :H2
            354.0 0.252 0.044 :H2
            253.0 0.483 0.071 :CO
            262.0 0.495 0.061 :CO
            270.6 0.477 0.062 :CO
            277.4 0.350 0.067 :CO
            284.9 0.308 0.021 :CO
            295.0 0.287 0.043 :CO
            304.2 0.216 0.036 :CO
            315.8 0.290 0.030 :CO
            317.0 0.414 0.058 :CO
            326.8 0.660 0.109 :CO
            339.8 0.639 0.154 :CO
            354.0 0.252 0.044 :CO
        ]



    elseif url ==  "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_MoortgatWarneck(1979)_298K_276-355nm.txt"
        data = [
            276.7 0.245 0.063 :H2
            284.1 0.232 NaN :H2
            295.0 0.240 0.067 :H2
            303.9 0.116 0.046 :H2
            316.3 0.240 NaN :H2
            317.0 0.347 0.044 :H2
            326.7 0.576 0.110 :H2
            340.2 0.572 0.150 :H2
            355.0 0.256 0.067 :H2
            276.7 0.245 0.063 :CO
            284.1 0.232 NaN :CO
            295.0 0.240 0.067 :CO
            303.9 0.116 0.046 :CO
            316.3 0.240 NaN :CO
            317.0 0.347 0.044 :CO
            326.7 0.576 0.110 :CO
            340.2 0.572 0.150 :CO
            355.0 0.256 0.067 :CO
        ]


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_Roeth(2015)_300K_250-360nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H2 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH2+CO%7D_SperlingToby(1973)_353-392K_313,334,366nm.txt"
        data = [
            313 0.20 0.2 :H2
            334 0.25 0.1 :H2
            366 0.45 0.1 :H2
            313 0.20 0.2 :CO
            334 0.25 0.1 :CO
            366 0.45 0.1 :CO
        ]


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3+CHO%7D_Horowitz(1986)_298K_300nm(zero%20press).txt"
        data = [
            300 0.92 NaN :CH3
            300 0.92 NaN :CHO
        ]


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3+CHO%7D_HorowitzCalvert(1982)_298K_290-331.2nm(1atm).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH3 for _ ∈ 1:length(λ)]
        species2 = [:CHO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3+CHO%7D_HorowitzCalvert(1982)_298K_290-331.2nm(zero%20press).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH3 for _ ∈ 1:length(λ)]
        species2 = [:CHO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3+CHO%7D_Weaver(1976-77)_298K_313nm.txt"
        data = [
            313.0 0.15 NaN :CH3
            313.0 0.15 NaN :CHO
        ]


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3CO+H%7D_Horowitz(1986)_298K_300nm(zero%20press).txt"
        data = [
            300 0.92 NaN :CH3CO
            300 0.92 NaN :H
        ]


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3CO+H%7D_HorowitzCalvert(1982)_298K_290-331.2nm(1atm).txt"
        data = [
            290.0 0.026 NaN :CH3CO
            300.0 0.009 NaN :CH3CO
            313.0 0.002 NaN :CH3CO
            320.0 0.00 NaN :CH3CO
            331.2 0.00 NaN :CH3CO
            290.0 0.026 NaN :H
            300.0 0.009 NaN :H
            313.0 0.002 NaN :H
            320.0 0.00 NaN :H
            331.2 0.00 NaN :H
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH3CO+H%7D_HorowitzCalvert(1982)_298K_290-331.2nm(zero%20press).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH3CO for _ ∈ 1:length(λ)]
        species2 = [:H for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4%7D_Meyrahn(1984)_298K_257-327nm(air).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH4 for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species1

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4%7D_Meyrahn(1984)_298K_257-327nm(N2).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH4 for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species1


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4+CO%7D_Horowitz(1986)_298K_300nm(zero%20press).txt"
        data = [
            300 0.015 NaN :CH4
            300 0.015 NaN :CO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4+CO%7D_HorowitzCalvert(1982)_298K_290-331.2nm(1atm).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH4 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4+CO%7D_HorowitzCalvert(1982)_298K_290-331.2nm(zero%20press).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CH4 for _ ∈ 1:length(λ)]
        species2 = [:CO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCH4+CO%7D_ParmenterNoyes(1963)_303K_257.3nm(1atm).txt"
        data = [
            253.7 0.63 NaN :CH4
            253.7 0.63 NaN :CO
        ]

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCO%7D_Meyrahn(1984)_298K_257-327nm(air).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CO for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species1


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH3CHO%7BCO%7D_Meyrahn(1984)_298K_257-327nm(N2).txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:CO for _ ∈ 1:length(λ)]
        return λ, Φ, ΔΦ, species1

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Armerding(1995)_298K_300-355nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Ball(1997)_298K_300-338nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_BallHancock(1995)_227K_300-322nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_BallHancock(1995)_298K_300-325nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Bauer(2000)_273K_304-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Bauer(2000)_295K_305-375nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Bauer(2000)_295K_305-375nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_BrockWatson(1980)_298K_297.5-325nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Cooper(1993)_295K_221.5-243.5nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_203K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_223K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_253K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_273K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_298K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Matsumi(2002)_321K_220-340nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Silvente(1997)_295K_301-336nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Smith(2000)_226K_312-337nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Smith(2000)_263K_312-337nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Smith(2000)_298K_295-338nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Takahashi(1996)_298K_308.07-325.90nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Takahashi(1998)_227K_304.56-313.53nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Takahashi(1998)_295K_304.19-328.95nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_203K_306-329nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_223K_306-329nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_253K_306-328nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_273K_306-329nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_298K_306-328nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)


    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Talukdar(1998)_321K_306-329nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_Taniguchi(2000)_295K_297-308nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Ozone/O3%7BO(1D)+O2%7D_TrolierWiesenfeld(1988)_295K_274-325nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [Symbol("O(1D)") for _ ∈ 1:length(λ)]
        species2 = [:O2 for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Organics%20(carbonyls)/Aldehydes(aliphatic)/CH2O%7BH+HCO%7D_Turco(1975)_298K_280-400nm.txt"
        λ = data[:, 1]
        Φ = data[:, 2]
        ΔΦ = [NaN for _ ∈ 1:length(λ)]
        species1 = [:H for _ ∈ 1:length(λ)]
        species2 = [:HCO for _ ∈ 1:length(λ)]
        return vcat(λ,λ), vcat(Φ,Φ), vcat(ΔΦ, ΔΦ), vcat(species1, species2)

    elseif url == "https://uv-vis-spectral-atlas-mainz.org//uvvis_data/quantum_yields/Nitrogen%20oxides/N2O5%7BNO3+NO2%7D_Swanson(1984)_298K_249-350nm.txt"
        data = [
            248.0 0.89 0.15 :NO3
            350.0 0.84 0.09 :NO3
            248.0 0.89 0.15 :NO2
            350.0 0.84 0.09 :NO2
        ]


    end

    λ = data[:, 1]
    Φ = data[:, 2]
    ΔΦ = data[:,3]
    species = data[:,4]
    return λ, Φ, ΔΦ, species
end


# data_url, url_page = data_and_summary(url_list[154])
# data_url
# url_page

# length(url_list)

# λ,Φ,ΔΦ,species = download_data(data_url)


summary_list

failed = String[]
p = Progress(length(summary_list))
Threads.@threads for i ∈ 1:length(summary_list)
    summary_df = CSV.File(summary_list[i]) |> DataFrame
    for row ∈ eachrow(summary_df)
        try
            λ,Φ,ΔΦ,species = download_data(row.download_url)
            out_df = DataFrame(λ=λ, Φ=Φ, ΔΦ=ΔΦ, species=species)
            CSV.write(joinpath(basepath, row.fname), out_df)
        catch e
            println("FAILED  ", row.download_url)
            println(e)
            push!(failed, row.download_url)
        end
    end
    next!(p)
end
