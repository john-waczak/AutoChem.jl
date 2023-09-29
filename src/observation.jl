function Obs!(h, u, idx_meas, idx_pos, idx_neg)
    n = length(idx_meas)
    h[1:n] = u[idx_meas]
    h[n+1] = sum(u[idx_pos])
    h[n+2] = sum(u[idx_neg])
end

function Obs(u, idx_meas, idx_pos, idx_neg)
    n = length(idx_meas)
    return vcat(u[idx_meas], sum(u[idx_pos]), sum(u[idx_neg]))
end


function JObs!(dhdu, u, idx_meas, idx_pos, idx_neg)
    n = length(idx_meas)

    # initialize to zero
    dhdu .= 0.0

    # now we update using direct measurements
    for i ∈ 1:n
        dhdu[i,idx_meas[i]] += 1.0
    end

    # now we update using positive ions
    for idx ∈ idx_pos
        dhdu[n+1, idx] += 1.0
    end

    # now we update using negative ions
    for idx ∈ idx_neg
        dhdu[n+1, idx] += 1.0
    end

end


function JObs(u, idx_meas, idx_pos, idx_neg)
    n = length(idx_meas)

    # initialize to zero
    dhdu = zeros(length(idx_meas)+2, length(u))

    # now we update using direct measurements
    for i ∈ 1:n
        dhdu[i,idx_meas[i]] += 1.0
    end

    # now we update using positive ions
    for idx ∈ idx_pos
        dhdu[n+1, idx] += 1.0
    end

    # now we update using negative ions
    for idx ∈ idx_neg
        dhdu[n+1, idx] += 1.0
    end

    return dhdu
end



function Rmat(idx_t::Int, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm((fudge_fac .* meas_ϵ[:,idx_t]) .^ 2)
end

function Rinv(idx_t::Int, meas_ϵ::AbstractArray{Float64}; fudge_fac::Float64=1.0)
    return diagm(1 ./ ((fudge_fac .* meas_ϵ[:,idx_t]) .^ 2))
end

