


function Obs!(h, u, idx_meas, idx_pos, idx_neg)
    n = length(idx_meas)
    h[1:n] = u[idx_meas]
    h[n+1] = sum(u[idx_pos])
    h[n+2] = sum(u[idx_neg])
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

