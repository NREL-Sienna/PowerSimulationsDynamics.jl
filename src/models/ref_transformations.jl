function dq_ri(δ)
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return [
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function ri_dq(δ)
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return [
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end
