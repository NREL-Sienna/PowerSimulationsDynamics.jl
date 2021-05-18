function dq_ri(δ::T) where {T <: Number}
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return T[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function ri_dq(δ::T) where {T <: Number}
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return T[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri(δ::Float64)
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return Float64[
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function ri_dq(δ::Float64)
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return Float64[
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end
