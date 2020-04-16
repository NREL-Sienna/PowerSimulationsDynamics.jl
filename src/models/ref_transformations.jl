@enum dq_ref begin
    q = 1
    d = 2
end
@enum RI_ref begin
    R = 1
    I = 2
end
function dq_ri_gen(δ::Real)
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return [
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function ri_dq_gen(δ::Real)
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return [
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri_gen(δ::Float64)
    ## Uses the referenceframe of the Kundur page 852 of dq to RI
    return [
        sin(δ) cos(δ)
        -cos(δ) sin(δ)
    ]
end

function ri_dq_gen(δ::Float64)
    #Uses the reference frame of the Kundur page 852 of RI to dq
    return [
        sin(δ) -cos(δ)
        cos(δ) sin(δ)
    ]
end

function dq_ri_inv(δ::Real)
    ## Uses the referenceframe of the D'Arco paper of dq to RI
    return [
        cos(δ) -sin(δ)
        sin(δ) cos(δ)
    ]
end

function ri_dq_inv(δ::Real)
    ## Uses the referenceframe of the D'Arco paper of RI to dq
    return [
        cos(δ) sin(δ)
        -sin(δ) cos(δ)
    ]
end

function dq_ri_inv(δ::Float64)
    ## Uses the referenceframe of the D'Arco paper of dq to RI
    return [
        cos(δ) -sin(δ)
        sin(δ) cos(δ)
    ]
end

function ri_dq_inv(δ::Float64)
    ## Uses the referenceframe of the D'Arco paper of RI to dq
    return [
        cos(δ) sin(δ)
        -sin(δ) cos(δ)
    ]
end
