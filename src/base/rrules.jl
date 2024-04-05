#RRules for ChainRulesCore
#= function ChainRulesCore.rrule(::typeof(_filter_kwargs), kwargs)
    #CRC.@ignore_derivatives @error "HIT"
    #@assert false 
    project_A = ProjectTo(A)
    y = _filter_kwargs(kwargs)
    function _filter_kwargs_pullback(ȳ)
        f̄ = NoTangent()
        a = NoTangent() 
        return f̄, project_A(a) 
    end
    return y, _filter_kwargs_pullback
end =#
#= function rrule(::typeof(*), A::AbstractMatrix, B::AbstractMatrix)
    project_A = ProjectTo(A)
    project_B = ProjectTo(B)
    function times_pullback(ȳ)
        dA = ȳ * B'
        dB = A' * ȳ
        return NoTangent(), project_A(dA), project_B(dB)
    end
    return A * B, times_pullback
end
 =#

#= function ChainRulesCore.rrule(::typeof(_filter_kwargs_simple), kwargs)
    #CRC.@ignore_derivatives @error "HIT"
    #@assert false 
    y = _filter_kwargs_simple(kwargs)
    function _filter_kwargs_pullback(ȳ)
        f̄ = NoTangent()
        a = NoTangent() 
        return f̄, a 
    end
    return y, _filter_kwargs_pullback
end =#

function _simple_f(x)
    return 2.0 * x
end

function ChainRulesCore.rrule(::typeof(_simple_f), x)
    y = _simple_f(x)
    function _simple_f_pullback(ȳ)
        f̄ = ChainRulesCore.NoTangent()
        a = 2.0 * ȳ
        return f̄, a
    end
    return y, _simple_f_pullback
end
