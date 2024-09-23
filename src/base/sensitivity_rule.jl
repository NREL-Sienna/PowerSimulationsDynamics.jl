function solve_up_with_callback(
    prob::Union{DiffEqBase.AbstractDEProblem, DiffEqBase.NonlinearProblem}, sensealg, u0, p,
    cb,
    args...; kwargs...)
    alg = DiffEqBase.extract_alg(
        args,
        kwargs,
        DiffEqBase.has_kwargs(prob) ? prob.kwargs : kwargs,
    )
    if isnothing(alg) || !(alg isa DiffEqBase.AbstractDEAlgorithm) # Default algorithm handling
        _prob = DiffEqBase.get_concrete_problem(prob, !(prob isa DiscreteProblem); u0 = u0,
            p = p, callback = cb, kwargs...)
        solve_call(_prob, args...; callback = cb, kwargs...)
    else
        _prob = DiffEqBase.get_concrete_problem(
            prob,
            DiffEqBase.isadaptive(alg);
            u0 = u0,
            p = p,
            callback = cb,
            kwargs...,
        )
        _alg = DiffEqBase.prepare_alg(alg, _prob.u0, _prob.p, _prob)
        DiffEqBase.check_prob_alg_pairing(_prob, alg) # use alg for improved inference
        if length(args) > 1
            DiffEqBase.solve_call(_prob, _alg, Base.tail(args)...; callback = cb, kwargs...)
        else
            DiffEqBase.solve_call(_prob, _alg; callback = cb, kwargs...)
        end
    end
end

function solve_with_callback(prob::DiffEqBase.AbstractDEProblem, cb, args...;
    sensealg = nothing, #add cb as second argument. 
    u0 = nothing, p = nothing, wrap = Val(true), kwargs...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    if wrap isa Val{true}
        DiffEqBase.wrap_sol(
            solve_up_with_callback(prob, sensealg, u0, p, cb, args...; kwargs...),
        )
    else
        solve_up_with_callback(prob, sensealg, u0, p, cb, args...; kwargs...)
    end
end

function Enzyme.EnzymeRules.augmented_primal(config::Enzyme.EnzymeRules.ConfigWidth{1},
    func::Enzyme.Const{typeof(solve_up_with_callback)},
    ::Type{<:Union{Enzyme.Duplicated{RT}, Enzyme.DuplicatedNoNeed{RT}}}, prob,
    sensealg::Union{
        Enzyme.Const{Nothing},
        Enzyme.Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
    },
    u0, p, cb, args...; kwargs...) where {RT}
    @inline function copy_or_reuse(val, idx)
        if Enzyme.EnzymeRules.overwritten(config)[idx] && ismutable(val)
            return deepcopy(val)
        else
            return val
        end
    end

    @inline function arg_copy(i)
        copy_or_reuse(args[i].val, i + 5)
    end

    res = DiffEqBase._solve_adjoint(
        copy_or_reuse(prob.val, 2), copy_or_reuse(sensealg.val, 3),
        copy_or_reuse(u0.val, 4), copy_or_reuse(p.val, 5),
        SciMLBase.EnzymeOriginator(), ntuple(arg_copy, Val(length(args)))...;
        callback = copy_or_reuse(cb.val, 6),
        kwargs...)

    dres = Enzyme.make_zero(res[1])::RT
    tup = (dres, res[2])
    return Enzyme.EnzymeRules.AugmentedReturn{RT, RT, Any}(res[1], dres, tup::Any)
end

function Enzyme.EnzymeRules.reverse(config::Enzyme.EnzymeRules.ConfigWidth{1},
    func::Enzyme.Const{typeof(solve_up_with_callback)},
    ::Type{<:Union{Enzyme.Duplicated{RT}, Enzyme.DuplicatedNoNeed{RT}}}, tape,
    prob,
    sensealg::Union{
        Enzyme.Const{Nothing},
        Enzyme.Const{<:DiffEqBase.AbstractSensitivityAlgorithm},
    },
    u0, p, cb, args...; kwargs...) where {RT}
    dres, clos = tape
    dres = dres::RT
    dargs = clos(dres)
    for (darg, ptr) in zip(dargs, (func, prob, sensealg, u0, p, cb, args...))
        if ptr isa Enzyme.Const
            continue
        end
        if darg == ChainRulesCore.NoTangent()
            continue
        end
        ptr.dval .+= darg
    end
    Enzyme.make_zero!(dres.u)
    return ntuple(_ -> nothing, Val(length(args) + 5))
end
