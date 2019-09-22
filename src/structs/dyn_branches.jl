mutable struct DynLine <: DynBranch
    name::Symbol
    available::Bool
    arc::PSY.Arc
    r::Float64  # System per-unit value
    x::Float64  # System per-unit value
    b::NamedTuple{(:from, :to), Tuple{Float64, Float64}}  # System per-unit value
    n_states::Int64
    states::Vector{Symbol}


                function DynLine(name::String,
                                 available::Bool,
                                 arc::PSY.Arc,
                                 r::Float64, # System per-unit value
                                 x::Float64, # System per-unit value
                                 b::NamedTuple{(:from, :to), Tuple{Float64, Float64}})

                    n_states = 2
                    states = Vector{Symbol}()
                    states = [:Il_R
                              :Il_I]

                    new(Symbol(name),
                        available,
                        arc,
                        r,
                        x,
                        b,
                        n_states,
                        states)

                end

end
