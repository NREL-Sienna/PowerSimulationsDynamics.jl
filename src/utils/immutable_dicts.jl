function Base.ImmutableDict(KV::Pair{K,V}, KVs::Pair{K,V}...) where {K,V}
    d = Base.ImmutableDict(KV)
    for p in KVs
        d = Base.ImmutableDict(d,p)
    end
    return d
end
