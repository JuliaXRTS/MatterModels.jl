@inline _internalize_density(ne::T) where {T <: Real} = ne
@inline function _internalize_density(ne::Quantity)
    return _density2me_dimless(ne)
end

@inline _internalize_temperature(temp::T) where {T <: Real} = temp
@inline function _internalize_temperature(temp::Quantity)
    return _energy2me_dimless(temp)
end
