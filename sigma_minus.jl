function sigma_minus(εᵉ::Array{T,2},Mat_ind,λ,G,planetype::String) where T<:Float64
    # s_minus = operator_minus(εᵉ[:,Mat_ind])
    # σ_minus = (kron( λ.*(mheaviside.(operator_tr(εᵉ[:,Mat_ind])))', [1, 1, 0]) +
    #     2.0*G .* ([s_minus[1,:] s_minus[2,:] 0.5*s_minus[3,:]]' ))#+ kron( 1/2 .*(-operator_tr(εᵉ[:,Mat_ind])+abs.(-operator_tr(εᵉ[:,Mat_ind])))', [1/3, 1/3, 0])
    if planetype=="plane-stress"
        s_minus = operator_minus(εᵉ[:,Mat_ind],"plane-stress")
        σ_minus = λ.* kron( (mheaviside.(operator_tr(εᵉ[:,Mat_ind],"plane-stress")))', [1, 1, 0]) +
            2.0*G .* ([s_minus[1,:] s_minus[2,:] 0.5*s_minus[3,:]]' )
    elseif planetype=="plane-strain"
        s_minus = operator_minus(εᵉ[:,Mat_ind],"plane-strain")
        σ_minus = λ.* kron( (mheaviside.(operator_tr(εᵉ[:,Mat_ind],"plane-strain")))', [1, 1, 0]) +
            2.0*G .* ([s_minus[1,:] s_minus[2,:] 0.5*s_minus[3,:]]' )
    end
    return σ_minus
end
