function sigma_plus(εᵉ::Array{T,2},λ,G,planetype::String) where T<:Float64
    # s_plus = operator_plus(εᵉ[:,Mat_ind])
    # σ_plus = (kron( λ.*(heaviside.(operator_tr(εᵉ[:,Mat_ind])))', [1, 1, 0]) +
    #     2.0*G .* ([s_plus[1,:] s_plus[2,:] 0.5*s_plus[3,:]]' ))#+ kron( 1/2 .*(-operator_tr(εᵉ[:,Mat_ind])+abs.(-operator_tr(εᵉ[:,Mat_ind])))', [1/3, 1/3, 0])
    if planetype=="plane-stress"
        s_plus = operator_plus(εᵉ[:,Mat_ind],"plane-stress")
        σ_plus = λ.* kron((heaviside.(operator_tr(εᵉ[:,Mat_ind],"plane-stress")))', [1, 1, 0]) +
            2.0*G .* ([s_plus[1,:] s_plus[2,:] 0.5*s_plus[3,:]]' )
    elseif planetype=="plane-strain"
        s_plus = operator_plus(εᵉ,"plane-strain")
        σ_plus =  kron(λ', [1,1,1]).* kron((heaviside.(operator_tr(εᵉ,"plane-strain")))', [1, 1, 0]) +
            2.0*kron(G', [1,1,1]) .* ([s_plus[1,:] s_plus[2,:] 0.5*s_plus[3,:]]' )
    end
    return σ_plus
end
