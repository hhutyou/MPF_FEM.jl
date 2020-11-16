function sigma_dev(εᵉ::Array{T,2},dg1::Array{T},Mat_ind,λ,G) where T<:Float64
    σ_dev = (kron( λ .*(heaviside.(operator_tr(εᵉ)))', [1, 1, 0]) +
        2.0*G .* ([εᵉ[1,:] εᵉ[2,:] 0.5*εᵉ[3,:]]' .+ kron((heaviside.(-operator_tr(εᵉ)))', [1/3, 1/3, 0])))#
    return σ_dev
end
