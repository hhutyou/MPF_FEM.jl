function Hn1_comp!(H0::Array{T},Hn1::Array{T,2},ψᵖ::Array{T},CC::Array{T}) where T<:Float64
    ##
    # σᴮ::Array{T,2} = σ
    # σᴮ::Array{T,2} = D0*(ε - εᵖ)
    # P, J2, s =FemBase.invariant(σᴮ,"plane-stress")
    # idg2=findall(dg2.<0.9)
    # ψᵉ = zeros(T,4*nel)
    # ψᵉ = gc/ls*d0*(1.0-(1.0-d0)^2)^2/2.0/(1.0-d0)*1.25e10 .+ zeros(4*nel)
    ##此处要保证Hn的初始值，不然第二次d循环求解时dg趋向于0，或者在初始化Hn的时候直接进行赋值
    # # β = zeros(T,4*nel)
    # # sigma_eq = zeros(T,4*nel)
    # σ₁, σ₃ = FemBase.principle(σᴮ)
    # # σc::Array{T,1}=((1.0.-dg2).^2 .+ k) .*(c .+h.*Eᵖ) ./(η/3.0-sqrt(2/3))
    # # β::Array{T,1} = (1-η) .*σc ./f0p .-(1+η)
    # function heaviside(x::Float64)
    #     if x>0
    #         x=1.0
    #     elseif x<=0
    #         x=0.0
    #     end
    #     return x
    # end
    # εᵉ_plus1 = operator_plus(εᵉ,"plane-strain")
    # ψᵉ= λ./2.0.*(heaviside.(operator_tr(εᵉ,"plane-strain"))).^2 .+
    #     G.*(εᵉ_plus1[1,:].^2 .+ εᵉ_plus1[2,:].^2 .+ 0.5*εᵉ_plus1[3,:].^2)
    # εᵉ_plus1 = operator_plus(εᵉ[:,Mat_ind0],"plane-strain")
    # ψᵉ[Mat_ind0].= λ[Mat_ind0]./2.0.*(heaviside.(operator_tr(εᵉ[:,Mat_ind0],"plane-strain"))).^2 .+
    #     G[Mat_ind0].*(εᵉ_plus1[1,:].^2 .+ εᵉ_plus1[2,:].^2 .+ 0.5*εᵉ_plus1[3,:].^2)
    # #
    # εᵉ_plus2 = operator_plus(εᵉ[:,Mat_ind12],"plane-strain")
    # ψᵉ[Mat_ind12] = λ2/2.0*(heaviside.(operator_tr(εᵉ[:,Mat_ind12],"plane-strain"))).^2 .+
    #     G2*(εᵉ_plus2[1,:].^2 .+ εᵉ_plus2[2,:].^2 .+ 0.5*εᵉ_plus2[3,:].^2)
    # ψᵉ[Mat_ind12] .= 1e19
    # #
    # ψᵉ = λ1/2.0*(heaviside.(operator_tr(εᵉ))).^2 .+
        # G1*((operator_dev(εᵉ)[1,:]).^2 .+ (operator_dev(εᵉ)[2,:]).^2 .+
        # 0.5.*operator_dev(εᵉ)[3,:].^2)
    # ψᵉ[Mat_ind1] = λ1/2.0*(heaviside.(operator_tr(εᵉ[:,Mat_ind1]))).^2 .+
    #    G1*((operator_dev(εᵉ[:,Mat_ind1])[1,:]).^2 .+ (operator_dev(εᵉ[:,Mat_ind1])[2,:]).^2 .+
    #    0.5.*operator_dev(εᵉ[:,Mat_ind1])[3,:].^2)
    # ψᵉ[Mat_ind2] = λ2/2.0*(heaviside.(operator_tr(εᵉ[:,Mat_ind2]))).^2 .+
    #    G2*((operator_dev(εᵉ[:,Mat_ind2])[1,:]).^2 .+ (operator_dev(εᵉ[:,Mat_ind2])[2,:]).^2 .+
    #    0.5.*operator_dev(εᵉ[:,Mat_ind2])[3,:].^2)
    D1::Array{T,1} = max.(0.0, H0 .+ ψᵖ./(gc1.*hc.(CC).+gc.*(1.0.-hc.(CC))))
    # D1 = zeros(size(ψᵉ))
    Hn1::Array{T,2} = max.(reshape(D1,4,nel),Hn1)
    return Hn1
end
##
function Hn1_comp!(H0::Array{T},Hn1::Array{T,2},ψᵖ::Array{T},CC::Array{T}) where T<:Float64
    D1::Array{T,1} = max.(0.0, H0 .+ ψᵖ./(hc.(CC).*gc1 .+ (1.0.-hc.(CC)).*gc))
    # D1 = zeros(size(ψᵉ))
    Hn1::Array{T,2} = max.(reshape(D1,4,nel),Hn1)
    return Hn1
end
#
function Hn2_comp!(ψᵖ::Array{T,1},Hn2::Array{T,2}) where T<:Float64
    D1::Array{T,1} = max.(0,ψᵖ)
    Hn2::Array{T,2} = max.(reshape(D1,4,size(element,1)),Hn2)
    return  Hn2
end
