
function LSFunc(LS::T,δu::Array{T},u_old::Array{T},εᵖ::Array{T},d1::Array{T},dg1::Array{T},H0::Array{T},Hn1::Array{T},εᵖ_com::Array{T},f_ext::Array{T},r_ref::Array{T}) where T<:Float64
    m0(x) = (1.0.-x).^2 ./ (ka .- (ka .- 1.0).*(1.0 .- x).^2)
    gd(x) = (1.0.-x).^2 ./ (ka .- ka.*(1.0 .- x).^2)
    σ = zeros(Float64,3,4*nel)
    εᵉ = zeros(Float64,3,4*nel)
    δεᵖ = zeros(Float64,3,4*nel)
    CC = zeros(Float64,4*nel)
    σᵗ = SharedArray{Float64,2}(12,nel)
    σcᵗ = SharedArray{Float64,2}(12,nel)
    σd_t = SharedArray{Float64,2}(12,nel)
    σc =zeros(Float64,3,4*nel)
    σd =zeros(Float64,3,4*nel)
    u = u_old .+ LS * δu
    # du = du_old + δu
    U = u[edofMat]
    ε = zeros(Float64,3,4*nel)
    for iel=1:nel
        ε[:,4*(iel-1)+1:4*iel]=reshape(Bu[:,8*(iel-1)+1:8*iel]*U[iel,:],3,4)
    end
    # εᵉ .= ε .- εᵖ
    @sync @distributed for iel=1:nel
        σᵗ[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
        reshape(ε[:,4*(iel-1)+1:4*iel] .- εᵖ[:,4*(iel-1)+1:4*iel],12,1)
        σcᵗ[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
        reshape(εᵖ[:,4*(iel-1)+1:4*iel],12,1)
        # σ[:,4*(iel-1)+1:4*iel] = reshape(σᵗ[:,iel],3,:) .+ kron((1.0 .-dg1[4*(iel-1)+1:4*iel]).^2 .- 1.0, ones(1,3))'.*
        #     sigma_plus(εᵉ, collect(4*(iel-1)+1:4*iel), λ[4*(iel-1)+1:4*iel],G[4*(iel-1)+1:4*iel], "plane-strain")
    end
    σ .= reshape(sdata(σᵗ),3,4*nel)
    ##计算局部应力
    σc .= σ .- reshape(kron(gd.(dg1),ones(3)) .* sdata(σcᵗ[:]),3,4*nel)
    CC .= sigma_P(σc,v,planetype)
    # σ1, maxfL = principle(σ,"plane-strain")
    idxTen=findall(CC.>=0.0)##拉伸断裂
    idxCom=findall(CC.<0.0)##压剪断裂
    if isempty(idxTen)==0  ##拉
        εᵖ[:,idxTen] .= kron(1.0 .- m0.(dg1[idxTen]), ones(1,3))' .* (ε[:,idxTen])
    end
    δεᵖ .= plasticity!(σc,idxCom)
    εᵖ[:,idxCom] .= εᵖ[:,idxCom] .+ δεᵖ[:,idxCom]
    εᵖ_com .= εᵖ_com .+ δεᵖ
    ##
    @sync @distributed for iel=1:nel
        σᵗ[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
        reshape(ε[:,4*(iel-1)+1:4*iel] .- εᵖ[:,4*(iel-1)+1:4*iel],12,1)
    end
    σ .= reshape(sdata(σᵗ),3,4*nel)
    f_int = FMat(node,element,Bu,detjacob,σ,edofMat)
    # S = (f_ext-f_int)'*δu
    S = norm(f_ext[freedofs] .- f_int[freedofs])/norm(r_ref)
    return S, u, σ, ε, f_int, d1, dg1, Hn1, εᵖ, εᵖ_com, CC
end
