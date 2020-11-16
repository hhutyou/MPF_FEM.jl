function plasticity!(σc::Array{T},sigma::Array{T},εᵖ::Array{T},dg1::Array{T},idxCom) where T<:Float64
    Ip=[1.0,1.0,0.0]
    itol=1e-6
    imaxit=100
    # m0(x) = (1.0.-x).^2 ./ (ka .- (ka .- 1.0).*(1.0 .- x).^2)
    gd(x) = (1.0.-x).^2 ./ (ka .- ka.*(1.0 .- x).^2)
    # Eᵖ_old = deepcopy(Eᵖ)
    # nel::Int=size(element,1)
    F = Array{Float64}(undef,4*nel)
    δλ = zeros(Float64,4*nel)
    δεᵖ = zeros(Float64,3,4*nel)
    ep = operator_dev(εᵖ,planetype)
    evp = operator_tr(εᵖ,planetype)
    P, J2, s = invariant(σc,v,planetype)
    Pm = sigma_P(sigma,v,planetype)
    F .= sqrt.(2.0*J2).+ AA.*P
    idxYield = findall(F.>itol)
    # UidxCY = intersect(idxYield)
    UidxCY = intersect(idxYield,idxCom)
    if isempty(UidxCY)
        @info "----未屈服"
        return δεᵖ
    else
        ##并行计算数组初始化
        # @info "----塑性调整"
        if !isempty(UidxCY)
            δλC::SharedArray = δλ[UidxCY]
            FC::SharedArray = F[UidxCY]
            J2C::SharedArray = J2[UidxCY]
            PC::SharedArray = P[UidxCY]
            sC::SharedArray = s[:,UidxCY]
            PmC::SharedArray = Pm[UidxCY]
            # J2mC::SharedArray = J2m[UidxCY]
            σC::SharedArray = σc[:,UidxCY]
            # σᵗC::SharedArray =  copy(σC)
            sigmaC::SharedArray = sigma[:,UidxCY]
            sigmaᵗC::SharedArray =  copy(sigmaC)
            δεᵖC::SharedArray = δεᵖ[:,UidxCY]
            εᵖC::SharedArray = εᵖ[:,UidxCY]
            epC = ep[:,UidxCY]
            evpC = evp[UidxCY]
            dg1C = dg1[UidxCY]
            GC = G[UidxCY]
            KvC = Kv[UidxCY]
            AC = AA[UidxCY]
            DKC = DK[:,UidxCY]
            sprime = sigma_s([σC;v*(σC[1,:]+σC[2,:])'],v,planetype)
            @sync @distributed for num = 1:size(UidxCY,1)#@distributed (+)
                δλC[num] = 0.0
                ierr=1.0
                init=0
                while ierr>itol && init<imaxit
                    init += 1
                    # @info ("---塑性迭代init=$init")
                    d = -2.0*(1.0+gd(dg1C[num]))*GC[num]-(1.0+gd(dg1C[num]))*KvC[num]*(AC[num])^2
                    ##残余屈服值对塑性因子的导数%%参考computation methods for plasticity-P332
                    δλC[num]=δλC[num]-FC[num]/d
                    ##检查收敛性
                    sigmaC[:,num] = sigmaᵗC[:,num]-δλC[num].*(2.0*GC[num].*sC[:,num]/sqrt(2.0*J2C[num]).+KvC[num]*AC[num]*Ip)
                    sigmaZ = v*(sigmaᵗC[1,num]+sigmaᵗC[2,num])-δλC[num]*(2.0*GC[num]*sprime[4,num]/sqrt(2.0*J2C[num])+KvC[num]*AC[num])
                    σC[:,num] = sigmaC[:,num] .- gd(dg1C[num]).*(2.0*GC[num]*(epC[:,num].+δλC[num].*sC[:,num]/sqrt(2.0*J2C[num]))+KvC[num]*(evpC[num]+δλC[num]*AC[num])*Ip)
                    sigmacZ = sigmaZ - gd(dg1C[num])*(2.0*GC[num]*(-evpC[num]/3.0+δλC[num]*sprime[4,num]/sqrt(2.0*J2C[num]))+KvC[num]*(evpC[num]+δλC[num]*AC[num]))
                    ss = sigma_s([σC[:,num];sigmacZ],v,planetype)
                    FC[num]=sqrt(sum([ss;ss[3]]'*[ss;ss[3]]))+AC[num]*(PmC[num]-(1.0+gd(dg1C[num]))*KvC[num]*AC[num]*δλC[num]-gd(dg1C[num])*KvC[num]*evpC[num])
                    ierr = abs(FC[num])
                end
                δεᵖC[:,num] = δλC[num]*(sC[:,num]/sqrt(2.0*J2C[num])+AC[num]/3.0*Ip)
            end
            # σc[:,UidxCY] .=   σC
            # sigma[:,UidxCY] .= sdata(sigmaC)
            δεᵖ[1:2,UidxCY] = sdata(δεᵖC)[1:2,:]
            δεᵖ[3,UidxCY] = 2.0*sdata(δεᵖC)[3,:]
        end
    end
    return δεᵖ
end
