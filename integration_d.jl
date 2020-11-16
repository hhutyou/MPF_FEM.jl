function integration_d(Hn1::Array{T},d1::Array{T},dg1::Array{T}) where T<:Float64 ##
    #
    function comput_r(d,dg,Hn1)
        dgdd = -2.0.*(1.0 .- reshape(dg,4,nel))./ka ./(1.0 .- (1.0 .- reshape(dg,4,nel)).^2).^2
        Hn = Hn1.*dgdd
        Nterm = kron(1.0/ls,ones(4,nel))
        Bterm = kron(1.0*ls,ones(4,nel))
        sKKd_B=SharedArray{Float64,2}(16,nel)
        sKKd_D=SharedArray{Float64,2}(16,nel)
        sFD=SharedArray{Float64,2}(4,nel)
        @sync @distributed for iel=1:nel
            sKKd_B[:,iel]=Bdp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*Bterm[:,iel]
            sKKd_D[:,iel]=Ndp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*Nterm[:,iel]
            sFD[:,iel]=Nd[:,4*(iel-1)+1:4*iel]'.*kron(detjacob[:,iel]',ones(4,1))*Hn[:,iel]
        end
        # sKd_B::Array{Float64,1} = reshape(sKKd_B,16*nel)
        # sKd_D::Array{Float64,1} = reshape(sKKd_D,16*nel)
        # sKd = sKd_B+sKd_D
        # sFd = reshape(sFD,4*nel)
        Kd = sparse(iKd,jKd,reshape(sKKd_B,16*nel)+reshape(sKKd_D,16*nel))#::SparseMatrixCSC
        # Kd::SparseMatrixCSC = sparse(iKd,jKd,reshape(sKKd_D,16*nel)) ::SparseVector
        Fd = sparse(iFd,jFd,reshape(sFD,4*nel))#::SparseVector
        r = Kd * d .+ Fd
        return Array(r)[:]
    end
    function comput_Kr(d,dg,Hn1)
        # ddgdd = (2.0/ka) ./(1.0 .- (1.0 .- reshape(dg,4,nel)).^2).^2
        ddgdd = (6.0.*(1.0 .- reshape(dg,4,nel)).^2 .+ 2.0) ./ka ./(1.0 .- (1.0 .- reshape(dg,4,nel)).^2).^3
        # Hn = Hn1./ka ./(1.0 .- (1.0 .- reshape(dg,4,nel)).^2).^2
        Nterm = kron(1.0/ls,ones(4,nel)) .+ ddgdd.*Hn1
        Bterm = kron(1.0*ls,ones(4,nel))
        sKKd_B=SharedArray{Float64,2}(16,nel)
        sKKd_D=SharedArray{Float64,2}(16,nel)
        @sync @distributed for iel=1:nel
            sKKd_B[:,iel]=Bdp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*Bterm[:,iel]
            sKKd_D[:,iel]=Ndp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*Nterm[:,iel]
        end
        # sKd_B::Array{Float64,1} = reshape(sKKd_B,16*nel)
        # sKd_D::Array{Float64,1} = reshape(sKKd_D,16*nel)
        # sKd = sKd_B+sKd_D
        # sFd = reshape(sFD,4*nel)
        Kr = sparse(iKd,jKd,reshape(sKKd_B,16*nel)+reshape(sKKd_D,16*nel))#::SparseMatrixCSC
        # Kd::SparseMatrixCSC = sparse(iKd,jKd,reshape(sKKd_D,16*nel)) ::SparseVector
        return Kr
    end
    ## 牛顿法迭代求d
    d1_old = deepcopy(d1)
    err_d = 1.0
    nit_d = 0
    while (err_d>1e-6) && (nit_d<30)
        nit_d += 1
        # @info "-d-迭代步nit_d=: $nit_d"
        # Hn = ka.*Hn1./(ka .- (ka.-1.0).*(1.0 .- reshape(dg1,4,nel)).^2).^2
        # Hn = Hn1./ka ./(1.0 .- (1.0 .- reshape(dg1,4,nel)).^2).^2
        # Hn = copy(Hn1)
        d1 = d1 .- comput_Kr(d1,dg1,Hn1)\comput_r(d1,dg1,Hn1)
        # replace!(x->x<d0 ? d0 : x, d1)
        err_d = maximum(abs.(d1.-d1_old))
        d1_old .= d1
        dg1=zeros(Float64,4*nel)
        for iel = 1:nel
            dg1[4*(iel-1)+1:4*iel]=Nd[:,4*(iel-1)+1:4*iel]*d1[dedofMat[iel,:]]
        end
        # @info "Hn1=$(maximum(Hn1)),d1=$(maximum(d1)),err_d=$err_d"
        # dg1 = sdata(dg1_t)
        # dg2 = sdata(dg2_t)
    end
    # replace!(x->x<d0 ? d0 : x, d1)
    # # replace!(x->x>1.0 ? 1.0 : x, d1)
    # dg1=zeros(Float64,4*nel)
    # for iel = 1:nel
    #     dg1[4*(iel-1)+1:4*iel]=Nd[:,4*(iel-1)+1:4*iel]*d1[dedofMat[iel,:]]
    # end
    return d1, dg1
end
