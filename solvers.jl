function solvers(u_inc::T, conf::T) where T<:Float64
    ## Initialization
    m0(x) = (1.0.-x).^2 ./ (ka .- (ka .- 1.0).*(1.0 .- x).^2)
    gd(x) = (1.0.-x).^2 ./ (ka .- ka.*(1.0 .- x).^2)
    u = zeros(Float64,2*nnode)
    u_old = zeros(Float64,2*nnode)
    f_ext=zeros(Float64,2*nnode)
    f_int=zeros(Float64,2*nnode)
    f_int_old=zeros(Float64,2*nnode)
    sigma=zeros(Float64,3,4*nel)
    σc =zeros(Float64,3,4*nel)
    CC = zeros(Float64,4*nel)
    σd =zeros(Float64,3,4*nel)
    σ_t = SharedArray{Float64,2}(12,nel) ## parallel computation of stress
    σc_t = SharedArray{Float64,2}(12,nel)
    σd_t = SharedArray{Float64,2}(3,4*nel)
    epsilon = zeros(Float64,3,4*nel)
    ψᵖ=zeros(Float64,4*nel)
    εᵖ = zeros(Float64,3,4*nel)
    εᵖ_old = zeros(Float64,3,4*nel)
    δεᵖ = zeros(Float64,3,4*nel)
    εᵖ_com = zeros(Float64,3,4*nel)
    εᵖ_com_old = zeros(Float64,3,4*nel)
    Hn1 = zeros(Float64,4,nel)
    Hn1_old=copy(Hn1)
    d1 = zeros(Float64,nnode).+d0; dg1=zeros(Float64,4*nel)
    for iel = 1:nel
        dg1[4*(iel-1)+1:4*iel]=Nd[:,4*(iel-1)+1:4*iel]*d1[dedofMat[iel,:]]
    end
    d1_old = deepcopy(d1)
    dg1_old = deepcopy(dg1)
    H0 = Initial_Hn(d1, dg1)
    Hn1 = Hn1_comp!(H0,Hn1,ψᵖ,zeros(Float64,4*nel)) ## prescribed initial history field
    err_storage = zeros(Float64,maxit,step_total+1) ##残差和收敛时间
    time_storge = zeros(Float64,step_total+1) ##收敛时间
    #Assemble whole Stiffness matrix
    begin
        @info "Formulating stiffness matrix takes"
        @time KK=Kmatrix(element, Bu, detjacob,DK,iK,jK)
    end
    init_total = 1 ## step number in the first load inc
    f_ext = FMat_ext(node,conf) ## external force
    alldofs = 1:2size(node,1)
    ###
    u_increment = u_inc
    STOnit=zeros(Float64,step_total)
    function monoli_initialStiff()
        for stp=1:step_total
            @info "step=: $stp"
            r_ref = zeros(Float64,2*nnode)
            u_ref = zeros(Float64,2*nnode)
            du = zeros(Float64,2*nnode)
            du[loaddofs] .= u_increment
            ## staggered iteration of u and d
            err_d = 1.0; nit_d = 0
            record = @timed while (err_d>1e-3) && (nit_d<maxit)
                nit_d+=1
                @info "nit_d=$nit_d"
                # compute displacement field u
                err=1; nit=0
                δu = zeros(Float64,2*nnode)
                rr = zeros(Float64,2*nnode)
                while (err>tol) && (nit<maxit)
                    nit+=1
                    itd = 0
                    @info "-迭代步=: $nit"
                    if nit_d==1 && nit==1
                        rr[freedofs] .= KK[freedofs,loaddofs]*du[loaddofs] # .- (f_ext[freedofs]-f_int_old[freedofs])
                        r_ref .= rr
                        δu[freedofs] .= -KK[freedofs,freedofs]\(rr[freedofs])
                        u .= u_old  .+ du .+ δu
                        u_ref .=  du .+ δu
                    elseif nit>1
                        δu[freedofs] .= KK[freedofs,freedofs]\(rr[freedofs])
                        u .= u_old .+ δu
                    else
                        itd = 1
                        nothing
                    end
                    UU = u[edofMat]
                    # epsilon = zeros(Float64,3,4*nel)
                    for iel=1:nel
                        epsilon[:,4*(iel-1)+1:4*iel]=reshape(Bu[:,8*(iel-1)+1:8*iel]*UU[iel,:],3,4)
                    end
                    @sync @distributed for iel=1:nel
                        σ_t[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
                        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
                        reshape((epsilon[:,4*(iel-1)+1:4*iel] .- εᵖ[:,4*(iel-1)+1:4*iel]),12,1)
                        σc_t[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
                        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
                        reshape(εᵖ[:,4*(iel-1)+1:4*iel],12,1)
                    end
                    sigma .=  reshape(sdata(σ_t),3,4*nel)
                    ##local stress
                    σc .= sigma .- reshape(kron(gd.(dg1),ones(3)) .* sdata(σc_t[:]),3,4*nel)
                    ###
                    CC .= sigma_P(σc,v,planetype) ## mean stress / microcrack open-closure criterion
                    idxTen = findall(CC.>=0.0)##拉伸断裂
                    idxCom = findall(CC.<0.0)##压剪断裂
                    if isempty(idxTen)==0  ##拉
                        # @info "拉伸"
                        εᵖ[:,idxTen] .= kron(1.0 .- m0.(dg1[idxTen]), ones(1,3))' .* (epsilon[:,idxTen])
                    end
                    # increment of plastic strian
                    δεᵖ .= plasticity!(σc,sigma,εᵖ,dg1,idxCom)
                    εᵖ[:,idxCom] .= εᵖ_old[:,idxCom] .+ δεᵖ[:,idxCom]
                    εᵖ_com .= εᵖ_com_old .+ δεᵖ # plastic strain due to friction
                    @sync @distributed for iel=1:nel
                        σ_t[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
                        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
                        reshape(epsilon[:,4*(iel-1)+1:4*iel] .- εᵖ[:,4*(iel-1)+1:4*iel],12,1)
                    end
                    sigma .= reshape(sdata(σ_t),3,4*nel)
                    ##
                    f_int = FMat(node,element,Bu,detjacob,sigma,edofMat)
                    # line search: secant line search
                    if abs((f_ext.-f_int)'*δu) > (0.8 * abs((f_ext.-f_int_old)' * δu)) && nit > 6
                        # line Search
                        # @info "-line Search"
                        nit_int = 0
                        LS = 0.0
                        LSS = 1.0
                        LSSS = 1.2
                        S = []
                        SS = []
                        SSS = []
                        id::Int64 = 0
                        err_ls=1
                        while err_ls > tol && nit_int < 15
                            nit_int += 1
                            LS =  LSS
                            LSS = LSSS
                            #
                            SS = LSFunc(LSS,δu,u_old,εᵖ_old,d1_old,dg1_old,H0,Hn1_old,εᵖ_com_old,f_ext,r_ref)
                            S = LSFunc(LS,δu,u_old,εᵖ_old,d1_old,dg1_old,H0,Hn1_old,εᵖ_com_old,f_ext,r_ref)
                            LSSS = LSS - SS[1] / (SS[1] - S[1]) * (LSS - LS)
                            LSSS>30 || LSSS<-10 || abs(LSSS-LSS)<tol  ? (@info "break" break) : nothing
                            LSSS>15.0 ? LSSS=15.0 : LSSS<0.1 ? LSSS=0.1 : nothing
                            ##计算line search目标函数
                            # f_int .= SS[5]
                            err_ls =  SS[1]/norm(r_ref)
                            @info "-line search err=$err_ls"
                        end
                        u .= SS[2]
                        sigma .= SS[3]
                        epsilon .= SS[4]
                        f_int .= SS[5]
                        d1 .= SS[6]
                        dg1 .= SS[7]
                        Hn1 .= SS[8]
                        εᵖ .= SS[9]
                        εᵖ_com .= SS[10]
                        CC .= SS[11]
                    end
                    rr[freedofs] .= f_ext[freedofs] .- f_int[freedofs]
                    # if itd == 1
                    #     err = 1.0
                    # else
                    #     err = abs((f_ext .- f_int)'*δu)/abs(r_ref' *u_ref)
                    # end
                    err = norm(rr[freedofs])/norm(r_ref)
                    @info "--err_u=$err"
                    f_int_old .= f_int
                    u_old .= u
                    d1_old .= d1
                    dg1_old .= dg1
                    Hn1_old .= Hn1
                    εᵖ_old .= εᵖ
                    εᵖ_com_old .= εᵖ_com
                end
                # compute damage field d
                @sync @distributed for iel=1:nel
                    σc_t[:,iel] = blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
                    sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*
                    reshape(εᵖ[:,4*(iel-1)+1:4*iel],12,1)
                end
                σd .= reshape(sdata(σc_t),3,4*nel)
                ψᵖ .= 1.0/2.0.*(σd[1,:].*εᵖ[1,:] .+ σd[2,:].*εᵖ[2,:] .+ σd[3,:].*εᵖ[3,:])
                ##计算d
                Hn1 = Hn1_comp!(H0,Hn1,ψᵖ,CC)
                d1, dg1 = integration_d(Hn1,d1,dg1)
                # replace!(x->x>1.0 ? (1.0-k) : x, dg1)
                # replace!(x->x<0.0 ? k : x, dg1)
                err_d = norm(d1.-d1_old)/norm(d1)
                # d1_old .= d1
                # dg1_old .= dg1
                @info "---err_d=$err_d"
                err_storage[nit_d,stp] = err_d
            end
            # "Storing the output data takes:"
            begin
                Fload1[stp+1] = sum(f_int[loaddofs])/(size(loaddofs,1)-1)
                Uload1[stp+1] = sum(u[loaddofs])/size(u[loaddofs],1)
                time_storge[stp+1] = record[2]
                if mod(stp,aa) == 0
                    n = Int(stp/aa)
                    numD[:,n].=d1[:]
                    numD2[:,:,n].=εᵖ_com
                    numD3[:,n].=u
                    numD4[:,:,n].=sigma
                    numD5[:,:,n].=εᵖ
                    numD6[:,:,n].= σc
                    numD7[:,:,n] .= Hn1
                    numD8[:,n] .= err_storage[:,n]
                end
            end
            # STOnit[stp]=nit
            # if nit==maxit && STOnit[stp-1] != maxit
            #     u_increment=0.25*u_increment
            # elseif  STOnit[stp]==maxit && STOnit[stp-1]==maxit && STOnit[stp-2]==maxit
            #     @info "Cannot converge"
            #     break
            # end
        end
        return Fload1, Uload1, time_storge, numD, numD2, numD3, numD4, numD5, numD6, numD7, numD8
    end
    return monoli_initialStiff()
end
