function shapeFunc(node::Array{T2},element::Array{T1}) where {T1<:Int, T2<:Float64}
    # export Assemble_shapefunc,indices_fields
    ##1.feisoq
    function feisoq(s::T2,t::T2)
        #求导derivation
        dNds::Array{T2} = [-(1-t)/4,(1-t)/4,(1+t)/4,-(1+t)/4]
        dNdt::Array{T2} = [-(1-s)/4,-(1+s)/4,(1+s)/4,(1-s)/4]
        return dNds, dNdt
    end
    ##2.fejacob
    function fejacob(nnel::T1,dNds::Array{T2},dNdt::Array{T2},Xcoor::Array{T2},Ycoor::Array{T2})
        jacob=zeros(T2,2,2)
        for g=1:nnel
            jacob[1,1]=jacob[1,1]+dNds[g]*Xcoor[g]
            jacob[1,2]=jacob[1,2]+dNds[g]*Ycoor[g]

            jacob[2,1]=jacob[2,1]+dNdt[g]*Xcoor[g]
            jacob[2,2]=jacob[2,2]+dNdt[g]*Ycoor[g]
        end
        return jacob
    end
    ##3.federiv
    function federiv(nnel::T1,dNds::Array{T2},dNdt::Array{T2},invjacob::Array{T2})
        dNdx=zeros(T2,nnel)
        dNdy=zeros(T2,nnel)
        for s=1:nnel
            dNdx[s]=invjacob[1,1]*dNds[s]+invjacob[1,2]*dNdt[s]
            dNdy[s]=invjacob[2,1]*dNds[s]+invjacob[2,2]*dNdt[s]
        end
        return dNdx, dNdy
    end
    ##Assemble
    nel::T1=size(element,1); nnel::T1=size(element,2)
    Nd = zeros(T2,4,4*nel); Bd = zeros(T2,8,4*nel)
    Nu = zeros(T2,8,8*nel); Bu = zeros(T2,12,8*nel)
    Ndp = zeros(T2,16,4*nel); Bdp = zeros(T2,16,4*nel)
    detjacob=zeros(T2,4,nel)
    for iel=1:nel
        #节点坐标
        nd=zeros(T1,nnel)
        Xcoor=zeros(T2,nnel)
        Ycoor=zeros(T2,nnel)
        for i=1:nnel
            nd[i]=element[iel,i]
            Xcoor[i]=node[nd[i],1]
            Ycoor[i]=node[nd[i],2]
        end
        gauss=[-0.5774 -0.5774;-0.5774 0.5774;0.5774 -0.5774;0.5774 0.5774]
        for i=1:4
            s::T2=gauss[i,1]; t::T2=gauss[i,2]
            N1::T2=0.25*(1-s)*(1-t); N2::T2=0.25*(1+s)*(1-t)
            N3::T2=0.25*(1+s)*(1+t); N4::T2=0.25*(1-s)*(1+t)
            Nd[i,4*(iel-1)+1:4*iel]=[N1 N2 N3 N4]
            dNds,dNdt =feisoq(s,t)
            jacob =fejacob(nnel,dNds,dNdt,Xcoor,Ycoor) #雅可比矩阵
            detjacob[i,iel]=det(jacob)
            invjacob =inv(jacob)
            dNdx,dNdy =federiv(nnel,dNds,dNdt,invjacob)

            # Bd1 = [dNdx[1]; dNdy[1]]; Bd2 = [dNdx[2]; dNdy[2]]
            # Bd3 = [dNdx[3]; dNdy[3]]; Bd4 = [dNdx[4]; dNdy[4]]
            Bd[2*i-1:2*i,4*(iel-1)+1:4*iel] = [[dNdx[1]; dNdy[1]] [dNdx[2]; dNdy[2]] [dNdx[3]; dNdy[3]] [dNdx[4]; dNdy[4]]]

            # Nu1 = [N1 0; 0 N1]; Nu2 = [N2 0; 0 N2]
            # Nu3 = [N3 0; 0 N3]; Nu4 = [N4 0; 0 N4]
            Nu[2*i-1:2*i,8*(iel-1)+1:8*iel] = [[N1 0; 0 N1] [N2 0; 0 N2] [N3 0; 0 N3] [N4 0; 0 N4]]

            # Bu1 = [dNdx[1] 0.0; 0.0 dNdy[1]; dNdy[1] dNdx[1]]; Bu2 = [dNdx[2] 0; 0 dNdy[2]; dNdy[2] dNdx[2]]
            # Bu3 = [dNdx[3] 0.0; 0.0 dNdy[3]; dNdy[3] dNdx[3]]; Bu4 = [dNdx[4] 0; 0 dNdy[4]; dNdy[4] dNdx[4]]
            Bu[3*i-2:3*i,8*(iel-1)+1:8*iel] = [[dNdx[1] 0.0; 0.0 dNdy[1]; dNdy[1] dNdx[1]] [dNdx[2] 0; 0 dNdy[2]; dNdy[2] dNdx[2]] [
            dNdx[3] 0.0; 0.0 dNdy[3]; dNdy[3] dNdx[3]] [dNdx[4] 0; 0 dNdy[4]; dNdy[4] dNdx[4]]]
        end
        Ndp1 = Nd[1,4*(iel-1)+1:4*iel]*Nd[1,4*(iel-1)+1:4*iel]'
        Ndp2 = Nd[2,4*(iel-1)+1:4*iel]*Nd[2,4*(iel-1)+1:4*iel]'
        Ndp3 = Nd[3,4*(iel-1)+1:4*iel]*Nd[3,4*(iel-1)+1:4*iel]'
        Ndp4 = Nd[4,4*(iel-1)+1:4*iel]*Nd[4,4*(iel-1)+1:4*iel]'
        Ndp[:,4*(iel-1)+1:4*iel] = [Ndp1[:] Ndp2[:] Ndp3[:] Ndp4[:]]

        Bdp1 = Bd[1:2,4*(iel-1)+1:4*iel]'*Bd[1:2,4*(iel-1)+1:4*iel]
        Bdp2 = Bd[3:4,4*(iel-1)+1:4*iel]'*Bd[3:4,4*(iel-1)+1:4*iel]
        Bdp3 = Bd[5:6,4*(iel-1)+1:4*iel]'*Bd[5:6,4*(iel-1)+1:4*iel]
        Bdp4 = Bd[7:8,4*(iel-1)+1:4*iel]'*Bd[7:8,4*(iel-1)+1:4*iel]
        Bdp[:,4*(iel-1)+1:4*iel] = [Bdp1[:] Bdp2[:] Bdp3[:] Bdp4[:]]
    end
    return Nd,Ndp,Bdp,Bd,Bu,detjacob
end
#5.自由度指标
function indices_fields(element::Array{T}) where T<:Int
    nel = size(element,1)
    edofMat=zeros(T,nel,2*size(element,2))
    for iel = 1:nel  #loop within element
        edofMat[iel,1]=2*element[iel,1]-1
        edofMat[iel,2]=2*element[iel,1]
        edofMat[iel,3]=2*element[iel,2]-1
        edofMat[iel,4]=2*element[iel,2]
        edofMat[iel,5]=2*element[iel,3]-1
        edofMat[iel,6]=2*element[iel,3]
        edofMat[iel,7]=2*element[iel,4]-1
        edofMat[iel,8]=2*element[iel,4]
    end
    iK::Array{T} = reshape(kron(edofMat,ones(8,1))',64*nel)
    jK::Array{T} = reshape(kron(edofMat,ones(1,8))',64*nel)
    #d的自由度
    # dedofMat = element
    iKd::Array{T} = reshape(kron(element,ones(4,1))',16*nel)
    jKd::Array{T} = reshape(kron(element,ones(1,4))',16*nel)
    iFd::Array{T} = reshape(element',4*nel)
    jFd::Array{T} = ones(4*nel)
    return iK, jK, iKd, jKd, iFd, jFd, edofMat, element
end
@info "Formulating the shape functions and indices takes"
@time begin
    Nd,Ndp,Bdp,Bd,Bu,detjacob = shapeFunc(node,element)
    iK, jK, iKd, jKd, iFd, jFd, edofMat, dedofMat = indices_fields(element)
    GC.gc()
end
# @time iK, jK, iKd, jKd, iFd, jFd, edofMat, dedofMat = indices_fields(element)
