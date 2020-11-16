# include("node_element.jl")
# include("Interpfunc.jl")
# Definition of global stiffness_matrix
#
function Kmatrix(element::Array{T1}, Bu::Array{T2}, detjacob::Array{T2},DK::Array{T2},iK::Array{T1},jK::Array{T1}) where {T1<:Int, T2<:Float64}
    sK=SharedArray{Float64,2}(64,size(element,1))
    @sync @distributed for iel=1:size(element,1)
        sK[:,iel] = reshape(kron(detjacob[:,iel]',ones(8,3)).*Bu[:,8*(iel-1)+1:8*iel]'*blockdiag(sparse(reshape(DK[:,4*(iel-1)+1],3,3)),sparse(reshape(DK[:,4*(iel-1)+2],3,3)),
        sparse(reshape(DK[:,4*(iel-1)+3],3,3)),sparse(reshape(DK[:,4*(iel-1)+4],3,3)))*Bu[:,8*(iel-1)+1:8*iel],64)
    end
    K=sparse(iK,jK,sdata(sK)[:])
    # GC.gc()
    return K
end
##根据内部应力计算内部等效节点力
function FMat(node::Array{T2},element::Array{T1},Bu::Array{T2},detjacob::Array{T2},σ::Array{T2},edofMat::Array{T1}) where {T1<:Int, T2<:Float64}
    f=SharedArray{Float64,1}(2*size(node,1))
    ff=SharedArray{Float64,2}(8,size(element,1))
    @sync @distributed for iel=1:size(element,1)
        ff[:,iel]=Bu[:,8*(iel-1)+1:8*iel]'.*kron(detjacob[:,iel]',ones(8,3))*reshape(σ[:,4*(iel-1)+1:4*iel],12)
        j=edofMat[iel,:]
        f[j]=f[j]+ff[:,iel]
    end
    # for iel=1:size(element,1)
        # j=edofMat[iel,:]
        # f[j]=f[j]+ff[:,iel]
    # end
    return sdata(f)
end
##根据面上分布荷载计算等效节点力
function FMat_ext(node::Array{T2},conf::T2) where {T1<:Int, T2<:Float64}
    f_ext = zeros(Float64,2*size(node,1))
    Lt_conf = conf
    Rt_conf = -conf
    ##按x坐标值对上部节点排序
    ymax=findall(node[:,2].==maximum(node[:,2])) #1.upper boundary
    xmax=findall(node[:,1].==maximum(node[:,1])) #2.right boundary
    xmin=findall(node[:,1].==minimum(node[:,1])) #3.left boundary
    ##
    elenum_top::Int64 = size(ymax,1).-1
    topNodes = sort(map(tuple, ymax, node[ymax,1]), by = x-> x[2])
    topNodes = [topNodes[i][1] for i=1:size(ymax,1)]
    ##按y坐标值对左部节点排序
    elenum_left::Int64 = size(xmin,1).-1
    leftNodes = sort(map(tuple, xmin, node[xmin,2]), by = x-> x[2])
    leftNodes = [leftNodes[i][1] for i=1:size(xmin,1)]
    ##按y坐标值对右部节点排序
    elenum_right::Int64 = size(xmax,1).-1
    rightNodes = sort(map(tuple, xmax, node[xmax,2]), by = x-> x[2])
    rightNodes = [rightNodes[i][1] for i=1:size(xmax,1)]
    ##
    # for i = 1:elenum_top
    #     tsctr = 2*[topNodes[i],topNodes[i+1]]
    #     f_ext[tsctr] = f_ext[tsctr] + [0.5,0.5]*conf*abs(node[topNodes[i+1],1]-node[topNodes[i],1])
    # end
    for i = 1:elenum_left
        lsctr = 2*[leftNodes[i],leftNodes[i+1]].-1
        f_ext[lsctr] = f_ext[lsctr] + [0.5,0.5]*Lt_conf*abs(node[leftNodes[i+1],2]-node[leftNodes[i],2])
    end
    for i = 1:elenum_right
        rsctr = 2*[rightNodes[i],rightNodes[i+1]].-1
        f_ext[rsctr] = f_ext[rsctr] + [0.5,0.5]*Rt_conf*abs(node[rightNodes[i+1],2]-node[rightNodes[i],2])
    end
    return f_ext
end
function Initial_Hn(d1::Array{T2}, dg1::Array{T2}) where {T1<:Int, T2<:Float64}
    Md = zeros(T2, 4,4,nel)
    H0 = zeros(Float64,4*nel)
    dgd(x) = -2.0 .* ka .* (1.0.-x) ./ (ka .- ka .* (1.0 .- x).^2).^2
    # iel = 1
    # delta=kron(ones(4,4),[1 1 0])
    for iel = 1:nel
        Md[:,:,iel] = reshape(Ndp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*kron(1.0/ls,ones(4)) .+ Bdp[:,4*(iel-1)+1:4*iel].*kron(detjacob[:,iel]',ones(16,1))*kron(1.0*ls,ones(4)),4,4)
        H0[4*(iel-1)+1:4*iel] = Md[:,:,iel]*d1[dedofMat[iel,:]]./(-Nd[:,4*(iel-1)+1:4*iel]'.*kron((dgd.(dg1[4*(iel-1)+1:4*iel]).*detjacob[:,iel])',ones(4,1))*ones(4))
    end
    return H0
end
