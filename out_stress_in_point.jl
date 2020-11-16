using Distributed, Plots,  HDF5, JLD, DelimitedFiles
addprocs(5)
@everywhere using LinearAlgebra, Distributed, SparseArrays, SharedArrays
# using DelimitedFiles
include("Mesh.jl")
include("FemBase.jl")
include("shapeFunc.jl")
include("out_stress.jl")
@load pwd()*"\\numD5.jld" numD5
@load pwd()*"\\numD6.jld" numD6
numd = size(numD2,3)
v = 0.2
nnode=size(node,1)
nel=size(element,1)
sigma = copy(numD2)
## solvle the stress in element node
Tmatrix = inv(Nd[:,1:4])
stress = Array{Float64,2}(undef,3,4*nel)
## matrix used in caculating node stress
ig = Int.(kron(ones(4*nel),[1, 2, 3]))
jg = Int.(kron(element'[:],ones(3)))
ixf = Int.(element'[:])
jxf = Int.(ones(4*nel))
###
freq = SharedArray{Float64}(nnode)
@time @sync @distributed for i = 1:nnode
    freq[i] = size(findall(element[:,:].==i),1)
end
# sum(freq[:].==6)
# @time freq = size(findall(element[:].==p),1)
# @time setdiff(element[:],1)
# @time replace(x->x==1 ? 1 : 0, element[:])
# @time filter(x->x==1, element[:])
# @time findall(element[:,:].==1)
# @time size(findall(element[:].==1))
stress = Array{Float64,2}(undef,3,4*nel)
stress_num = Array{Float64,3}(undef,3,nnode,numd)
xf_num = Array{Float64,2}(undef,nnode,numd)
@time for st = 1:1:1
    for iel = 1:nel
        stress[:,4*(iel-1)+1:4*(iel-1)+4] = numD2[:,4*(iel-1)+1:4*(iel-1)+4,st] * Tmatrix'
    end
    # method 1
    stress_num[:,:,st] = sparse(ig,jg,stress[:]) ./ kron(freq',ones(3))
    # xf_num[:,st] = sparse(ixf,jxf,numD6[:,st]) ./ freq
    # method 2
    # @time stress_num[:,:,st] = sdata(out_stress(sdata(stress)))
end
####
@save "D:\\Columbia_University\\precrack\\alfa=45-conf=10\\stress_num.jld" stress_num
## 单元1 编号 63994
## 节点1 编号 54692 18215
# sigma_p1 = copy(sigma[:,Int(4*63993+1),:])
sigma_p1 = copy(stress_num[:,:,1])
p1_p, p1_J2,  = invariant(sigma_p1,"plane-strain")
p1_s1,   = principle(sigma_p1,"plane-strain")
# output s-p
    # A=[node sqrt.(2 .* p1_J2)]
    A=[node p1_p]
    fid=open("p1p-data.dat","w")
    writedlm(fid,A)
    close(fid)
    # output s1-p
        A=[p1_p p1_s1]
        fid=open("p1s1-data.dat","w")
        writedlm(fid,A)
        close(fid)
        # output chif
            A=[collect(1:numd) xf_num[18215,:]]
            fid=open("p1xf-data.dat","w")
            writedlm(fid,A)
            close(fid)
## 单元2 编号 46409 42764
## 节点2 编号 33013 34558 77674 13459(flaw point)
# sigma_p2 = copy(sigma[:,Int(4*42763+2),:])
sigma_p2 = copy(stress_num[:,13459,:])
p2_p, p2_J2,  = invariant(sigma_p2,"plane-strain")
p2_s1,   = principle(sigma_p2,"plane-strain")
# output s-p
    A=[p2_p sqrt.(2 .* p2_J2)]
    fid=open("p2p-data.dat","w")
    writedlm(fid,A)
    close(fid)
    # output s1-p
        A=[p2_p p2_s1]
        fid=open("p2s1-data.dat","w")
        writedlm(fid,A)
        close(fid)
        # output chif
            A=[collect(1:numd) xf_num[13459,:]]
            fid=open("p2xf-data.dat","w")
            writedlm(fid,A)
            close(fid)
## 单元3 编号 87564 87596
## 节点3 编号 47643
# sigma_p3 = copy(sigma[:,Int(4*87595+2),:])
sigma_p3 = copy(stress_num[:,47643,:])
p3_p, p3_J2,  = invariant(sigma_p3,"plane-strain")
p3_s1,   = principle(sigma_p3,"plane-strain")
# output s-p
    A=[p3_p sqrt.(2 .* p3_J2)]
    fid=open("p3p-data.dat","w")
    writedlm(fid,A)
    close(fid)
    # output s1-p
        A=[p3_p p3_s1]
        fid=open("p3s1-data.dat","w")
        writedlm(fid,A)
        close(fid)
        # output chif
            A=[collect(1:numd) xf_num[47643,:]]
            fid=open("p3xf-data.dat","w")
            writedlm(fid,A)
            close(fid)
# ## 单元4 编号 79331
# sigma_p4 = copy(sigma[:,Int(4*79330+1),:])
# # sigma_p4 = copy(stress_num[element[79331,1]])
# p4_p, p4_J2,  = invariant(sigma_p4,"plane-strain")
# p4_s1,   = principle(sigma_p4,"plane-strain")
# # output s-p
#     A=[p4_p sqrt.(2 .* p4_J2)]
#     fid=open("p4p-data.dat","w")
#     writedlm(fid,A)
#     close(fid)
#     # output s1-p
#         A=[p4_p p4_s1]
#         fid=open("p4s1-data.dat","w")
#         writedlm(fid,A)
#         close(fid)
#         # output chif
#             A=[collect(1:200) numD6[Int(4*79330+1),:]]
#             fid=open("p4xf-data.dat","w")
#             writedlm(fid,A)
#             close(fid)
# ## 单元5 编号 63994
# sigma_p5 = copy(sigma[:,Int(4*63993+1),:])
# # sigma_p5 = copy(stress_num[element[63994,1]])
# p5_p, p5_J2,  = invariant(sigma_p5,"plane-strain")
# p5_s1,   = principle(sigma_p5,"plane-strain")
# # output s-p
#     A=[p5_p sqrt.(2 .* p5_J2)]
#     fid=open("p5p-data.dat","w")
#     writedlm(fid,A)
#     close(fid)
#     # output s1-p
#         A=[p5_p p5_s1]
#         fid=open("p5s1-data.dat","w")
#         writedlm(fid,A)
#         close(fid)
#         # output chif
#             A=[collect(1:200) numD6[Int(4*63993+1),:]]
#             fid=open("p5xf-data.dat","w")
#             writedlm(fid,A)
#             close(fid)
