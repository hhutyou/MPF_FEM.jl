#
using Distributed, HDF5, JLD
addprocs(5)
@everywhere using LinearAlgebra, Distributed, SparseArrays, SharedArrays, DelimitedFiles
include("Mesh.jl") #include functions:node, element
@everywhere include("FemBase.jl")
include("boundary.jl")
include("shapeFunc.jl")
include("solvers.jl")
include("K_f_matrix.jl")
include("plasticity.jl")
include("Hn.jl")
include("LSFunc.jl")
include("sigma_plus.jl")
include("sigma_minus.jl")
include("sigma_dev.jl")
include("integration_d.jl")
# Parameters for a case of Indiana Limestone
# elastic parameters
const E0, v = 34000.0, 0.2
const E12 =  E0
const λ0, μ0 = E0*v/((1.0+v)*(1.0-2.0v)), E0/(2.0*(1.0+v))
const G0, Kv0=μ0, λ0+2/3*μ0
# phase field parameters-pure numerical paramers
const ls, k = 2.0, 1e-16
# Initial yield function parameters
const A0, ka, d0 = 0.68, 1.0, 1e-4
const A01 = 0.6
const XX = A0^2*Kv0 + 2.0*μ0
const gc = (75.64*16*(sqrt(6.0)-A0)/9.0/sqrt(3.0))^2 / XX *ka*ls
const gc1 = gc
## Initialization of integrative parameters
  const maxit=30
  const tol=1.0e-3
  nnode=size(node,1)
  nel=size(element,1)
  ##
  planetype = "plane-strain"
  # const λ1, μ1 = E1*v/((1.0+v)*(1.0-2.0v)), E1/(2.0*(1.0+v))
  # const λ2, μ2 = E2*v/((1.0+v)*(1.0-2.0v)), E2/(2.0*(1.0+v))
  # const G1, Kv1=μ1, λ1+2/3*μ1
  # const G2, Kv2=μ2, λ2+2/3*μ2
  # # const D1, D2=E1/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2], E2/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2]##plane-stress
  # const D1, D2=[λ1+2.0μ1 λ1 0.0;λ1 λ1+2.0μ1 0.0;0.0 0.0 μ1], [λ2+2.0μ2 λ2 0.0;λ2 λ2+2.0μ2 0.0;0.0 0.0 μ2]  ##plane-strain
  Mat_1 = map(x->collect(4*(x-1)+1:4*x), Mat_set1)
  Mat_2 = map(x->collect(4*(x-1)+1:4*x), Mat_set2)
  Mat_ind1 = Array{Int64}([])
  Mat_ind2 = Array{Int64}([])
  for i=1:size(Mat_1,1)
      append!(Mat_ind1,getindex(Mat_1,i))
  end
  for i=1:size(Mat_2,1)
      append!(Mat_ind2,getindex(Mat_2,i))
  end
  Mat_ind12 = union(Mat_ind1,Mat_ind2)
  Mat_ind0 = setdiff(1:4*nel,Mat_ind12)
  Mat_set12 = union(Mat_set1,Mat_set2)
  Mat_set0 = setdiff(1:nel,Mat_set12)
  DK = Array{Float64,2}(undef,9,4*size(element,1))
  E = Array{Float64,1}(undef,4*nel)
  AA = Array{Float64,1}(undef,4*nel) ## frictional coefficient
  λ = Array{Float64,1}(undef,4*nel)
  μ = Array{Float64,1}(undef,4*nel)
  Kv = Array{Float64,1}(undef,4*nel)
  for iel in Mat_set0 ## Matrix material
      AA[4*(iel-1)+1:4*iel] .= A0
      E[4*(iel-1)+1:4*iel] .= E0
      λ[4*(iel-1)+1:4*iel] .= E[4*(iel-1)+1:4*iel]*v/((1.0+v)*(1.0-2.0v))
      μ[4*(iel-1)+1:4*iel] .= E[4*(iel-1)+1:4*iel]/(2.0*(1.0+v))
      Kv[4*(iel-1)+1:4*iel] .= λ[4*(iel-1)+1:4*iel] .+ 2.0/3.0 .* μ[4*(iel-1)+1:4*iel]
      DK[:,4*(iel-1)+1:4*iel] .= [λ[4*(iel-1)+1:4*iel]'.+2.0μ[4*(iel-1)+1:4*iel]'; λ[4*(iel-1)+1:4*iel]';
           zeros(1,4); λ[4*(iel-1)+1:4*iel]'; λ[4*(iel-1)+1:4*iel]'.+2.0μ[4*(iel-1)+1:4*iel]';
           zeros(1,4); zeros(1,4); zeros(1,4); μ[4*(iel-1)+1:4*iel]'] ##平面应变
      # DK[:,4*(iel-1)+1:4*iel] = kron(E[4*(iel-1)+1:4*iel]'./(1-v^2), [1; v; 0;v; 1; 0;0; 0; (1-v)/2]) ##平面应力
  end
  for iel in Mat_set12 ## Weak inclusion
      AA[4*(iel-1)+1:4*iel] .= A01
      E[4*(iel-1)+1:4*iel] .= E0 #.*(-log.(1.0.-rand(1))).^(1.0/m)
      λ[4*(iel-1)+1:4*iel] = E[4*(iel-1)+1:4*iel]*v/((1.0+v)*(1.0-2.0v))
      μ[4*(iel-1)+1:4*iel] = E[4*(iel-1)+1:4*iel]/(2.0*(1.0+v))
      Kv[4*(iel-1)+1:4*iel] = λ[4*(iel-1)+1:4*iel] .+ 2.0/3.0 .* μ[4*(iel-1)+1:4*iel]
      DK[:,4*(iel-1)+1:4*iel] = [λ[4*(iel-1)+1:4*iel]'.+2.0μ[4*(iel-1)+1:4*iel]'; λ[4*(iel-1)+1:4*iel]';
           zeros(1,4); λ[4*(iel-1)+1:4*iel]'; λ[4*(iel-1)+1:4*iel]'.+2.0μ[4*(iel-1)+1:4*iel]';
           zeros(1,4); zeros(1,4); zeros(1,4); μ[4*(iel-1)+1:4*iel]']
      # DK[:,4*(iel-1)+1:4*iel] = kron(E[4*(iel-1)+1:4*iel]'./(1-v^2), [1; v; 0;v; 1; 0;0; 0; (1-v)/2]) ##平面应力
  end
  G = deepcopy(μ)
  ##🎺 DK需考虑不均质性💚
## output
u_inc, step_total = -0.005, 120
const numd=step_total ##output number
const aa=Int.(step_total/numd)
begin ##初始化结果储存矩阵
    numD=Array{Float64,2}(undef,size(node,1),numd); numD2=Array{Float64,3}(undef,3,4*size(element,1),numd)# d1,ε
    numD4=Array{Float64,3}(undef,3,4*size(element,1),numd); numD3=Array{Float64,2}(undef,2*size(node,1),numd)#σ, u
    numD5=Array{Float64,3}(undef,3,4*size(element,1),numd);numD6=Array{Float64,3}(undef,3,4*size(element,1),numd)# εᴾ, σc
    numD7=Array{Float64,3}(undef,4,nel,numd); numD8=Array{Float64,2}(undef,maxit,step_total+1)
    Fload1=Array{Float64}(undef,step_total+1); Uload1=Array{Float64}(undef,step_total+1)
    Fload2=Array{Float64}(undef,step_total+1); Uload2=Array{Float64}(undef,step_total+1)
end
##Solution
    ## Confinement
        conf = 0.0
    ## Sovling
    Fload1, Uload1, time_storge, numD, numD2, numD3, numD4, numD5, numD6, numD7, numD8 = solvers(u_inc,conf)
## Output
    for opt=1:aa:step_total
        # A=[node numD3[2:2:end,opt]] ## Displacement
        A=[node numD[:,opt]] ## Phase-field variable
        fid=open("d_data$opt.dat","w")
        StringVariable="TITLE=\"2Dmodel\" VARIABLES=\"X\",\"Y\",\"d1\" ZONE N=$nnode,E=$nel,F=FEPOINT,ET=QUADRILATERAL, "
        write(fid,StringVariable)
        writedlm(fid,A)
        writedlm(fid,element)
        close(fid)
    end
    fid2=open("data_U1_F1.dat","w")
    writedlm(fid2,[Uload1 Fload1])
    close(fid2)
    ##data storage
    @save pwd()*"\\numD.jld" numD
    @save pwd()*"\\numD2.jld" numD2
    @save pwd()*"\\numD3.jld" numD3
    @save pwd()*"\\numD4.jld" numD4
    @save pwd()*"\\numD5.jld" numD5
    @save pwd()*"\\numD6.jld" numD6
    @save pwd()*"\\numD7.jld" numD7
    @save pwd()*"\\numD8.jld" numD8
