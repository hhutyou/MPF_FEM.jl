##1.计算节点自由度
function xdirect(x::T)::T where T<: Array{Int}
    x=2 .*x .-1
end   ##x方向的自由度
function ydirect(x::T) where T<: Array{Int}
    x=2 .*x
end   ##y方向的自由度
##2.计算主应力
function principle(σ::Array{T},planetype::String) where T<:Float64
    ##
    P, J2, s = invariant(σ,planetype)
    J3 = zeros(Float64,size(J2))
    # J3 = s[1,:].*s[2,:].*stressz-stressz.*stress[3,:].^2
    if planetype=="plane-stress"
        nothing
    elseif planetype=="plane-strain"
        sz =   v.*(σ[1,:].+σ[2,:])-P
        J3 = s[1,:].*s[2,:].*sz.-sz.*σ[3,:].^2
    end
    nu=-3.0*sqrt(3.0)/2.0.*J3./(J2.^1.5)
    function threshold(x)
        if x>=1
            x=1
        elseif x<=-1
            x=-1
        end
        return x
    end
    θ = 1.0/3.0*asin.(threshold.(nu))
    σp = [2.0/sqrt(3.0).*sqrt.(J2).*sin.(θ.+2.0/3.0*π)+P, 2.0/sqrt(3.0).*sqrt.(J2).*sin.(θ.+2.0/3.0*π)+P, 2.0/sqrt(3)*sqrt.(J2).*sin.(θ.+4.0/3.0*π)+P]
    sort!(σp, rev=true)
    # σ₁ = max.(2.0/sqrt(3.0).*sqrt.(J2).*sin.(θ.+2.0/3.0*π).+P, 2.0/sqrt(3.0)*sqrt.(J2).*sin.(θ).+P, 2.0/sqrt(3.0)*sqrt.(J2).*sin.(θ.+4.0/3.0*π).+P)
    # yt = (σ₁ .- P)./(2.0/sqrt(3.0).*sqrt.(J2))
    ##
    return σp[1], σp[2], σp[3]  #sqrt(6.0)/2.0./sin.(θ.+2.0/3.0*π) #max.(sqrt(6.0)/2.0./sin.(θ.+2.0/3.0*π),sqrt(6.0)/2.0./sin.(θ),sqrt(6.0)/2.0./sin.(θ.+4.0/3.0*π))
end
##3.计算应力不变量
function sigma_P(σ::Array{T},v::T,planetype::String) where T<:Float64
    #
    P=Array{T}(undef,size(σ,2))
    if planetype=="plane-stress"
        P = 1.0/3.0.*(σ[1,:]+σ[2,:])
    elseif planetype=="plane-strain"
        stressz =   v.*(σ[1,:]+σ[2,:])
        P = 1.0/3.0.*(σ[1,:]+σ[2,:]+stressz)
    end
    return P
end
function sigma_J2(σ::Array{T},v::T,planetype::String) where T<:Float64
    # stress = σ
    J2=Array{T}(undef,size(σ,1),size(σ,2))
    if planetype=="plane-stress"
        J2 = 1/6*((σ[1,:]-σ[2,:]).^2+σ[1,:].^2+σ[2,:].^2+6.0*σ[3,:].^2)
    elseif planetype=="plane-strain"
        stressz =   v.*(σ[1,:]+σ[2,:])
        J2 = 1/6*((σ[1,:]-σ[2,:]).^2+(σ[1,:]-stressz).^2+
        (σ[2,:]-stressz).^2+6.0*σ[3,:].^2)
    end
    return J2
end
function sigma_s(σ::Array{T},v::T,planetype::String) where T<:Float64
    # stress = σ
    I=[1.0,1.0,0.0]
    # J2=Array{T}(undef,size(σ,2))
    # P=Array{T}(undef,size(σ,2))
    s=Array{T}(undef,4,size(σ,2))
    if planetype=="plane-stress"
        P = 1.0/3.0.*(σ[1,:]+σ[2,:])
        s[1:3,:] = σ[1:3,:].-kron(P', I)
        s[4,:] = -P
        # s[5,:] = s[3,:]
    elseif planetype=="plane-strain"
        # stressz =   v.*(σ[1,:]+σ[2,:])
        P = 1.0/3.0.*(σ[1,:]+σ[2,:]+σ[4,:])
        s[1:3,:] = σ[1:3,:].-kron(P',I)
        s[4,:] = σ[4,:] .- P
        # s[5,:] = s[3,:]
    end
    return s
end
function invariant(σ::Array{T},v::T,planetype::String) where T<:Float64
    # stress = σ
    I=[1.0,1.0,0.0]
    J2=Array{T}(undef,size(σ,2))
    P=Array{T}(undef,size(σ,2))
    s=Array{T}(undef,3,size(σ,2))
    if planetype=="plane-stress"
        P = 1.0/3.0.*(σ[1,:]+σ[2,:])
        s = σ-kron(P', I)
        J2 = 1/6*((σ[1,:]-σ[2,:]).^2+σ[1,:].^2+σ[2,:].^2+6.0*σ[3,:].^2)
        # J3 = zeros(Float64,para[2]*para[3])
    elseif planetype=="plane-strain"
        stressz =   v.*(σ[1,:]+σ[2,:])
        P = 1.0/3.0.*(σ[1,:]+σ[2,:]+stressz)
        # println("typeofP=$(typeof(kron(P',I)))")
        # σ[end,:]=2*σ[end,:]
        s = σ-kron(P',I)
        J2 = 1/6*((σ[1,:]-σ[2,:]).^2+(σ[1,:]-stressz).^2+
        (σ[2,:]-stressz).^2+6.0*σ[3,:].^2)
        # J3 = stress[1,:].*stress[2,:].*stressz-stressz.*stress[3,:].^2
    end
    return P, J2, s
end
#4. 计算应变球量
function operator_tr(epsilon_gauss::Array{T,2},planetype::String) where T<:Float64
    epsilon_tr=Array{T}(undef,size(epsilon_gauss,2))
    if planetype=="plane-stress"
        epsilon_tr= epsilon_gauss[1,:] + epsilon_gauss[2,:] - v/(1-v).*(epsilon_gauss[1,:] + epsilon_gauss[2,:])
    elseif planetype=="plane-strain"
        epsilon_tr = epsilon_gauss[1,:] + epsilon_gauss[2,:]
    end
    return epsilon_tr
end
#5. 计算偏应变(平面应变)
function operator_dev(epsilon_gauss::Array{T,2},planetype::String) where T<:Float64
    # eval(:(using Distributed))
    epsilon_dev = Array{T,2}(undef,size(epsilon_gauss))
    if planetype=="plane-stress"
        epsilon_dev=epsilon_gauss -kron(operator_tr(epsilon_gauss,"plane-stress")',[1/3,1/3,0.0])
    elseif planetype=="plane-strain"
        epsilon_dev=epsilon_gauss -kron(operator_tr(epsilon_gauss,"plane-strain")',[1/3,1/3,0.0])
    end
    epsilon_dev[3,:] .= 0.5*epsilon_dev[3,:]
    return epsilon_dev
end
#5. 计算正应变
function operator_plus(epsilon_gauss::Array{T,2},planetype::String) where T<:Float64
    # eval(:(using Distributed))
    ε = SharedArray{T,2}(size(epsilon_gauss))
    # epsilon = Array{T,2}(undef,2,2)
    if planetype=="plane-stress"
        @sync @distributed for iel = 1:size(epsilon_gauss,2)
            # epsilon = epsilon_gauss[1,:] + epsilon_gauss[2,:]
            V, D = eigen([epsilon_gauss[1,iel] epsilon_gauss[3,iel]/2 0.0; epsilon_gauss[3,iel]/2 epsilon_gauss[2,iel] 0.0; 0.0 0.0 -v/(1.0-v)*(epsilon_gauss[1,iel]+epsilon_gauss[2,iel])])
            # epsilon = (V[1]+abs(V[1]))/2*D[:,1]*D[:,1]' + (V[2]+abs(V[2]))/2*D[:,2]*D[:,2]'
            epsilon = D*diagm([(V[1]+abs(V[1]))/2.0, (V[2]+abs(V[2]))/2.0, (V[3]+abs(V[3]))/2.0])*inv(D)
            ε[:,iel] = [epsilon[1,1], epsilon[2,2], 2*epsilon[1,2]]
        end
    elseif  planetype=="plane-strain"
        @sync @distributed for iel = 1:size(epsilon_gauss,2)
            # epsilon = epsilon_gauss[1,:] + epsilon_gauss[2,:]
            # V, D = eigen([epsilon_gauss[1,iel] epsilon_gauss[3,iel]/2 0.0; epsilon_gauss[3,iel]/2 epsilon_gauss[2,iel] 0.0; 0.0 0.0 0.0])
            V, D = eigen([epsilon_gauss[1,iel] epsilon_gauss[3,iel]/2; epsilon_gauss[3,iel]/2 epsilon_gauss[2,iel]])
            # epsilon = (V[1]+abs(V[1]))/2*D[:,1]*D[:,1]' + (V[2]+abs(V[2]))/2*D[:,2]*D[:,2]'
            epsilon = D*diagm([(V[1]+abs(V[1]))/2.0, (V[2]+abs(V[2]))/2.0])*D'
            ε[:,iel] = [epsilon[1,1], epsilon[2,2], 2*epsilon[1,2]]
        end
    end
    epsilon_gauss = Array(sdata(ε))
    return epsilon_gauss
end
#6. 计算负应变
function operator_minus(epsilon_gauss::Array{T,2},planetype::String) where T<:Float64
    # eval(:(using Distributed))
    ε = SharedArray{T,2}(size(epsilon_gauss))
    # epsilon = Array{T,2}(undef,2,2)
    if planetype=="plane-stress"
        @sync @distributed for iel = 1:size(epsilon_gauss,2)
            # epsilon = epsilon_gauss[1,:] + epsilon_gauss[2,:]
            V, D = eigen([epsilon_gauss[1,iel] epsilon_gauss[3,iel]/2 0.0; epsilon_gauss[3,iel]/2 epsilon_gauss[2,iel] 0.0; 0.0 0.0 -v/(1.0-v)*(epsilon_gauss[1,iel]+epsilon_gauss[2,iel])])
            epsilon = (V[1]-abs(V[1]))/2*D[:,1]*D[:,1]' + (V[2]-abs(V[2]))/2*D[:,2]*D[:,2]'
            # epsilon = D*diagm([(V[1]-abs(V[1]))/2.0, (V[2]-abs(V[2]))/2.0, (V[3]-abs(V[3]))/2.0])*inv(D)
            ε[:,iel] = [epsilon[1,1], epsilon[2,2], 2*epsilon[1,2]]
        end
    elseif  planetype=="plane-strain"
        @sync @distributed for iel = 1:size(epsilon_gauss,2)
            # epsilon = epsilon_gauss[1,:] + epsilon_gauss[2,:]
            V, D = eigen([epsilon_gauss[1,iel] epsilon_gauss[3,iel]/2; epsilon_gauss[3,iel]/2 epsilon_gauss[2,iel]])
            # epsilon = (V[1]-abs(V[1]))/2*D[:,1]*D[:,1]' + (V[2]-abs(V[2]))/2*D[:,2]*D[:,2]'
            epsilon = D*diagm([(V[1]-abs(V[1]))/2.0, (V[2]-abs(V[2]))/2.0])*D'
            ε[:,iel] = [epsilon[1,1], epsilon[2,2], 2*epsilon[1,2]]
        end
    end
    epsilon_gauss = Array(sdata(ε))
    return epsilon_gauss
end
##6. heaviside 函数
function heaviside(x::Float64)
    x = (x + abs(x))/2.0
    return x
end
function mheaviside(x::Float64)
    x = (x - abs(x))/2.0
    return x
end
function hc(x::Float64)
    if x<0.0
        x = 0.0
    else
        x = 1.0
    end
    return x
end
# function xinc!(ret::T) where T<:Int
#    ret = ret+1
# end
# function loopinc_prealloc()
#    ret = Vector{Int}(undef, 3)
#    y = 0
#    for i = 1:10^7
#       xinc!(ret, i)
#       y += ret[2]
#    end
#    return y
# end;
# ret=1
# xinc!(ret)
# a1=rand(10000,10000);b1=rand(10000,10000)
# @time a1.+b1
