function boundary(node::Array{T2},element::Array{T1}) where {T1<:Int, T2<:Float64}
    # ymax=findall(node[:,2].==maximum(node[:,2])) #1.upper boundary
    # ymin=findall(node[:,2].==minimum(node[:,2])) #2.bottom boundary
    # xmax=findall(node[:,1].==maximum(node[:,1])) #3.right boundary
    # xmin=findall(node[:,1].==minimum(node[:,1])) #4.left boundary
    # xmintop::Array{T1}=findall((node[:,2].>=10.1) .& (node[:,1].==minimum(node[:,1])))## 5.left-top
    # xminlow::Array{T1}=findall((node[:,2].<=9.9) .& (node[:,1].==minimum(node[:,1])))## 6.left-lower
    # xmaxtop::Array{T1}=findall((node[:,2].>=10.1) .& (node[:,1].==maximum(node[:,1]))) ## 7.up-right
    # xmaxlow::Array{T1}=findall((node[:,2].<=9.9) .& (node[:,1].==maximum(node[:,1]))) ## 7.low-right
    # eval(:(using DelimitedFiles))
    fix_x = try
        union(readdlm("fix_x.txt",',',Int)[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    fix_y = try
        union(readdlm("fix_y.txt",',',Int)[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    load_x = try
        union(readdlm("load_x.txt",',',Int)[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    load_y = try
        union(readdlm("load_y.txt",',',Int)[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    ##
    fixeddofs=union(ydirect(fix_y),xdirect(fix_x))
    # fixeddofs=union(ydirect(ymin),ydirect(ymax),xdirect(ymin))
    ##
    # function xdirect(x::T) where T<: Nothing
    #     x=Int64[]
    # end
    loaddofsx::Array{T1}=xdirect(load_x)
    loaddofsy::Array{T1}=ydirect(load_y)
    loaddofs::Array{T1}=union(loaddofsx,loaddofsy)
    # loaddofs::Array{T1}=union(xdirect(ymax))
    freedofs::Array{T1}=setdiff(1:2size(node,1),fixeddofs,loaddofs)
    loaddofs_d::Array{T1} = Node_set1
    freedofs_d::Array{T1} = setdiff(1:size(node,1), Node_set1)
#     # xmaxtop=intersect(find(node[:,2].>=102.5),find((node[:,1].==maximum(node[:,1])))) ## 8.top-right
    return fixeddofs, loaddofs, freedofs, loaddofs_d, freedofs_d
end
@info "Formulating the arrays of freedom degree takes"
@time fixeddofs, loaddofs, freedofs, loaddofs_d, freedofs_d=boundary(node,element)
