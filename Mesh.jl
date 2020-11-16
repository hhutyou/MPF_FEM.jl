
function Mesh()
    eval(:(using DelimitedFiles))
    node::Array{Float64,2}=readdlm("node.txt",',')[:,2:3]
    # fid2=open("data_U1_F1.dat","w")
    # writedlm(fid2,readdlm("node.txt",',')[:,2:3])
    element::Array{Int64,2}=readdlm("element.txt",',')[:,2:5]
    Mat_set1 = try
        union((readdlm("Mat_set1.txt",',',Int))[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    Mat_set2 = try
        union((readdlm("Mat_set2.txt",',',Int))[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    Node_set1 = try
        union((readdlm("loadpoint_d.txt",',',Int))[:])
    catch y
        if isa(y, ArgumentError)
            Int64[]
        end
    end
    # Mat_set1::Array{Int64,1}=union((readdlm("loadele_Hn.txt",','))[:])
    # Mat_set2::Array{Int64,1}=union(readdlm("Mat_set2.txt",',')[:])
    # Node_set1::Array{Int64,1}=union((readdlm("loadpoint_d.txt",','))[:])
    # Mat_set1::Array{Int64,1} = collect(69644:69892)
    # Mat_set2::Array{Int64,1} = collect(69893:70137)
    return node, element, Mat_set1, Mat_set2, Node_set1
end
@info "Formulating the arrays of 'node', 'element', 'Node_set1'  takes"
@time node, element, Mat_set1, Mat_set2, Node_set1=Mesh()
