"""
2.3 -> "2.3"
"""
function _scriptpov(x::Real)
    if isfinite(x)
        return repr(x)
        # return @sprintf "%1.24f" x
    else
        error("numerical value must be finite")
    end
end

"""
[2.3, -5.2] -> "<2.3,-5.2>"
"""
function _scriptpov(x::AbstractVector{<:Real})
    return "<"*join(_scriptpov.(x), ", ")*">"
end

"""
[[2.3, -5.2], [-1.1, 8.2]] -> "<2.3,-5.2>, <-1.1,8.2>"
"""
function _scriptpov(x::AbstractVector{<:AbstractVector{<:Real}})
    return join(_scriptpov.(x), ", ")
end


function _mesh2(M::AbstractBSplineManifold; mesh=(10,10), smooth=true, preindent=0)
    P = bsplinespaces(M)
    p1, p2 = p = degree.(P)
    k1, k2 = k = knots.(P)
    l1, l2 = l = length.(k)
    m1, m2 = m = mesh

    ts1 = unique!(vcat([range(k1[i],k1[i+1],length = m1+1) for i in p1+1:l1-p1-1]...))
    ts2 = unique!(vcat([range(k2[i],k2[i+1],length = m2+1) for i in p2+1:l2-p2-1]...))

    N1 = length(ts1)-1
    N2 = length(ts2)-1

    ts = [[ts1[i], ts2[j]] for i in eachindex(ts1), j in eachindex(ts2)]
    tc = [mean([ts[i1,i2], ts[i1+1,i2], ts[i1,i2+1], ts[i1+1,i2+1]]) for i1 in 1:N1, i2 in 1:N2]

    # M′(t)=ForwardDiff.jacobian(M,t)
    # 𝒆(t) = normalize(cross(M′(t)[1:3,1],M′(t)[1:3,2]))

    𝒑s = M.(ts)
    𝒑c = M.(tc)

    # 𝒆s = 𝒆.(ts)
    # 𝒆c = 𝒆.(tc)

    Ns(i1, i2) = i1 + (N1+1) * (i2-1)
    Nc(i1, i2) = i1 + N1 * (i2-1)

    F1 = [[Ns(i1,i2)-1, Ns(i1+1,i2)-1, (N1+1)*(N2+1)+Nc(i1,i2)-1] for i1 in 1:N1, i2 in 1:N2]
    F2 = [[Ns(i1,i2)-1, Ns(i1,i2+1)-1, (N1+1)*(N2+1)+Nc(i1,i2)-1] for i1 in 1:N1, i2 in 1:N2]
    F3 = [[Ns(i1+1,i2+1)-1, Ns(i1+1,i2)-1, (N1+1)*(N2+1)+Nc(i1,i2)-1] for i1 in 1:N1, i2 in 1:N2]
    F4 = [[Ns(i1+1,i2+1)-1, Ns(i1,i2+1)-1, (N1+1)*(N2+1)+Nc(i1,i2)-1] for i1 in 1:N1, i2 in 1:N2]

    np = (N1+1)*(N2+1) + N1*N2
    nf = 4*N1*N2

    script = "  "^(preindent)
    script *= "mesh2{\n" * "  "^(preindent)
    script *= "  vertex_vectors{\n" * "  "^(preindent)
    script *= "    " * _scriptpov(np) * ", \n" * "  "^(preindent)
    script *= "    " * _scriptpov([𝒑s...]) * ", \n" * "  "^(preindent)
    script *= "    " * _scriptpov([𝒑c...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)

    # if smooth
    #     script *= "  normal_vectors{\n" * "  "^(preindent)
    #     script *= "    " * _scriptpov(np) * ", \n" * "  "^(preindent)
    #     script *= "    " * _scriptpov([𝒆s...]) * ", \n" * "  "^(preindent)
    #     script *= "    " * _scriptpov([𝒆c...]) * "\n" * "  "^(preindent)
    #     script *= "  }\n" * "  "^(preindent)
    # end

    script *= "  face_indices{\n" * "  "^(preindent)
    script *= "    " * _scriptpov(nf) * ", \n" * "  "^(preindent)
    script *= "    " * _scriptpov([F1...]) * "\n" * "  "^(preindent)
    script *= "    " * _scriptpov([F2...]) * "\n" * "  "^(preindent)
    script *= "    " * _scriptpov([F3...]) * "\n" * "  "^(preindent)
    script *= "    " * _scriptpov([F4...]) * "\n" * "  "^(preindent)
    script *= "  }\n" * "  "^(preindent)
    script *= "}\n"

    return script
end

function _spherecylinder(M::AbstractBSplineManifold; preindent=0)
    _spherecylinder(controlpoints(M), preindent=preindent)
end

function _spherecylinder(a::AbstractArray{<:Real}; preindent=0)
    s = size(a)
    d̂ = s[end]
    N = s[1:end-1]
    a_flat = reshape(a,prod(N),d̂)
    a_vec = [a_flat[i,:] for i in 1:prod(N)]
    a′ = reshape(a_vec,N)
    _spherecylinder(a′, preindent=preindent)
end

function _spheres(a::AbstractArray{<:AbstractVector{<:Real}}; preindent=0)
    script = "  "^(preindent)
    script *= "union{\n" * "  "^(preindent)
    for ai in a
        script *= "  sphere{" * _scriptpov(ai) * "," * "0.05" * "}\n" * "  "^(preindent)
    end
    script *= "}\n"
end

function _cylinders(a::AbstractVector{<:AbstractVector{<:Real}}; preindent=0)
    n = length(a)
    script = "  "^(preindent)
    script *= "union{\n" * "  "^(preindent)
    for i in 1:n-1
        script *= "  cylinder{" * _scriptpov(a[i]) * "," * _scriptpov(a[i+1]) * "," * "0.025" * "}\n" * "  "^(preindent)
    end
    script *= "}\n"
end

function _cylinders(a::AbstractArray{<:AbstractVector{<:Real}}; preindent=0)
    n1, n2 = size(a)
    script = "  "^(preindent)
    script *= "union{\n" * "  "^(preindent)
    for i in 1:n1
        script *= _cylinders(a[i,:], preindent=preindent+1)*"\n" * "  "^(preindent)
    end
    for i in 1:n2
        script *= _cylinders(a[:,i], preindent=preindent+1)*"\n" * "  "^(preindent)
    end
    script *= "}\n"
end


function _spherecylinder(a::AbstractArray{<:AbstractVector{<:Real}}; preindent=0)
    script = "  "^(preindent)
    script *= "union{\n"
    script *= "  " * _spheres(a, preindent=preindent+1)
    script *= "  " * _cylinders(a, preindent=preindent+1)
    script *= "}\n"
end

function _save_povray_2d3d(name::String, M::AbstractBSplineManifold; mesh=(10,10), points=true, thickness=1, linecolor=RGB(1,0,0))
    script = "union{\n"
    script *= _mesh2(M,mesh=mesh, preindent=1)
    script *= _spherecylinder(M, preindent=1)
    script *= "}"

    open(name, "w") do f
        write(f,script)
    end

    return nothing
end
