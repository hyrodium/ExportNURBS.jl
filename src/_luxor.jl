"""
control points of Bézier curve from given function and range
"""
function BézPts(𝒑,a,b)
    𝒑(a),
    3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)-3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27)/2,
    -3*(𝒑(2*a/3+b/3)-(8*𝒑(a)+𝒑(b))/27)/2+3*(𝒑(a/3+2*b/3)-(𝒑(a)+8*𝒑(b))/27),
    𝒑(b)
end

function LxrPt(p::Array{T,1},step) where T<:Real
    Point(step*[1,-1].*p...)
end

"""
export svg file
"""
function save_svg(name::String, M::AbstractBSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true, thickness=1, backgroundcolor=RGB(1,1,1), linecolor=RGB(1,0,0))
    if split(name,'.')[end] ≠ "svg"
        name = name * ".svg"
    end
    if dim(M) == 1
        _save_luxor_1d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points, thickness=thickness, backgroundcolor=backgroundcolor, linecolor=linecolor)
    elseif dim(M) == 2
        _save_luxor_2d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points, thickness=thickness, backgroundcolor=backgroundcolor, linecolor=linecolor)
    else
        error("the dimension of B-spline manifold must be 2 or less")
    end
end

"""
export png file
"""
function save_png(name::String, M::AbstractBSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true, thickness=1, backgroundcolor=RGB(1,1,1), linecolor=RGB(1,0,0))
    if split(name,'.')[end] ≠ "png"
        name = name * ".png"
    end
    if dim(M) == 1
        _save_luxor_1d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points, thickness=thickness, backgroundcolor=backgroundcolor, linecolor=linecolor)
    elseif dim(M) == 2
        _save_luxor_2d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points, thickness=thickness, backgroundcolor=backgroundcolor, linecolor=linecolor)
    else
        error("the dimension of B-spline manifold must be 2 or less")
    end
end

"""
export png file
"""
function save_png(name::String, M::AbstractBSplineManifold, colors::Array{<:Colorant,2}; up=5, down=-5, right=5, left=-5, zoom=1, unitlength=100)
    if split(name,'.')[end] ≠ "png"
        name = name * ".png"
    end

    if dim(M) == 2
        _save_luxor_2d2d_color(name, M, colors, up=up, down=down, right=right, left=left, zoom=zoom, unitlength=unitlength)
    else
        error("the dimension of B-spline manifold must be 2")
    end
end

"""
export png file
"""
function save_png(name::String, M::AbstractBSplineManifold, colorfunc::Function; up=5, down=-5, right=5, left=-5, zoom=1, unitlength=100)
    if split(name,'.')[end] ≠ "png"
        name = name * ".png"
    end

    if dim(M) == 2
        _save_luxor_2d2d_color(name, M, colorfunc, up=up, down=down, right=right, left=left, zoom=zoom, unitlength=unitlength)
    else
        error("the dimension of B-spline manifold must be 2")
    end
end


function _save_luxor_2d2d(name::String, M::AbstractBSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true, thickness=1, backgroundcolor=RGB(1,1,1), linecolor=RGB(1,0,0))
    step = unitlength
    P1, P2 = P = collect(bsplinespaces(M))
    p¹, p² = p = degree.(P)
    k¹, k² = k = knots.(P)
    𝒂 = controlpoints(M)
    n¹, n² = n = length.(k)-p.-1
    𝒑(u) = M(u)

    K¹,K² = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N¹,N² = length.(K).-1
    m¹,m² = mesh

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    setline(thickness)
    background(backgroundcolor)

    setcolor(1,.5,.5) # Pale Red
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[1]]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[end],u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[end]]),K¹[end-i+1],K¹[end-i]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[1],u²]),K²[end-i+1],K²[end-i]))...) for i ∈ 1:N²]
    )),:fill,close=true)

    setcolor(linecolor) # Red
    for u¹ ∈ range(K¹[1],stop=K¹[end],length=m¹+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([u¹,u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²]),:stroke)
    end
    for u² ∈ range(K²[1],stop=K²[end],length=m²+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,u²]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)
    end

    if points
        CtrlPts = [LxrPt(𝒂[i,j,:],step) for i ∈ 1:size(𝒂)[1], j ∈ 1:size(𝒂)[2]]

        setcolor(RGB(.3,.3,.3)) # Light Gray
        setline(thickness)
        for i ∈ 1:n¹
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n²
            poly(CtrlPts[:,j], :stroke)
        end

        setcolor(RGB(.1,.1,.1)) # Dark Gray
        map(p->circle(p,3*thickness,:fill), CtrlPts)
    end
    finish()
    return nothing
end

function _save_luxor_1d2d(name::String, M::AbstractBSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=10, unitlength=100, points=true, thickness=1, backgroundcolor=RGB(1,1,1), linecolor=RGB(1,0,0))
    step = unitlength
    P1, = P = collect(bsplinespaces(M))
    p¹, = p = degree.(P)
    k¹, = k = knots.(P)
    𝒂 = controlpoints(M)
    n¹, = n = length.(k)-p.-1
    𝒑(u) = M(u)

    K¹, = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:1]
    N¹, = length.(K).-1
    m¹, = mesh

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    setline(2*thickness)
    background(backgroundcolor)

    setcolor(linecolor)
    drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)

    if points
        CtrlPts = [LxrPt(𝒂[i,:],step) for i ∈ 1:size(𝒂)[1]]

        setcolor(RGB(.3,.3,.3)) # Light Gray
        setline(thickness)
        poly(CtrlPts[:], :stroke)

        setcolor(RGB(.1,.1,.1)) # Dark Gray
        map(p->circle(p,3*thickness,:fill), CtrlPts)
    end
    finish()
    return nothing
end

function _save_luxor_2d2d_color(name::String, M::AbstractBSplineManifold, colors::Array{T,2} where T <: Colorant; up=5, down=-5, right=5, left=-5, zoom=1, unitlength=100)
    P = collect(bsplinespaces(M))
    colorfunc(u) = sum(bsplinebasis(P,u).*colors)
    _save_luxor_2d2d_color(name, M, colorfunc; up=up, down=down, right=right, left=left, zoom=zoom, unitlength=unitlength)
end

function _save_luxor_2d2d_color(name::String, M::AbstractBSplineManifold, colorfunc::Function; up=5, down=-5, right=5, left=-5, zoom=1, unitlength=100)
    mesh = 10

    step = unitlength
    P = collect(bsplinespaces(M))
    p¹, p² = p = degree.(P)
    k¹, k² = k = knots.(P)
    𝒂 = controlpoints(M)
    n¹, n² = n = length.(k)-p.-1
    𝒑(u) = M(u)

    D = [k[i][1+p[i]]..k[i][end-p[i]] for i in 1:2]

    K¹ = unique(vcat([collect(range(k¹[i], k¹[i+1], length=mesh+1)) for i in 1+p¹:length(k¹)-p¹-1]...))
    K² = unique(vcat([collect(range(k²[i], k²[i+1], length=mesh+1)) for i in 1+p²:length(k²)-p²-1]...))

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    setline(2*zoom)
    background(RGBA(0.0,0.0,0.0,0.0))

    for I₁ ∈ 1:length(K¹)-1, I₂ ∈ 1:length(K²)-1
        BézPth=BezierPath([
                BezierPathSegment(map(p->LxrPt(p,step),BézPts(t->𝒑([t,K²[I₂]]),K¹[I₁],K¹[I₁+1]))...),
                BezierPathSegment(map(p->LxrPt(p,step),BézPts(t->𝒑([K¹[I₁+1],t]),K²[I₂],K²[I₂+1]))...),
                BezierPathSegment(map(p->LxrPt(p,step),BézPts(t->𝒑([t,K²[I₂+1]]),K¹[I₁+1],K¹[I₁]))...),
                BezierPathSegment(map(p->LxrPt(p,step),BézPts(t->𝒑([K¹[I₁],t]),K²[I₂+1],K²[I₂]))...)])
        mesh1 = Luxor.mesh(BézPth, [
            colorfunc([K¹[I₁], K²[I₂]]),
            colorfunc([K¹[I₁+1], K²[I₂]]),
            colorfunc([K¹[I₁+1], K²[I₂+1]]),
            colorfunc([K¹[I₁], K²[I₂+1]])
            ])
        setmesh(mesh1)
        box(LxrPt([right+left,up+down]/2,step), (right-left)*step,(up-down)*step,:fill)
    end
    finish()
    return nothing
end
