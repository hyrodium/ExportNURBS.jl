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
function save_svg(name::String, M::BSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true)
    if split(name,'.')[end] ≠ "svg"
        name = name * ".svg"
    end
    if dim(M) == 1
        _save_1d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points)
    elseif dim(M) == 2
        _save_2d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points)
    else
        error("the dimension of B-spline manifold must be 2 or less")
    end
end

"""
export png file
"""
function save_png(name::String, M::BSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true)
    if split(name,'.')[end] ≠ "png"
        name = name * ".png"
    end
    if dim(M) == 1
        _save_1d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points)
    elseif dim(M) == 2
        _save_2d2d(name, M, up=up, down=down, right=right, left=left, zoom=zoom, mesh=mesh, unitlength=unitlength, points=points)
    else
        error("the dimension of B-spline manifold must be 2 or less")
    end
end


function _save_2d2d(name::String, M::BSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=(10,10), unitlength=100, points=true)
    step = unitlength
    p¹,p² = p = degree.(M.bsplinespaces)
    k¹,k² = k = knots.(M.bsplinespaces)
    𝒂 = M.controlpoints
    n¹,n² = n = length.(k)-p.-1
    𝒑(u) = mapping(M,u)

    K¹,K² = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:2]
    N¹,N² = length.(K).-1
    m¹,m² = mesh

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    setline(zoom)
    background("white")

    sethue(1,.5,.5) # Pale Red
    drawbezierpath(BezierPath(vcat(
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[1]]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[end],u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,K²[end]]),K¹[end-i+1],K¹[end-i]))...) for i ∈ 1:N¹],
        [BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([K¹[1],u²]),K²[end-i+1],K²[end-i]))...) for i ∈ 1:N²]
    )),:fill,close=true)

    sethue("red") # Red
    for u¹ ∈ range(K¹[1],stop=K¹[end],length=m¹+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u²->𝒑([u¹,u²]),K²[i],K²[i+1]))...) for i ∈ 1:N²]),:stroke)
    end
    for u² ∈ range(K²[1],stop=K²[end],length=m²+1)
        drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹,u²]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)
    end

    if points
        sethue(.1,.1,.1) # Dark Gray
        setline(zoom)
        CtrlPts = [LxrPt(𝒂[i,j,:],step) for i ∈ 1:size(𝒂)[1], j ∈ 1:size(𝒂)[2]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)

        sethue(.3,.3,.3) # Light Gray
        for i ∈ 1:n¹
            poly(CtrlPts[i,:], :stroke)
        end
        for j ∈ 1:n²
            poly(CtrlPts[:,j], :stroke)
        end
    end
    finish()
end


function _save_1d2d(name::String, M::BSplineManifold; up=5, down=-5, right=5, left=-5, zoom=1, mesh=10, unitlength=100, points=true)
    step = unitlength
    p¹, = p = degree.(M.bsplinespaces)
    k¹, = k = knots.(M.bsplinespaces)
    𝒂 = M.controlpoints
    n¹, = n = length.(k)-p.-1
    𝒑(u) = mapping(M,u)

    K¹, = K = [unique(k[i][1+p[i]:end-p[i]]) for i ∈ 1:1]
    N¹, = length.(K).-1
    m¹, = mesh

    Drawing(step*(right-left),step*(up-down),name)
    Luxor.origin(-step*left,step*up)
    setline(2*zoom)
    background("white")

    sethue("red") # Red
    drawbezierpath(BezierPath([BezierPathSegment(map(p->LxrPt(p,step),BézPts(u¹->𝒑([u¹]),K¹[i],K¹[i+1]))...) for i ∈ 1:N¹]),:stroke)

    if points
        sethue(.1,.1,.1) # Dark Gray
        setline(zoom)
        CtrlPts = [LxrPt(𝒂[i,:],step) for i ∈ 1:size(𝒂)[1]]
        map(p->circle(p,3*zoom,:fill), CtrlPts)

        sethue(.3,.3,.3) # Light Gray
        poly(CtrlPts[:], :stroke)
    end
    finish()
end
