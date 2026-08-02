import Base: +
import Base: *
import LinearAlgebra: norm
using GLMakie


struct Quarternion{T<:Real}
    w::T
    i::T
    j::T
    k::T
end

+(a::Quarternion, b::Quarternion) = Quarternion(
    a.w + b.w,
    a.i + b.i,
    a.j + b.j,
    a.k + b.k
)

*(a::Quarternion, b::Quarternion) = Quarternion(
    a.w*b.w - a.i*b.i - a.j*b.j - a.k*b.k,
    a.w*b.i + a.i*b.w + a.j*b.k - a.k*b.j,
    a.w*b.j - a.i*b.k + a.j*b.w + a.k*b.i,
    a.w*b.k + a.i*b.j - a.j*b.i + a.k*b.w
)

sqnorm(q::Quarternion) = q.w^2 + q.i^2 + q.j^2 + q.k^2


fman(q::Quarternion,c::Quarternion) = q*q + c

z0 = Quarternion(0.0,0.0,0.0,0.0)

function mandelRecur(c; level=50, start=z0, limit = 100)
    out=start
    limsq=limit^2
    for i in 1:level
        out = fman(out,c)
        if sqnorm(out)>limsq
            return i
        end
    end
    return 0
end

function valtoquartMan(x,i,j; level = 50, start=z0, limit = 100, zCo=0.0)
    fx,fi,fj = convert(Float64,x), convert(Float64,i), convert(Float64,j)
    q = Quarternion(fx,fi,fj,zCo)
    return mandelRecur(q; level, start, limit)
end    

function plotMan(st, en; grain=50, level = 50, start=z0, limit=100, zCo=0.0)
    x = range(st[1], en[1]; length =grain)
    i = range(st[2], en[2]; length =grain)
    j = range(st[3], en[3]; length =grain)
    f(x,i,j) = valtoquartMan(x,i,j; level, start, limit, zCo)
    vol = [f(ix,iy,iz) for ix in x, iy in i, iz in j]
    fig, _ = volume(st[1]..en[1], st[2]..en[2], st[3]..en[3], vol,colorrange = (minimum(vol), maximum(vol)),
        axis=(; type=Axis3, perspectiveness = 0.5,  azimuth = 7.19, elevation = 0.57,  
            aspect = (1,1,1)))

    fig
end

function plotManSlice(st, en; grain=50, level = 50, start=z0, limit=100, Coef=[0.0,0.0,0.0,0.0])
    x = range(st[1], en[1]; length =grain)
    i = range(st[2], en[2]; length =grain)
    j = range(st[3], en[3]; length =grain)
    k(x,i,j) = (Coef[1]*x) + (Coef[2]*i) + (Coef[3]*j) + Coef[4] 
    f(x,i,j) = valtoquartMan(x,i,j; level, start, limit, zCo=k(x,i,j))
    vol = [f(ix,iy,iz) for ix in x, iy in i, iz in j]
    fig, _ = volume(st[1]..en[1], st[2]..en[2], st[3]..en[3], vol,colorrange = (minimum(vol), maximum(vol)),
        axis=(; type=Axis3, perspectiveness = 0.5,  azimuth = 7.19, elevation = 0.57,  
            aspect = (1,1,1)))

    fig
end

plotManSlice([-0.55,-0.4,-0.4], [0.35,0.4,0.4]; grain=600, limit=5, Coef=[0.0,0.0,0.0,1.0])

