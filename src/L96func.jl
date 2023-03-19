
# using DifferentialEquations
# using Random, Distributions
using StatsBase, Statistics


function lorenz96!(dy,y,p,t) 
    F, K = p
    
    dy[1]=(y[2]-y[K-1])*y[K]-y[1]+F;
    dy[2]=(y[3]-y[K])*y[1]-y[2]+F;
    dy[K]=(y[1]-y[K-2])*y[K-1]-y[K]+F;
    for j=3:K-1
        dy[j]=(y[j+1]-y[j-2])*y[j-1]-y[j]+F;
    end
end

@doc raw"""
circshift_left(a, shifts)


array index shift, reverse diection of in-built circshift function to make sign of index same as in formula.  
"""
function circshift_left(a, shifts)

    return circshift(a, -shifts)
end

@doc raw"""
Lorenz96One_shift!(dy,y,p,t)


One-layer Lorenz 96 equation 

```math

\begin{aligned}
\frac{d x_i}{d t}= & x_{i-1}\left(x_{i+1}-x_{i-2}\right) x_i+F-\frac{h c}{b} \sum_{j=1}^J y_{(j, i)}, \\
\frac{d y_{(j, i)}}{d t}= & c b y_{(j+1, i)}\left(y_{(j-1, i)}-y_{(j+2, i)}\right)  -c y_{(j, i)}+\frac{h c}{b} x_i
\end{aligned}
```


Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
function Lorenz96One_shift!(dy,y,p,t)

    F, K = p;
    dy[:] = circshift(y,-1).*(circshift(y, 1)-circshift(y,-2)) - y .+ F;
end

@doc raw"""
Lorenz96Two_shift_LR!(dy,y,p,t) 


specific formualtion of One-layer Lorenz 96 for lienar regression.
```math 
\begin{aligned}
 \frac{d X_k}{d t}&={-X_{k-1}\left(X_{k-2}-X_{k+1}\right)} {-X_k} {+F} {-h c \bar{Y}_k} \\
 \frac{1}{c} \frac{d Y_{j, k}}{d t}&={-b Y_{j+1, k}\left(Y_{j+2, k}-Y_{j-1, k}\right)}_{-Y_{j, k}}+\frac{h}{J} X_k\\ 
\bar{Y}_k &= \frac{1}{J}\sum^{J}_{j=1}{Y_{k,j}}
\end{aligned}
```
Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
function Lorenz96Two_shift_LR!(dy,y,p,t) 

    # dy_new = zero(y)
    K, F, J, b, p_subgrid = p
    # println(p_subgrid)
    c = p_subgrid.c
    h = p_subgrid.h
    X = y[1:K]
    Y = y[K+1:end]
    Y_mean = dropdims(mean(reshape(Y, K, J), dims=2), dims=2)
    dX = circshift_left(X,-1).*(circshift_left(X, 1)-circshift_left(X,-2)) - X .+ F .- h*c.*Y_mean
    # println(size(Y_mean))
    dY = (circshift_left(Y,1).*(circshift_left(Y, -1)-circshift_left(Y,2)).*b - Y .+ h/J .* repeat(X, J)) .* c
    dy[1:K] = dX
    dy[K+1:end] = dY
    # return dy_new
end

@doc raw"""
Lorenz96Two_shift!(dy,y,p,t) 


Two-layer Lorenz 96

```math

\begin{aligned}
\frac{d x_i}{d t}= & x_{i-1}\left(x_{i+1}-x_{i-2}\right) x_i+F-\frac{h c}{b} \sum_{j=1}^J y_{(j, i)}, \\
\frac{d y_{(j, i)}}{d t}= & c b y_{(j+1, i)}\left(y_{(j-1, i)}-y_{(j+2, i)}\right)  -c y_{(j, i)}+\frac{h c}{b} x_i
\end{aligned}
```
Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
function Lorenz96Two_shift!(dy,y,p,t) 

    F, K, J, h, b, c = p
    # println(p_subgrid)
    dimX = K
    dimY = K*J

    X = y[1:dimX]
    Y = y[dimX+1:end]

    X_subgrid = sum(reshape(Y, :,dimX), dims=1)[:]

    dX = circshift_left(X,-1).*(circshift_left(X, 1)-circshift_left(X,-2)) - X .+ F .- h*c/b * X_subgrid

    X_repeat = repeat(X', J)[:]
    dY = b*c*circshift_left(Y,1).*(circshift_left(Y, -1)-circshift_left(Y,2)) - c*Y .+ h*c/b .* X_repeat
    dy[1:K] = dX
    dy[K+1:end] = dY
    # return dy_new
end

@doc raw"""
Lorenz96Three_shift!(dy,y,p,t) 


Three-layer Lorenz 96
"""
function Lorenz96Three_shift!(dy,y,p,t) 

    F, K, J, I, h, b, c, e, d, g = p
    # println(p_subgrid)
    # c = p_subgrid.c
    # h = p_subgrid.h
    dimX = K
    dimY = K*J
    dimZ = dimY * I
    X = y[1:dimX]
    Y = y[dimX+1:dimY+dimX]
    Z = y[(dimX+dimY+1):end]
    
    X_subgrid = sum(reshape(Y, :, dimX), dims=1)[:]

    Y_subgrid = sum(reshape(Z, :, dimY), dims=1)[:]
    # println(size(Y_subgrid))
    dX = circshift_left(X,-1).*(circshift_left(X, 1)-circshift_left(X,-2)) - X .+ F .- h*c/b.*X_subgrid 
    # println(size(Y_mean))
    dY = (circshift_left(Y,1).*(circshift_left(Y, -1)-circshift_left(Y,2)).*b*c - Y.*c .+ h*c/b .* repeat(X, J)) - h*e/d .* Y_subgrid
    dZ = e*d .*  circshift_left(Z,-1).*(circshift_left(Z, 1)-circshift_left(Z,-2)) - g*e .* Z + h*e/d .* repeat(Y, I)
    dy[1:dimX] = dX
    dy[dimX+1:dimY+dimX] = dY
    dy[dimY+dimX+1: end] = dZ
    # return dy_new
end

@doc raw"""
Lorenz96Two_polyB!(dy,y,p,t)

Subgrid estimation of Two-layer Lorenz 96 by one parameter for all dimension
```math
\begin{aligned}
 \frac{d X_k}{d t}&={-X_{k-1}\left(X_{k-2}-X_{k+1}\right)} {-X_k} {+F} + Poly(X_k) \\
 \frac{1}{c} \frac{d Y_{j, k}}{d t}&={-b Y_{j+1, k}\left(Y_{j+2, k}-Y_{j-1, k}\right)}_{-Y_{j, k}}+\frac{h}{J} X_k
\end{aligned} 
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
function Lorenz96Two_polyB!(dy,y,p,t) 

    K, F, J, b, p_subgrid, B_poly = p

    Bk = B_poly.(y)
    dX = circshift_left(y,-1).*(circshift_left(y, 1)-circshift_left(y,-2)) - y .+ F .+ Bk
 
    dy[:] = dX
end


@doc raw"""
Lorenz96Two_polyBk!(dy,y,p,t)


Subgrid estimation of Two-layer Lorenz 96 by k parameters for k dimension
```math
\begin{aligned}
 \frac{d X_k}{d t}&={-X_{k-1}\left(X_{k-2}-X_{k+1}\right)} {-X_k} {+F} + Poly(X) \\
 \frac{1}{c} \frac{d Y_{j, k}}{d t}&={-b Y_{j+1, k}\left(Y_{j+2, k}-Y_{j-1, k}\right)}_{-Y_{j, k}}+\frac{h}{J} X_k
\end{aligned} 
```

Objects of this type are callable with the signature `(du, u, p, t)` which
performs an in-place update of the derivatives `du` of the system.
"""
function Lorenz96Two_polyBk!(dy,y,p,t) 
    
    K, F, J, b, p_subgrid, B_poly = p

    Bk = zeros(K)
    for i = 1:K
        Bk[i] = B_poly[i](y[i])
    end

    dX = circshift_left(y,-1).*(circshift_left(y, 1)-circshift_left(y,-2)) - y .+ F .+ Bk
    dy[:] = dX

end