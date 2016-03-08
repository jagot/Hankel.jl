# Guizar-Sicairos2004JotOSoAA - Computation of quasi-discrete Hankel transforms of integer order for propagating optical wave fields

#=
Calculates an approximation to
\[f_2(\nu) = 2\pi\int\limits_0^\infty dr\; r f_1(r)J_p(2\pi \nu r)\]

as \[\vec{F}_2=\mathsf{T}\vec{F}_1\]

where \[\vec{F}_2(m)=\left\{\frac{f_2\left(\frac{\alpha_{pm}}{2\pi R}\right)V}{|J_{p+1}(\alpha_{pm})|}\right\}\],
\[\vec{F}_1(n)=\left\{\frac{f_1\left(\frac{\alpha_{pn}}{2\pi V}\right)R}{|J_{p+1}(\alpha_{pn})|}\right\}\],
and \[\mathsf{T}_{mn}=\frac{2J_p(\alpha_{pn}\alpha_{pm}/S)}{|J_{p+1}(\alpha_{pn})||J_{p+1}(\alpha_{pm})|S}\].

\(R\) is maximum radius, \(V\) the maximum frequency and \(S=2\pi RV\).

\(f_1\) is assumed to be sampled on these gridpoints.
=#

export Guizar, call, guizar

using GSL

type Guizar
    p::Int
    R
    N::Int
    r
    ν
    V
    Jp1
    T
end

function Guizar(p::Int, R, N::Int)
    α = sf_bessel_zero_Jnu(p, 1:N)
    αNp1 = sf_bessel_zero_Jnu(p, N+1)
    r = α*R/αNp1
    ν = α/(2π*R)
    V = αNp1/(2π*R)
    S = 2π*R*V
    Jp1 = abs(besselj(p+1,α))
    T = 2besselj(p,(α*α')/S)./(S*Jp1*Jp1')
    Guizar(p,R,N,r,ν,V,Jp1,T)
end

import Base.call
function call(gt::Guizar, v::AbstractArray)
    F1 = (v./gt.Jp1)*gt.R
    F2 = gt.T*F1
    gt.ν,F2.*gt.Jp1/gt.V
end

function call(gt::Guizar, f1::Function)
    gt(f1(gt.r))
end

function guizar(f1, p::Int, R, N::Int)
    Guizar(p,R,N)(f1)
end
