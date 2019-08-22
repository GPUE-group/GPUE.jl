mutable struct Operators

  V::CuArray{Complex{Float64}}
  K::CuArray{Complex{Float64}}
  wfc::CuArray{Complex{Float64}}

  Ax::CuArray{Complex{Float64}}
  Ay::CuArray{Complex{Float64}}
  Az::CuArray{Complex{Float64}}
end

function Operators(par::Params)
 
  wfc = CuArray{Complex{Float64}}(undef, (par.xDim, par.yDim, par.zDim)[1:par.dimnum])
  Ax = CuArray{Complex{Float64}}(undef, size(wfc))
  Ay = CuArray{Complex{Float64}}(undef, size(wfc))
  Az = CuArray{Complex{Float64}}(undef, size(wfc))
  V = CuArray{Complex{Float64}}(undef, size(wfc))
  K = CuArray{Complex{Float64}}(undef, size(wfc))

  @. Ax = par.y * par.omega * par.omegaX
  @. Ay = -par.x * par.omega * par.omegaY
  @. Az = 0

  @. wfc = CUDAnative.exp(-((par.x / par.Rxy / par.a0x) * (par.x / par.Rxy / par.a0x)
                +(par.y / par.Rxy / par.a0y) * (par.y / par.Rxy / par.a0y)
                +(par.z / par.Rxy / par.a0z) * (par.z / par.Rxy / par.a0z)))

  @. wfc *= CUDAnative.exp(1.0im * ((par.winding * CUDAnative.atan(par.y, par.x)) % (2π)))

  ωx2 = par.omegaX ^ 2
  ωy2 = par.omegaY ^ 2
  ωz2 = par.omegaZ ^ 2

  @. K = 0.5ħ / par.mass * (par.px * par.px + par.py * par.py + par.pz * par.pz)
  @. K = CUDAnative.exp(1.0im * par.dt * K)
  
  @. V = 0.5 * par.mass * (((ωx2 * par.x * par.x) + (ωy2 * par.y * par.y) + (ωz2 * par.z * par.z)) + (Ax * Ax + Ay * Ay + Az * Az))
  @. V = CUDAnative.exp(0.5im * par.dt / ħ * V)

  @. Ax = CUDAnative.exp(1.0im * par.dt * Ax * par.px)
  @. Ay = CUDAnative.exp(1.0im * par.dt * Ay * par.py)
  @. Az = CUDAnative.exp(1.0im * par.dt * Az * par.pz)

  return Operators(V, K, wfc, Ax, Ay, Az)
end

