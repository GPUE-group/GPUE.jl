"""
    Operators

GPUE.jl structure for all operator values:
- `v`: Real-valued real-space trapping potential, before exponentiation
- `V`: Real-space trapping potential operator
- `k`: Real-valued momentum-space values, before exponentiation
- `K`: Momentum-space operator
- `wfc`: Wave function
- `ax`: Gauge field values in the X dimension, before exponentiation
- `ay`: Gauge field values in the Y dimension, before exponentiation
- `az`: Gauge field values in the Z dimension, before exponentiation
- `Ax`: Gauge field operator in the X dimension
- `Ay`: Gauge field operator in the Y dimension
- `Az`: Gauge field operator in the Z dimension
"""
mutable struct Operators

  v::Array{Complex{Float64}}
  V::CuArray{Complex{Float64}}

  k::Array{Complex{Float64}}
  K::CuArray{Complex{Float64}}

  wfc::CuArray{Complex{Float64}}

  ax::Array{Complex{Float64}}
  ay::Array{Complex{Float64}}
  az::Array{Complex{Float64}}
  Ax::CuArray{Complex{Float64}}
  Ay::CuArray{Complex{Float64}}
  Az::CuArray{Complex{Float64}}
end

"""
  Operators(par::Params)

Constructs the `Operators` structure.
Initializes the operators given the simulation parameters.
"""
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
  @. V = 0.5 * par.mass * (((ωx2 * par.x * par.x) + (ωy2 * par.y * par.y) + (ωz2 * par.z * par.z)) + (ax * ax + ay * ay + az * az))

  ax = Array(Ax)
  ay = Array(Ay)
  az = Array(Az)

  k = Array(k)
  v = Array(v)

  @. V = CUDAnative.exp(0.5im * par.dt / ħ * V)
  @. K = CUDAnative.exp(1.0im * par.dt * K)
  @. Ax = CUDAnative.exp(1.0im * par.dt * ax * par.px)
  @. Ay = CUDAnative.exp(1.0im * par.dt * ay * par.py)
  @. Az = CUDAnative.exp(1.0im * par.dt * az * par.pz)

  return Operators(v, V, k, K, wfc, ax, ay, az, Ax, Ay, Az)
end

