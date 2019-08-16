mutable struct Operators

  V::CuArray{Complex{Float64}}
  K::CuArray{Complex{Float64}}
  wfc::CuArray{Complex{Float64}}

  Ax::CuArray{Complex{Float64}}
  Ay::CuArray{Complex{Float64}}
  Az::CuArray{Complex{Float64}}
end

function Operators(par::Params)
  
  ϕ = @. (par.winding * atan(par.y, par.x)) % (2π)

  wfc = @. exp(-((par.x / par.Rxy / par.a0x) * (par.x / par.Rxy / par.a0x)
                +(par.y / par.Rxy / par.a0y) * (par.y / par.Rxy / par.a0y)
                +(par.z / par.Rxy / par.a0z) * (par.z / par.Rxy / par.a0z)))
  wfc = @. wfc * cos(ϕ) + wfc * im * sin(ϕ)

  Ax = CuArray{Complex{Float64}}(undef, size(wfc))
  Ay = CuArray{Complex{Float64}}(undef, size(wfc))
  Az = CuArray{Complex{Float64}}(undef, size(wfc))

  @. Ax = par.y * par.omega * par.omegaX
  @. Ay = -par.x * par.omega * par.omegaY
  @. Az = 0

  Vx = @. par.omegaX * par.x
  Vy = @. par.omegaY * par.y
  Vz = @. par.omegaZ * par.z

  V = @. 0.5 * par.mass * ((Vx * Vx + Vy * Vy + Vz * Vz) + (Ax * Ax + Ay * Ay + Az * Az))
  
  V = @. exp(0.5im * par.dt / ħ * V)
  K = @. exp(1.0im * par.dt / ħ * par.k)
  Ax = @. exp(1.0im * par.dt * Ax * par.px)
  Ay = @. exp(1.0im * par.dt * Ay * par.py)
  Az = @. exp(1.0im * par.dt * Az * par.pz)

  return Operators(V, K, wfc, Ax, Ay, Az)
end

