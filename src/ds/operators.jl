mutable struct Operators

  V::CuArray{Complex{Float64}}
  K::CuArray{Complex{Float64}}
  wfc::CuArray{Complex{Float64}}

  function Operators(par::Params)
    Vx = @. par.omegaX * par.x
    Vy = @. par.omegaY * par.y
    Vz = @. par.omegaZ * par.z

    V = @. 0.5 * par.mass * (Vx * Vx + Vy * Vy + Vz * Vz)
    V = @. exp(-0.5 * par.dt / ħ * V)

    K = @. exp(-par.dt / ħ * par.k)

    ϕ = @. (par.winding * atan(par.y, par.x)) % (2π)

    wfc = @. exp(-((par.x / par.Rxy / par.a0x) * (par.x / par.Rxy / par.a0x)
                  +(par.y / par.Rxy / par.a0y) * (par.y / par.Rxy / par.a0y)
                  +(par.z / par.Rxy / par.a0z) * (par.z / par.Rxy / par.a0z))
                  + ϕ * im)

    return new(V, K, wfc)
  end

end
