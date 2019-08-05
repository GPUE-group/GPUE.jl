mutable struct Operators

  V::CuArray{Complex{Float64}}
  K::CuArray{Complex{Float64}}
  WFC::CuArray{Complex{Float64}}

  function Operators(par::Params)
    Vx = @. par.omegaX * par.x
    Vy = @. par.omegaY * par.y
    Vz = @. par.omegaZ * par.z

    V = @. 0.5 * par.mass * (Vx^2 + Vy^2 + Vz^2)
    V = @. exp(-0.5 * par.dt / ħ * V)

    K = @. exp(-par.dt / ħ * par.k)

    winding = 0 # Unsure where to set this
    ϕ = @. (winding * atan(par.y, par.x)) % (2π)

    WFC = @. exp(-((par.x / par.Rxy / par.a0x) ^ 2
                  +(par.y / par.Rxy / par.a0y) ^ 2
                  +(par.z / par.Rxy / par.a0z) ^ 2)) * exp(ϕ * 1im)

    return new(V, K, WFC)
  end

end
