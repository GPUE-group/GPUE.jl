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

    # WFC?

    return new(V, K, WFC)
  end

end
