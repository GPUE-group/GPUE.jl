struct Params
  dimnum::Integer

  xDim::Integer
  yDim::Integer
  zDim::Integer

  dx::Float64
  dy::Float64
  dz::Float64

  xMax::Float64
  yMax::Float64
  zMax::Float64

  omegaX::Float64
  omegaY::Float64
  omegaZ::Float64

  x::CuArray{Complex{Float64}}
  y::CuArray{Complex{Float64}}
  z::CuArray{Complex{Float64}}

  mass::Float64

  k::CuArray{Complex{Float64}}

  dt::Complex{Float64}

  gstate::Bool

  # kwargs?
  function Params(xDim=128, yDim=128, zDim=1, xMax=1.67731e-05, yMax=1.67731e-05, zMax=1.67731e-05, omegaX=1.0, omegaY=1.0, omegaZ=1.0, mass=1.4431607e-25, dt=1e-4)
    if yDim == zDim == 1
      dimnum = 1
    elseif zDim == 1
      dimnum = 2
    else
      dimnum = 3
    end

    dx = 2 * xMax / xDim
    dy = 2 * yMax / yDim
    dz = 2 * zMax / zDim

    x = (dimnum > 0 ? 1 : 0) * collect(-xMax + dx : dx : xMax)
    y = (dimnum > 1 ? 1 : 0) * collect(-yMax + dy : dy : yMax)
    z = (dimnum > 2 ? 1 : 0) * collect(-zMax + dz : dz : zMax)

    dpx = pi / xMax
    dpy = pi / yMax
    dpz = pi / zMax

    pxMax = dpx * (xDim / 2)
    pyMax = dpy * (yDim / 2)
    pzMax = dpz * (zDim / 2)

    px = (dimnum > 0 ? 1 : 0) * vcat(0 : dpx : pxMax - dpx, -pxMax : dpx : -dpx / 2)
    py = (dimnum > 1 ? 1 : 0) * vcat(0 : dpy : pyMax - dpy, -pyMax : dpy : -dpy / 2)
    pz = (dimnum > 2 ? 1 : 0) * vcat(0 : dpz : pzMax - dpz, -pzMax : dpz : -dpz / 2)
    
    println("dpx: $dpx, dpy: $dpy, dpz: $dpz, pxMax: $pxMax, pyMax: $pyMax, pzMax: $pzMax")
 
    k = @. (ħ * ħ / (2.0 * mass)) * (px^2 + py^2 + pz^2)
 
    gstate = real(dt) == dt

    return new(dimnum, xDim, yDim, zDim,
              dx, dy, dz,
              xMax, yMax, zMax,
              omegaX, omegaY, omegaZ,
              CuArray(complex(x)), CuArray(complex(y)), CuArray(complex(z)),
              mass,
              CuArray(complex(k)),
              dt,
              gstate)
  end
end
