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

  x::CuArray{Float64}
  y::CuArray{Float64}
  z::CuArray{Float64}

  nAtoms::Integer
  mass::Float64
  scatterLen::Float64

  a0x::Float64
  a0y::Float64
  a0z::Float64

  Rxy::Float64

  winding::Float64

  k::CuArray{Complex{Float64}}

  dt::Complex{Float64}

  gstate::Bool

  iterations::Integer

end

function Params(xDim=256, yDim=256, zDim=1, boxSize=0.0, omegaX=2*pi, omegaY=2*pi, omegaZ=2*pi, winding=0.0, dt=1e-4, iterations=100)
  if yDim == zDim == 1
    dimnum = 1
  elseif zDim == 1
    dimnum = 2
  else
    dimnum = 3
  end

  nAtoms = 1
  mass = 1.4431607e-25
  scatterLen = 4.76e-9

  a0x = sqrt(ħ / (2 * mass * omegaX))
  a0y = sqrt(ħ / (2 * mass * omegaY))
  a0z = sqrt(ħ / (2 * mass * omegaZ))

  Rxy = (15.0 * nAtoms * scatterLen * sqrt(mass * omegaZ / ħ)) ^ 0.2

  if boxSize > 0
    xMax = yMax = zMax = boxSize
  else
    xMax = 6 * Rxy * omegaX
    yMax = 6 * Rxy * omegaY
    zMax = 6 * Rxy * omegaZ
  end

  dx = 2 * xMax / xDim
  dy = 2 * yMax / yDim
  dz = 2 * zMax / zDim

  x = dimnum < 1 ? [0] : collect(-xMax + dx/2 : dx : xMax)
  y = dimnum < 2 ? [0] : collect(-yMax + dy/2 : dy : yMax)
  z = dimnum < 3 ? [0] : collect(-zMax + dz/2 : dz : zMax)

  dpx = pi / xMax
  dpy = pi / yMax
  dpz = pi / zMax

  pxMax = dpx * (xDim / 2)
  pyMax = dpy * (yDim / 2)
  pzMax = dpz * (zDim / 2)

  px = dimnum < 1 ? [0] : vcat(0 : dpx : pxMax - dpx, -pxMax : dpx : -dpx / 2)
  py = dimnum < 2 ? [0] : vcat(0 : dpy : pyMax - dpy, -pyMax : dpy : -dpy / 2)
  pz = dimnum < 3 ? [0] : vcat(0 : dpz : pzMax - dpz, -pzMax : dpz : -dpz / 2)
  
  k = @. (ħ * ħ / (2.0 * mass)) * (px * px + py * py + pz * pz)

  gstate = real(dt) == dt

  return new(dimnum, xDim, yDim, zDim,
            dx, dy, dz,
            xMax, yMax, zMax,
            omegaX, omegaY, omegaZ,
            CuArray(x), CuArray(y), CuArray(z),
            nAtoms, mass, scatterLen,
            a0x, a0y, a0z, Rxy, winding,
            CuArray(complex(k)),
            dt,
            gstate,
            iterations)
end

