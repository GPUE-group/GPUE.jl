"""
    Aux

GPUE.jl structure for miscellaneous values which can be modified dynamically:
- `i`: The iteration counter (starts at 0)
- `forward_plan{all,x,y,z} CuFFT FFT plans along each and every dimension
- `inverse_plan{all,x,y,z} CuFFT iFFT plans along each and every dimension
"""
mutable struct Aux
  
  i::Integer

  forward_plan_all
  forward_plan_x
  forward_plan_y
  forward_plan_z

  inverse_plan_all
  inverse_plan_x
  inverse_plan_y
  inverse_plan_z
end

"""
    Aux(f::FileData, par::Params, opr::Operators)

Constructs the `Aux` structure.
Initializes auxiliary variables, and creates plans for the FFTs.
Will load scalars (`i`) from file if FileData was created from an existing file.
"""
function Aux(f::FileData, par::Params, opr::Operators)
  forward_plan_all = plan_fft!(opr.wfc)
  inverse_plan_all = plan_ifft!(opr.wfc)
  
  if par.dimnum > 0
    forward_plan_x = plan_fft!(opr.wfc, 1)
    inverse_plan_x = plan_ifft!(opr.wfc, 1)
  else
    forward_plan_x = inverse_plan_x = nothing
  end

  if par.dimnum > 1
    forward_plan_y = plan_fft!(opr.wfc, 2)
    inverse_plan_y = plan_ifft!(opr.wfc, 2)
  else
    forward_plan_y = inverse_plan_y = nothing
  end
  
  if par.dimnum > 2
    forward_plan_z = plan_fft!(opr.wfc, 3)
    inverse_plan_z = plan_ifft!(opr.wfc, 3)
  else
    forward_plan_z = inverse_plan_z = nothing
  end

  aux = Aux(0,
            forward_plan_all,
            forward_plan_x, forward_plan_y, forward_plan_z,
            inverse_plan_all,
            inverse_plan_x, inverse_plan_y, inverse_plan_z)

  loadAux!(f, aux)

  return aux
end

