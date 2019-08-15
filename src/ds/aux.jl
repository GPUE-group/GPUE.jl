mutable struct Aux
  
  i::Integer
  density::Float64

  forward_plan_all
  forward_plan_x
  forward_plan_y
  forward_plan_z

  inverse_plan_all
  inverse_plan_x
  inverse_plan_y
  inverse_plan_z
end

function Aux(par::Params, opr::Operators)
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
  
  return Aux(0, -1,
            forward_plan_all,
            forward_plan_x, forward_plan_y, forward_plan_z,
            inverse_plan_all,
            inverse_plan_x, inverse_plan_y, inverse_plan_z)
end
