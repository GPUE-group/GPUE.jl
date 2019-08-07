mutable struct Aux
  
  i::Integer

  forward_plan
  inverse_plan

end

function Aux(opr::Operators)
  forward_plan = plan_fft!(opr.wfc)
  inverse_plan = plan_ifft!(opr.wfc)
  return Aux(0, forward_plan, inverse_plan)
end
