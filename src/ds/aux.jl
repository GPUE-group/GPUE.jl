mutable struct Aux
  
  i::Integer

  forward_plan
  inverse_plan

  function Aux(opr::Operators)
    forward_plan = plan_fft!(opr.wfc)
    inverse_plan = plan_ifft!(opr.wfc)
    return new(1, forward_plan, inverse_plan)
  end

end
