function applyX!(par::Params, opr::Operators, aux::Aux)
  aux.forward_plan_x * opr.wfc
  opr.wfc .*= opr.Ax
  aux.inverse_plan_x * opr.wfc
end

function applyY!(par::Params, opr::Operators, aux::Aux)
  aux.forward_plan_y * opr.wfc
  opr.wfc .*= opr.Ay
  aux.inverse_plan_y * opr.wfc
end

function applyZ!(par::Params, opr::Operators, aux::Aux)
  aux.forward_plan_z * opr.wfc
  opr.wfc .*= opr.Az
  aux.inverse_plan_z * opr.wfc
end

function apply_gauge!(par::Params, opr::Operators, aux::Aux)
  transformations = [applyX!, applyY!, applyZ!][1 : par.dimnum]

  if isodd(aux.i)
    reverse!(transformations)
  end

  for t in transformations
    t(par, opr, aux)
  end
end

