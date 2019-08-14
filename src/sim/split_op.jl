function split_op!(par::Params, opr::Operators, aux::Aux)
  # Half-step in real space
  opr.wfc .*= opr.V

  # FFT to momentum space
  aux.forward_plan_all * opr.wfc

  # Full step in momentum space
  opr.wfc .*= opr.K

  # iFFT back
  aux.inverse_plan_all * opr.wfc

  # Final half-step in real space
  opr.wfc .*= opr.V

  # Renormalize if we are in imaginary time
  if (par.gstate)
    opr.wfc ./= sum(abs2.(opr.wfc)) * imag(par.dt)
  end
end
