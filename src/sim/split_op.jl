function split_op!(par::Params, opr::Operators, aux::Aux)
  # Half-step in real space
  opr.wfc .*= opr.V

  # FFT to momentum space
  aux.forward_plan * opr.wfc

  # Full step in momentum space
  opr.wfc .*= opr.K

  # iFFT back
  aux.inverse_plan * opr.wfc

  # Final half-step in real space
  opr.wfc .*= opr.V
end
