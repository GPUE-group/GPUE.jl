"""
    evolve!(f::FileData, par::Params, opr::Operators, aux::Aux)

Start running the main evolution of the system, according to the
settings defined in the Params.
"""
function evolve!(f::FileData, par::Params, opr::Operators, aux::Aux)
  normalize!(par, opr, aux)

  iterations = par.iterations
  while aux.i <= iterations
    if aux.i % par.printSteps == 0 && par.writeOut
      writeWfc(f, par, opr, aux)
      writeAux(f, aux)
      if aux.i == 0
        writeAttributes(f, par)
        writeV(f, par, opr, aux)
        writeK(f, par, opr, aux)
        writeGauge(f, par, opr, aux)
        writeAxes(f, par)
      end
    end

    split_op!(par, opr, aux)

    apply_gauge!(par, opr, aux)

    aux.i += 1
  end
end

