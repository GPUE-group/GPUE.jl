function evolve(f::FileData, par::Params, opr::Operators, aux::Aux)
  normalize!(par, opr)

  writeAttributes(f, par)

  iterations = par.iterations
  while aux.i < iterations
    if aux.i % par.printSteps == 0
      writeWfc(f, par, opr, aux)
      if aux.i == 0
        writeV(f, par, opr, aux)
        writeK(f, par, opr, aux)
        writeGauge(f, par, opr, aux)
        writeAxes(f, par)
      end
    end

    split_op!(par, opr, aux)

    # apply_gauge!(par, opr, aux)

    aux.i += 1
  end
end
