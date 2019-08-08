function evolve(f::FileData, par::Params, opr::Operators, aux::Aux)
  writeAttributes(f, par)

  iterations = par.iterations
  while aux.i < iterations
    split_op!(par, opr, aux)

    apply_gauge!(par, opr, aux)

    if aux.i % par.printSteps == 0
      writeWfc(f, par, opr, aux)
      if aux.i == 0
        writeV(f, par, opr, aux)
        writeK(f, par, opr, aux)
        writeAxes(f, par)
      end
    end

    aux.i += 1
  end
end