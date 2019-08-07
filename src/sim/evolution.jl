function evolve(f::FileData, par::Params, opr::Operators, aux::Aux)
  writeAttributes(f, par)

  iterations = par.iterations
  while (aux.i < iterations)
    split_op!(par, opr, aux)
    aux.i += 1
  end
end
