"""
    terminate!(f::FileData, par::Params, opr::Operators, aux::Aux)

Central function to be called at the end of simulation.
Performs deallocation, and calls finalizers and destructors.
"""
function terminate!(f::FileData, par::Params, opr::Operators, aux::Aux)
  closeFile(f)
end
