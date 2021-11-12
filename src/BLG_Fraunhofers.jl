module BLG_Fraunhofers

module EncapsulatedBLG

# Use README as the docstring of the module:
using ProgressMeter: showprogress
using Quantica
using Dates: CONVERSION_SPECIFIERS
using DataStructures: left_rotate
using DataFrames: _broadcast_unalias_helper
@doc read(joinpath(dirname(@__DIR__), "README.md"), String) EncapsulatedBLG

using Requires

function __init__()
    #   @require VegaLite = "112f6efa-9a02-5b7d-90c0-432ed331239a" include("vegalite_plots.jl")
    #   @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0" include("cairomakie_plots.jl")
end

using Interpolations, Parameters, DataFrames, CSV, 
    SharedArrays,  Random, StaticArrays, LinearAlgebra
using Arpack, ArnoldiMethod, Distributed
using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B, c_0

# using KrylovKit, DataFrames, Dates, CSV, SharedArrays, Statistics
# using Arpack, ArnoldiMethod
# using Parameters, LinearAlgebra, Random, ProgressMeter, Random, StaticArrays
# using FFTW, JLD, DataStructures, SymDict, ParallelDataTransfer, DelimitedFiles, ArnoldiMethod
# using Clustering
# using Dierckx
# using LsqFit, Statistics, Random, Optim, Waveforms

export Params
export paramhams


include("regions.jl")
include("model.jl") 
include("hamiltonians.jl")
include("transport.jl")
include("adaptive_maxfinder.jl")

end

end

# include("densityplot.jl")
# include("spin.jl")
# include("KPMdensity.jl")
# include("semiinfinite.jl")
# include("sevmethodfunctions.jl")
# include("KPMtransport.jl")