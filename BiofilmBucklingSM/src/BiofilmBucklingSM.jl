module BiofilmBucklingSM
using LinearAlgebra
using ElasticArrays
using DifferentialEquations
using CairoMakie
using Base

include("Model/DataStructs.jl")
include("Model/GeneralEquations.jl")
include("Model/ODEproblem.jl")
include("Model/ODECallbacks.jl")
include("Model/ProblemSetup.jl")
include("Model/SubstrateSurface.jl")
end # module BiofilmBucklingSM
