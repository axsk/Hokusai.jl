module Hokusai

using LinearAlgebra
using Distances
using Clustering
using Optim
import Statistics
import PyPlot

include("pccap.jl")
include("cluster.jl")
include("data.jl")
include("automatedMethod.jl")
include("reversibility.jl")

end
