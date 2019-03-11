# initialize data
using DataFrames, CSV, PyPlot

datapath = joinpath(@__DIR__, "..", "data")
#using RData

global DATA

function __init__()
    try
        global DATA = CSV.read(joinpath(datapath, "sallsac_Hokusai.seq"), delim = '\t')
    catch
        @warn("could not load data")
    end
end


# TODO: think about a clean filtering interface
# mirrored: 0- not mirrored, 1 - mirrored
# startPos: 0-all, 52-left, 69-right
function filterdata(data::DataFrame, image, mirrored::Int=0, startPos::Int=0)
    # select the corresponding image
    if mirrored == 0
        data = data[data[:image] .== "$image.jpg", :]
    else
        data = data[data[:image] .== string(image,"_mirrored.jpg"), :]
    end

    if startPos != 0
        data = data[data[:fixcrosspos] .== startPos,:]
    end

    width, height = data[1,:width], data[1,:height]
    dx, dy = (1280-width)/2, (1024-height)/2

    data[:fposx] = data[:fposx] .- dx
    data[:fposy] = data[:fposy] .- dy
    data = data[((data[:fposx] .> 0) .& (data[:fposy] .> 0) .& (data[:fposx] .< width) .& (data[:fposy] .< height)) , :]


    # experimental: "disable" grouping by subjects (MUTATING!)
    # data[:subj] = 0

    # TODO: this is not the best way to handle this
    # - it would be nice to be able to just overlay un- and mirrored data in a comparable way, which this does not do

    # mirror x coordinates of mirrored version to make it comparable
    if occursin("mirrored", string(image)) || (mirrored == 1)
        data[:fposx] = width .- data[:fposx]
    end

    return data
end

function TimeSeries(data::DataFrame)
    ts = TimeSeries[]
    by(data, :subj) do d
        times = Array(cumsum(d[:fixdur])) ::Vector
        points = Array(hcat(d[:fposx], d[:fposy])) ::Matrix{Float64}
        push!(ts, TimeSeries(times, points))
        nothing # return value for by construct / hotfix
    end
    ts
end

## convenience wrapper for plotting
function run(ts::Union{TimeSeries, Vector{TimeSeries}}, n, sigma, tau, person; kwargs...)
    timeseries = person == 0 ? ts : ts[person]
    ass, chi = cluster(timeseries, n, sigma, tau; kwargs...)
    figure()
    plot(timeseries, ass)
    ass, chi
end

run(d::DataFrame, n, sigma, tau, person; kwargs...) = run(TimeSeries(d), n, sigma, tau, person; kwargs...)

# person: 0- all corresponding subjects, else choose a subj from 1:length(corresponding subjects)
# mirrored: 0- not mirrored, 1 - mirrored
# startPos: 0-all, 52-left, 69-right
function run(img::Integer, n, sigma, tau, person::Int = 0, mirrored::Int=0, startPos::Int=0; kwargs...)
    run(filterdata(DATA, img, mirrored, startPos), n, sigma, tau, person; kwargs...)
    plotimg(img)
end

## plot functions
function plot(ts::Union{TimeSeries, Vector{TimeSeries}}, ass)
    ps = points(ts)
    for i = 1:maximum(ass)
        j = find(ass .== i)
        PyPlot.scatter(ps[j,1], ps[j,2], label = i)
    end
    PyPlot.legend()
    gcf()
end

plotimg(img::Integer) = plotimg(joinpath(datapath,"$img.jpg"))

function plotimg(path::String)
    PyPlot.imread(path) |> PyPlot.imshow
    PyPlot.gcf()
end
