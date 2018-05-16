type HokusaiResult
    data
    assignments
    n
    tau
    sigma
    P
end

function readdata!(datafile = joinpath(datapath, "sallsac_Hokusai.seq"))
    global data = readtable(datafile, separator = '\t')
end

# initialize data
datapath = joinpath(@__DIR__, "..", "data")
global data
try
    readdata!()
catch
    warn("could not load data")
end

# TODO: add filters for mirrored, starting location, ..
function filterdata(data::DataFrame, image)
    # select the corresponding image
    data = data[data[:image] .== "$image.jpg", :]

    # consider only one starting position
    const left = 52
    const right = 69
    #data = data[data[:fixcrosspos] .== left]

    # "disable" grouping by subjects
    #data[:subj] = 0

    # scale/filter coordinates to 0->width, 0->height
    width, height = data[1,:width], data[1,:height]
    dx, dy = (1280-width)/2, (1024-height)/2
    data[:fposx] = data[:fposx] - dx
    data[:fposy] = data[:fposy] - dy
    data = data[((data[:fposx] .> 0) .& (data[:fposy] .> 0) .& (data[:fposx] .< width) .& (data[:fposy] .< height)) , :]

    # mirror x coordinates of mirrored version to make it comparable
    if ismatch(r"mirrored", string(image))
        data[:fposx] = width - data[:fposx]
    end

    return data
end;

# convenience function for clustering, plotting and saving a run (written for batch use)
function savecl(i, n, tau, sigma; method=:scaling, precl=0, overwrite=false, folder="out", caption=true, symmetrize=false)
    kmeans = method == :kmeans
    path = kmeans ? "$folder/img$i n$n kmeans.png" : "$folder/img$i n$n tau$tau sigma$sigma method$method precl$precl symm$symmetrize.png"

    if isfile(path) && !overwrite
        print("Clustering already existing. Skipping...")
        return
    end

    print("computing $path \n")
    imgdata = filterdata(data,i)

    PyPlot.figure(figsize=(10, 5))

    result = Hokusai.cluster(imgdata[[:fposx, :fposy, :fixdur, :subj]], n, tau=tau, sigma=sigma, precluster=precl, sort=:size, method=method, symmetrize=symmetrize)

    imgfile = string(i)[1]

    PyPlot.figure()
    PyPlot.axis(:off)
    Hokusai.plot(result, joinpath(datapath, "$imgfile.jpg"), imgdata[1,:width], imgdata[1,:height])

    if caption
        PyPlot.text(0,0,"img=$i n=$n tau=$tau sigma=$sigma\nmethod=$method precl=$precl")
    end

    println("saving image")
    mkpath(folder)
    PyPlot.savefig(path, bbox_inches="tight")

    # was used to store metadata
    #push!(df, [i n tau sigma method precl metastab minstab])
    return result
end

# legacy cluster wrapper
function cluster(data::DataFrame, n; tau=50, sigma=100, precluster=0, sort=:size, method=:scaling, symmetrize=false)
    ts = TimeSeries[]
    by(data, :subj) do d
        times = Array(cumsum(d[:fixdur])) ::Vector
        points = Array(hcat(d[:fposx], d[:fposy])) ::Matrix{Float64}
        push!(ts, TimeSeries(times, points))
        0 # since the return value from push crashes the "by" construction
    end
    ass = cluster(ts, n, tau=tau, sigma=sigma, precluster=precluster, sort=sort, method=method, symmetrize=symmetrize)

    data[:cluster] = ass
    HokusaiResult(data, ass, n, tau, sigma, nothing)
end

function plot(result::HokusaiResult, imagepath::String, width, height)
    data = result.data
    img = PyPlot.imread(imagepath)

    PyPlot.imshow(img, extent=(0,width,height,0))

    colors = distinguishable_colors(maximum(Array(data[:cluster])), lchoices=linspace(60, 255, 5))
    #lchoices guarantees a miminmal luminance here

    i=0
    by(data, :cluster) do data
	i+=1
	PyPlot.plot(Array(data[:fposx]), Array(data[:fposy]),
		    "o", markerfacecolor=(colors[i].r, colors[i].g, colors[i].b), alpha=1)
    end

    PyPlot.autoscale(tight=true)

    return plt
end