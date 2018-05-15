const left = 52
const right = 69
data = readtable("sallsac_Hokusai.seq", separator = '\t');

function loadresults()
    #df = DataFrame(img=Any[], n=Any[], tau=Any[], simga=Any[], method=Any[], precl=Any[], metastab=Any[], minstab=Any[])
    df = readtable("results.csv")
    df[:method] = DataArray{Symbol}(df[:method]) # necassary to write symbols into column
    df
end

function filterdata(data, image)
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

function savecl()
    savecl(5,6,50,50,:scaling,200,overwrite=true)
end

function savecl(i, n, tau, sigma; method=:scaling, precl=0, overwrite=false, folder="out", caption=true, symmetrize=false)
    kmeans = method == :kmeans
    path = kmeans ? "$folder/img$i n$n kmeans.png" : "$folder/img$i n$n tau$tau sigma$sigma method$method precl$precl symm$symmetrize.png"

    if isfile(path) && !overwrite
        print(".")
        #PyPlot.imshow(path)
        return
    end

    print("computing $path \n")
    imgdata = filterdata(data,i)

    PyPlot.figure(figsize=(10, 5))

    result = Hokusai.cluster(imgdata[[:fposx, :fposy, :fixdur, :subj]], n, tau=tau, sigma=sigma, precluster=precl, sort=:size, method=method, symmetrize=symmetrize)

    println("saving image")

    imgfile = string(i)[1]

    PyPlot.figure()
    PyPlot.axis(:off)
    Hokusai.plot(result, "$imgfile.jpg", imgdata[1,:width], imgdata[1,:height])

    W = Hokusai.coupling(result)
    metastab = round(trace(W)*100/n, 2)
    minstab  = round(minimum(diag(W))*100, 2)

    if caption
        PyPlot.text(0,0,"img=$i n=$n tau=$tau sigma=$sigma\nmethod=$method precl=$precl metastab=$metastab% minstab=$minstab")
    end

    PyPlot.savefig(path, bbox_inches="tight")

    # 18.4.18: commented this out, dont know what it did
    #push!(df, [i n tau sigma method precl metastab minstab])
    return result
end;

# WIP
function sample(i, n)
    data = df[df[:img] == i && df[:n] == n]
end

