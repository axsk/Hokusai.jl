using Colors
using PyPlot

function plot(result::HokusaiResult, imagepath::String, width, height)
    data = result.data
    img = PyPlot.imread(imagepath)

	PyPlot.imshow(img, extent=(0,width,height,0))

	colors = distinguishable_colors(maximum(Array(data[:cluster])), lchoices=linspace(60, 255, 5))
		#lchoices guarantees a miminmal luminance here

	i=0
    by(data, :cluster) do data
		i+=1
		PyPlot.plot(Array(data[:x]), Array(data[:y]),
			"o", markerfacecolor=(colors[i].r, colors[i].g, colors[i].b), alpha=1)
    end
	
	PyPlot.autoscale(tight=true)

	return plt
end
