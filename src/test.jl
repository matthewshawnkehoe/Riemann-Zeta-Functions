using Plots
M = [ 1 2 3 ; 3 2 1; 4 8 2];
# Create a heatmap with a grayscale colormap
heatmap(M, c=:grays)

using ColorSchemes, Images

function imagesc(data::AbstractMatrix{<:Real};
                 colorscheme::ColorScheme=ColorSchemes.viridis,
                 maxsize::Integer=500, rangescale=:extrema)
    s = maximum(size(data))
    if s > maxsize
        return imagesc(imresize(data, ratio=maxsize/s);   # imresize from Images.jl
                       colorscheme, maxsize, rangescale)
    end
    return get(colorscheme, data, rangescale) # get(...) from ColorSchemes.jl
end

x = range(0,1,length=2*10^4)
y = x'
img = atan.(x, y); # a 20000x20000 Float64 array
@time imagesc(img)

# using PyPlot
# Set the backend to Qt5Agg
# ENV["MPLBACKEND"] = "Qt5Agg"
# using PyCall

M = [1 2 3; 3 2 1; 4 8 2]

# Plot the matrix with a grayscale colormap
imshow(M, cmap="gray")
colorbar()