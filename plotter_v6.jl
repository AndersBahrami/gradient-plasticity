using Pkg
using Plots
Pkg.add("LaTeXStrings")
using LaTeXStrings

include("variables.jl")
using .Variables

function main()
    NUM_TIMESTEPS = Variables.__NUM_TIMESTEPS
    DELTA_X = Variables.__DELTA_X
    DELTA_Y = Variables.__DELTA_Y
    modulus = Variables.__modulus
    OUTPATH_TXT = Variables.__OUTPATH_TXT
    OUTPATH_VTK = Variables.__OUTPATH_VTK
    OUTPATH_GIF = Variables.__OUTPATH_GIF
    OUTPATH_PNG = Variables.__OUTPATH_PNG
    GIF_NAME = Variables.__GIF_NAME
    DOMAIN_X = Variables.__X_DOMAIN
    DOMAIN_Y = Variables.__Y_DOMAIN

    x = [0:DELTA_X:DOMAIN_X]
    y = [0:DELTA_Y:DOMAIN_Y]

    # Open the file and read line by line
    eps11_lines = readlines(OUTPATH_TXT*"eps11.txt")
    eps12_lines = readlines(OUTPATH_TXT*"eps12.txt")
    eps22_lines = readlines(OUTPATH_TXT*"eps22.txt")
    eps33_lines = readlines(OUTPATH_TXT*"eps33.txt")
    eps_e11_lines = readlines(OUTPATH_TXT*"eps_e11.txt")
    eps_e12_lines = readlines(OUTPATH_TXT*"eps_e12.txt")
    eps_e22_lines = readlines(OUTPATH_TXT*"eps_e22.txt")
    eps_e33_lines = readlines(OUTPATH_TXT*"eps_e33.txt")
    sigma11_lines = readlines(OUTPATH_TXT*"sigma11.txt")
    sigma12_lines = readlines(OUTPATH_TXT*"sigma12.txt")
    sigma22_lines = readlines(OUTPATH_TXT*"sigma22.txt")
    sigma33_lines = readlines(OUTPATH_TXT*"sigma33.txt")
    v1_lines = readlines(OUTPATH_TXT*"v1.txt")
    v2_lines = readlines(OUTPATH_TXT*"v2.txt")
    p_lines = readlines(OUTPATH_TXT*"p.txt")
    p_til_lines = readlines(OUTPATH_TXT*"p_til.txt")
    w_lines = readlines(OUTPATH_TXT*"w.txt")
    p1_lines = readlines(OUTPATH_TXT*"p1.txt")
    p2_lines = readlines(OUTPATH_TXT*"p2.txt")

    # Convert each line into an array (assuming space delimits the elements within each line)
    eps11_arrays = [parse.(Float64, split(line)) for line in eps11_lines]
    eps12_arrays = [parse.(Float64, split(line)) for line in eps12_lines]
    eps22_arrays = [parse.(Float64, split(line)) for line in eps22_lines]
    eps_e11_arrays = [parse.(Float64, split(line)) for line in eps_e11_lines]
    eps_e12_arrays = [parse.(Float64, split(line)) for line in eps_e12_lines]
    eps_e22_arrays = [parse.(Float64, split(line)) for line in eps_e22_lines]
    sigma11_arrays = [parse.(Float64, split(line)) for line in sigma11_lines]
    sigma12_arrays = [parse.(Float64, split(line)) for line in sigma12_lines]
    sigma22_arrays = [parse.(Float64, split(line)) for line in sigma22_lines]
    v1_arrays = [parse.(Float64, split(line)) for line in v1_lines]
    v2_arrays = [parse.(Float64, split(line)) for line in v2_lines]
    p_arrays = [parse.(Float64, split(line)) for line in p_lines]
    p_til_arrays = [parse.(Float64, split(line)) for line in p_til_lines]
    w_arrays = [parse.(Float64, split(line)) for line in w_lines]
    p1_arrays = [parse.(Float64, split(line)) for line in p1_lines]
    p2_arrays = [parse.(Float64, split(line)) for line in p2_lines]

    anim = @animate for i in 1:floor(Int,NUM_TIMESTEPS/modulus)
        yeps11 = eps11_arrays[i]
        yeps12 = eps12_arrays[i]
        yeps22 = eps22_arrays[i]
        yeps_e11 = eps_e11_arrays[i]
        yeps_e12 = eps_e12_arrays[i]
        yeps_e22 = eps_e22_arrays[i]
        ysigma11 = sigma11_arrays[i]
        ysigma12 = sigma12_arrays[i]
        ysigma22 = sigma22_arrays[i]
        yv1 = v1_arrays[i]
        yv2 = v2_arrays[i]
        yp = p_arrays[i]
        yp_til = p_til_arrays[i]
        yw = w_arrays[i]
        yp1 = p1_arrays[i]
        yp2 = p2_arrays[i]

        p = plot(x, yeps_e11, seriestype=:line, label=L"\epsilon^e_{11}", color=:black, markersize=0.1, layout=9)
        plot!(p[2], x, yeps_e12, seriestype=:line, label=L"\epsilon^e_{12}", color=:black, markersize=0.1)
        plot!(p[3], x, yeps_e22, seriestype=:line, label=L"\epsilon^e_{22}", color=:black, markersize=0.1)
        plot!(p[4], x, ysigma11, seriestype=:line, label=L"\sigma_{11}", color=:black, markersize=0.1)
        plot!(p[5], x, yv1, seriestype=:line, label=L"v_1", color=:black, markersize=0.1)
        plot!(p[6], x, yv2, seriestype=:line, label=L"v_2", color=:black, markersize=0.1)
        plot!(p[7], x, yp, seriestype=:line, label=L"p", color=:black, markersize=0.1)
        plot!(p[8], x, yp_til, seriestype=:line, label=L"\tilde{p}", color=:black, markersize=0.1)
        plot!(p[9], x, yp1, seriestype=:line, label="p1", color=:black, markersize=0.1)

        if i == floor(Int,NUM_TIMESTEPS/modulus)
            savefig(p, OUTPATH_PNG)
        end

        # plasticity = plot(x, yeps11, seriestype=:scatter, label="ε11", color=:red, markersize=0.5, layout=7)

        # plot!(p[4], x, yeps_e22, seriestype=:scatter, label="ε11", color=:purple, markersize=1)
        # plot!(p[5], x, yeps12, seriestype=:scatter, label="ε12", color=:orange, markersize=1)
        # plot!(p[6], x, yeps22, seriestype=:scatter, label="ε22", color=:yellow, markersize=1)
        # plot!(p[7], x, ysigma11, seriestype=:scatter, label="σ_11", color=:magenta, markersize=1)
        # plot!(p[8], x, ysigma12, seriestype=:scatter, label="σ_12", color=:pink, markersize=1)
        # plot!(p[9], x, ysigma22, seriestype=:scatter, label="σ_22", color=:grey, markersize=1)
        # plot!(p[10], x, yv1, seriestype=:scatter, label="v_1", color=:grey, markersize=1)
        # plot!(p[11], x, yv2, seriestype=:scatter, label="v_2", color=:grey, markersize=1)
        # plot!(p[12], x, yp, seriestype=:scatter, label="p", color=:grey, markersize=1)
    end
    gif(anim, OUTPATH_GIF*GIF_NAME, fps=10)
end

main()