include("variables.jl")
include("ics_bcs_v6.jl")
using DifferentialEquations
using WriteVTK
using Plots
using Printf
using .Variables: eps_e11, eps_e12, eps_e22, eps_e33, eps11, eps12, eps22, eps33, u1, u2, p, p_til, w, p1, p2


### DEFINITIONS
#----------------------------------------------------------------------------------#
# U = (eps_e11, eps_e12, eps_e22, eps_e33, eps11, eps12, eps22, eps33, u1, u2, p)
# F_x = (-v1, -1/2*v2, 0, -v1, -1/2*v2, 0, -1/p_0*s11, -1/p_0*s12, 0)
# F_y = (0, -1/2v1, -v2, 0, -1/2*v1, -v2, -1/p_0*s12, -1/p_0*s22, 0)
# S = (...)
#----------------------------------------------------------------------------------#


function compute_flux_riemann(U, Fx, Fy, Fx_shift, Fy_shift, v1, v2, sigma11, sigma12, sigma22, sigma33, eta, tr_eps, Mu, lambda, RHO, lc, beta, NUM_X_CELLS, NUM_Y_CELLS, p0, oneD)
    ### Flux term variables
    tr_eps[:,:] = U[:,:,eps_e11] + U[:,:,eps_e22] + U[:,:,eps_e33]

    if oneD          # -- SWITCHED TO 1D --
        sigma11[:,:] = @. 2*Mu*U[:,:,eps_e11] + lambda*U[:,:,eps_e11]
    else
        sigma11[:,:] = @. 2*Mu*U[:,:,eps_e11] + lambda*tr_eps[:,:]
    end
    sigma12[:,:] = @. 2*Mu*U[:,:,eps_e12]
    sigma22[:,:] = @. 2*Mu*U[:,:,eps_e22] + lambda*tr_eps[:,:]
    sigma33[:,:] = @. 2*Mu*U[:,:,eps_e33] + lambda*tr_eps[:,:]
    v1[:,:] = @. U[:,:,u1]/RHO
    v2[:,:] = @. U[:,:,u2]/RHO

    ### Compute Flux
    Fx[:,:,eps_e11]= @. -v1[:,:]
    Fx[:,:,eps_e12]= @. -0.5*v2[:,:]
    Fx[:,:,eps_e22].= 0
    Fx[:,:,eps_e33].= 0
    Fx[:,:,eps11]= @. -v1[:,:]
    Fx[:,:,eps12]= @. -0.5*v2[:,:]
    Fx[:,:,eps22].= 0
    Fx[:,:,eps33].= 0
    Fx[:,:,u1]= @. -sigma11[:,:]
    Fx[:,:,u2]= @. -sigma12[:,:]
    Fx[:,:,p].= 0
    Fx[:,:,p_til].= 0
    Fx[:,:,w].= lc/beta * U[:,:,p1]
    Fx[:,:,p1].= -U[:,:,w]
    Fx[:,:,p2].= 0

    Fy[:,:,eps_e11].= 0
    Fy[:,:,eps_e12]= @. -0.5*v1[:,:]
    Fy[:,:,eps_e22]= @. -v2[:,:]
    Fy[:,:,eps_e33].= 0
    Fy[:,:,eps11].= 0
    Fy[:,:,eps12]= @. -0.5*v1[:,:]
    Fy[:,:,eps22]= @. -v2[:,:]
    Fy[:,:,eps33].= 0
    Fy[:,:,u1]= @. -sigma12[:,:]
    Fy[:,:,u2]= @. -sigma22[:,:]
    Fy[:,:,p].= 0
    Fy[:,:,p_til].= 0
    Fy[:,:,w].= lc/beta * U[:,:,p2]
    Fy[:,:,p1].= 0
    Fy[:,:,p2].= -U[:,:,w]

    for i in 1:NUM_X_CELLS+1
        for j in 1:NUM_Y_CELLS+1
            Fx_shift[i,j,:] = @. 0.5 * ((Fx[i+1,j,:]+Fx[i,j,:]) -  eta[i,j] * (U[i+1,j,:] - U[i,j,:])) # NTS: eta mult is ambiguous but if compiler is good it should be fine.
            if oneD          # -- SWITCHED TO 1D --
                Fy_shift[i,j,:] .= 0
            else
                Fy_shift[i,j,:] .= @. 0.5 * ((Fy[i,j+1,:]+Fy[i,j,:]) -  eta[i,j] * (U[i,j+1,:] - U[i,j,:]))
            end
        end
    end
end


function evolve_U_godunov(U, Fx_shift, Fy_shift, eta, cL, cR, RHO, NUM_X_CELLS, NUM_Y_CELLS, DELTA_T, DELTA_X, DELTA_Y)
    ### Conservative scheme variables
    eta[:,:] = @. max(cL[:,:], cR[:,:]) # needs fixing.
    for i in 2:NUM_X_CELLS+1
        for j in 2:NUM_Y_CELLS+1
            U[i,j,:]= @. U[i,j,:] - DELTA_T*(1/DELTA_X * (Fx_shift[i,j,:] - Fx_shift[i-1,j,:]) + 1/DELTA_Y * (Fy_shift[i,j,:] - Fy_shift[i,j-1,:]))
        end
    end
end


function solve_source_term(U, Mu, lambda, RHO, sigma_y, kappa, k, beta, YOUNG_MOD_E, NUM_X_CELLS, NUM_Y_CELLS, DELTA_T, nsteps, oneD; solver=Euler())
    for i in 1:NUM_X_CELLS+2
        for j in 1:NUM_Y_CELLS+2
            u0 = vec(U[i,j,:])
            tspan = (0.0, 1*DELTA_T)
            p = (Mu, lambda, RHO, sigma_y, kappa, k, YOUNG_MOD_E, beta, oneD)
            prob = ODEProblem(f!, u0, tspan, p)
            sol = solve(prob, solver, dt=DELTA_T/nsteps, adaptive=true)
            # print(sol.u)
            U[i,j,:] = sol.u[nsteps+1]
        end
    end
end


function f!(du, u, prm, t)

    Mu = prm[1]
    lambda = prm[2]
    RHO = prm[3]
    sigma_y = prm[4]
    kappa = prm[5]
    k = prm[6]
    YOUNG_MOD_E = prm[7]
    beta = prm[8]
    oneD = prm[9]

    eps_e11_val = u[eps_e11]
    eps_e12_val = u[eps_e12]
    eps_e22_val = u[eps_e22]
    eps_e33_val = u[eps_e33]
    eps11_val   = u[eps11]
    eps12_val   = u[eps12]
    eps22_val   = u[eps22]
    eps33_val   = u[eps33]
    u1_val      = u[u1]
    u2_val      = u[u2]
    p_val       = u[p]
    p_til_val   = u[p_til]
    w_val       = u[w]
    p1_val      = u[p1]
    p2_val      = u[p2]

    tr_eps = eps_e11_val + eps_e22_val + eps_e33_val
    
    if oneD          # -- SWITCHED TO 1D --
        sigma11 = 2 * Mu * eps_e11_val + lambda * eps_e11_val
    else
        sigma11 = 2 * Mu * eps_e11_val + lambda * tr_eps
    end
    sigma12 = 2 * Mu * eps_e12_val
    sigma22 = 2 * Mu * eps_e22_val + lambda * tr_eps
    sigma33 = 2 * Mu * eps_e33_val + lambda * tr_eps
    tr_sigma = sigma11 + sigma22 + sigma33
    one_third_tr_sigma = (1/3) * tr_sigma
    if oneD          # -- SWITCHED TO 1D --
        sigma11_prime = sigma11
    else
        sigma11_prime = sigma11 - one_third_tr_sigma
    end
    sigma12_prime = sigma12
    sigma22_prime = sigma22 - one_third_tr_sigma
    sigma33_prime = sigma33 - one_third_tr_sigma
    
    if Variables.__oneD          # -- SWITCHED TO 1D --
        sigmaP_dot_sigmaP = sigma11_prime^2
    else
        sigmaP_dot_sigmaP = sigma11_prime^2 + 2*sigma12_prime^2 + sigma22_prime^2 + sigma33_prime^2
    end
    norm_sigmaP = sqrt(1.5 * sigmaP_dot_sigmaP)

    if Variables.__SOFTENING
        R_val = max(YOUNG_MOD_E / 20 * (p_val - 5000 * p_val^2),-0.95 * sigma_y)
    else
        R_val = YOUNG_MOD_E / 20 * p_val
    end

    phi = norm_sigmaP - sigma_y - R_val + kappa * (p_til_val - p_val)
    
    if phi > 0
        lambda_dot = phi / (k * sigma_y)
        src_eps_e11 = -lambda_dot * sigma11_prime / norm_sigmaP
        src_eps_e12 = -lambda_dot * sigma12_prime / norm_sigmaP
        src_eps_e22 = -lambda_dot * sigma22_prime / norm_sigmaP
        src_eps_e33 = -lambda_dot * sigma33_prime / norm_sigmaP
    else
        lambda_dot = 0.0
        src_eps_e11 = 0.0
        src_eps_e12 = 0.0
        src_eps_e22 = 0.0
        src_eps_e33 = 0.0
    end

    du[eps_e11] = src_eps_e11
    du[eps_e12] = src_eps_e12
    du[eps_e22] = src_eps_e22
    du[eps_e33] = src_eps_e33
    du[eps11]   = 0.0
    du[eps12]   = 0.0
    du[eps22]   = 0.0
    du[eps33]   = 0.0
    du[u1]      = 0.0
    du[u2]      = 0.0
    du[p]       = lambda_dot
    du[p_til]   = w_val
    du[w]       = kappa * (p_val - p_til_val) * 1/beta
    du[p1]      = 0.0
    du[p2]      = 0.0
end


function compute_stress_strain_curve(U, half_NXC, half_NYC, Mu, lambda, sigma11, sigma12, sigma22, sigma33, stress, strain, t)
    # not yet configured for 1D problem
    offset=Variables.__OFFSET
    t_tr_sigma = sigma11[half_NXC-offset, half_NYC-offset] + sigma22[half_NXC-offset, half_NYC-offset] + sigma33[half_NXC-offset, half_NYC-offset]
    t_one_third_tr_sigma = (1/3) * t_tr_sigma
    t_sigma11_prime = sigma11[half_NXC-offset, half_NYC-offset] - t_one_third_tr_sigma
    t_sigma12_prime = sigma12[half_NXC-offset, half_NYC-offset]
    t_sigma22_prime = sigma11[half_NXC-offset, half_NYC-offset] - t_one_third_tr_sigma
    t_sigma33_prime = sigma11[half_NXC-offset, half_NYC-offset] - t_one_third_tr_sigma
    t_sigmaP_dot_sigmaP = t_sigma11_prime^2 + 2 * t_sigma12_prime^2 + t_sigma22_prime^2 + t_sigma33_prime^2
    t_norm_sigmaP = sqrt(1.5 * t_sigmaP_dot_sigmaP)

    stress[t+1] = t_norm_sigmaP
    strain[t+1] = (U[half_NXC-offset, half_NYC-offset, eps11] + U[half_NXC-offset, half_NYC-offset, eps22] + U[half_NXC-offset, half_NYC-offset, eps33])
end

function prepare_output_folders()
    # creates folder where simulation output is set to be stored
    out_folder = Variables.__OUTPATH
    if !isdir(out_folder)
        mkdir(out_folder)
        println("Folder '$out_folder' created.")
    else
        println("Folder '$out_folder' already exists.")
    end

    txt_folder = Variables.__OUTPATH_TXT
    if !isdir(txt_folder)
        mkdir(txt_folder)
        println("Folder '$txt_folder' created.")
    else
        println("Folder '$txt_folder' already exists.")
    end

    vtk_folder = Variables.__OUTPATH_VTK
    if !isdir(vtk_folder)
        mkdir(vtk_folder)
        println("Folder '$vtk_folder' created.")
    else
        # Delete all .vti files in the vtk directory
        for file in readdir(Variables.__OUTPATH_VTK)
            if endswith(file, ".vti")
                rm(Variables.__OUTPATH_VTK * file)
            end
        end
        println("Folder '$vtk_folder' already exists.")
    end
end






function main()
    ### DEFINITIONS
    #----------------------------------------------------------------------------------#
    # U = (eps_e11, eps_e12, eps_e22, eps_e33, eps11, eps12, eps22, eps33, u1, u2, p)
    # F_x = (-v1, -1/2*v2, 0, -v1, -1/2*v2, 0, -1/p_0*s11, -1/p_0*s12, 0)
    # F_y = (0, -1/2v1, -v2, 0, -1/2*v1, -v2, -1/p_0*s12, -1/p_0*s22, 0)
    #----------------------------------------------------------------------------------#

    ### CONSTANTS defined in variables.jl
    NUM_EQUATIONS = 15
    OUTPATH_TXT = Variables.__OUTPATH_TXT
    OUTPATH_VTK = Variables.__OUTPATH_VTK
    VTK_NAME = Variables.__VTK_NAME
    OUTPATH_STRESS_STRAIN = Variables.__OUTPATH_STRESS_STRAIN
    NUM_X_CELLS = Variables.__NUM_X_CELLS
    NUM_Y_CELLS = Variables.__NUM_Y_CELLS
    half_NXC = Variables.__half_NXC
    half_NYC = Variables.__half_NYC
    DOMAIN_X = Variables.__X_DOMAIN
    DOMAIN_Y = Variables.__Y_DOMAIN
    RHO = Variables.__RHO
    YOUNG_MOD_E = Variables.__YOUNG_MOD_E
    MAX_WAVESPEED = Variables.__MAX_WAVESPEED
    CFL = Variables.__CFL
    DELTA_X = Variables.__DELTA_X
    DELTA_Y = Variables.__DELTA_Y
    DELTA_T = Variables.__DELTA_T
    NUM_TIMESTEPS = Variables.__NUM_TIMESTEPS
    oneD = Variables.__oneD
    modulus = Variables.__modulus
    INIT_CONDITIONS = Variables.__INIT_CONDITIONS
    BOUND_CONDITIONS = Variables.__BOUND_CONDITIONS

    ### INITIALIZE conservative scheme
    U = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2, NUM_EQUATIONS)
    Fx = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2, NUM_EQUATIONS)
    Fy = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2, NUM_EQUATIONS)
    Fx_shift = zeros(Float32, NUM_X_CELLS+1, NUM_Y_CELLS+1, NUM_EQUATIONS)
    Fy_shift = zeros(Float32, NUM_X_CELLS+1, NUM_Y_CELLS+1, NUM_EQUATIONS)
    x = 0:DELTA_X:DOMAIN_X
    y = 0:DELTA_Y:DOMAIN_Y

    ### VARIABLES and CONSTANTS used to define flux F and F_shift
    Nu = 0.3
    Mu = YOUNG_MOD_E/(2*(1+Nu))
    lambda = YOUNG_MOD_E*Nu/(1+Nu)/(1-2*Nu)
    cR = [(YOUNG_MOD_E*(1) / RHO)^0.5 for i in 1:NUM_X_CELLS+1, j in 1:NUM_Y_CELLS+1]
    cL = [(YOUNG_MOD_E*(1) / RHO)^0.5 for i in 1:NUM_X_CELLS+1, j in 1:NUM_Y_CELLS+1]
    sigma11 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    sigma12 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    sigma22 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    sigma33 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    v1 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    v2 = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    eta = zeros(Float32, NUM_X_CELLS+1, NUM_Y_CELLS+1)
    tr_eps = zeros(Float32, NUM_X_CELLS+2, NUM_Y_CELLS+2)
    strain = zeros(Float32, NUM_TIMESTEPS)
    stress = zeros(Float32, NUM_TIMESTEPS)
    
    ### CONSTANTS used for the source term S
    sigma_y=2.00e6
    p0 = sigma_y/YOUNG_MOD_E
    lc = 0.01
    beta_0 = 0.9
    beta = beta_0 * RHO/YOUNG_MOD_E * lc
    kappa = 10



    k=Variables.__k
    solver=ImplicitEuler()
    nsteps=1

    #----- set the extra variables we are plotting -----#
    v1[:,:] .= U[:,:,u1]/RHO
    v2[:,:] .= U[:,:,u2]/RHO






    prepare_output_folders()

    open(OUTPATH_TXT * "eps11.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps12.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps22.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps33.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps_e11.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps_e12.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps_e22.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "eps_e33.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "sigma11.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "sigma12.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "sigma22.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "sigma33.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "v1.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "v2.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "p.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "p_til.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "w.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "p1.txt", "w") do file write(file, "") end
    open(OUTPATH_TXT * "p2.txt", "w") do file write(file, "") end
 

    open(OUTPATH_TXT * "eps11.txt", "a") do eps11_file
    open(OUTPATH_TXT * "eps12.txt", "a") do eps12_file
    open(OUTPATH_TXT * "eps22.txt", "a") do eps22_file
    open(OUTPATH_TXT * "eps33.txt", "a") do eps33_file
    open(OUTPATH_TXT * "eps_e11.txt", "a") do eps_e11_file
    open(OUTPATH_TXT * "eps_e12.txt", "a") do eps_e12_file
    open(OUTPATH_TXT * "eps_e22.txt", "a") do eps_e22_file
    open(OUTPATH_TXT * "eps_e33.txt", "a") do eps_e33_file
    open(OUTPATH_TXT * "sigma11.txt", "a") do sigma11_file
    open(OUTPATH_TXT * "sigma12.txt", "a") do sigma12_file
    open(OUTPATH_TXT * "sigma22.txt", "a") do sigma22_file
    open(OUTPATH_TXT * "sigma33.txt", "a") do sigma33_file
    open(OUTPATH_TXT * "v1.txt", "a") do v1_file
    open(OUTPATH_TXT * "v2.txt", "a") do v2_file
    open(OUTPATH_TXT * "p.txt", "a") do p_file
    open(OUTPATH_TXT * "p_til.txt", "a") do p_til_file
    open(OUTPATH_TXT * "w.txt", "a") do w_file
    open(OUTPATH_TXT * "p1.txt", "a") do p1_file
    open(OUTPATH_TXT * "p2.txt", "a") do p2_file

    xs = 0:DELTA_X:DOMAIN_X
    ys = 0:DELTA_Y:DOMAIN_Y

    @time begin
    initial_conditions(U, xs, RHO, INIT_CONDITIONS)
    for t in 0:NUM_TIMESTEPS-1
        compute_stress_strain_curve(U, half_NXC, half_NYC, Mu, lambda, sigma11, sigma12, sigma22, sigma33, stress, strain, t)
        boundary_conditions(U, BOUND_CONDITIONS)
        if t % modulus == 0
            ### write data to files         ### gif command : ffmpeg -framerate 10 -i frame%04d.png output.gif
            # joined = join(U[1:NUM_X_CELLS+2, half_NYC, eps11], " ") # USE A FOR LOOP AND THIS?
            println(eps11_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps11], " "))
            println(eps12_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps12], " "))
            println(eps22_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps22], " "))
            println(eps33_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps33], " "))
            println(eps_e11_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps_e11], " "))
            println(eps_e12_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps_e12], " "))
            println(eps_e22_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps_e22], " "))
            println(eps_e33_file, join(U[1:NUM_X_CELLS+2, half_NYC, eps_e33], " "))
            println(sigma11_file, join(sigma11[1:NUM_X_CELLS+2, half_NYC], " "))
            println(sigma12_file, join(sigma12[1:NUM_X_CELLS+2, half_NYC], " "))
            println(sigma22_file, join(sigma22[1:NUM_X_CELLS+2, half_NYC], " "))
            println(sigma33_file, join(sigma33[1:NUM_X_CELLS+2, half_NYC], " "))
            println(v1_file, join(v1[1:NUM_X_CELLS+2, half_NYC], " "))
            println(v2_file, join(v2[1:NUM_X_CELLS+2, half_NYC], " "))
            println(p_file, join(U[1:NUM_X_CELLS+2, half_NYC, p], " "))
            println(p_til_file, join(U[1:NUM_X_CELLS+2, half_NYC, p_til], " "))
            println(w_file, join(U[1:NUM_X_CELLS+2, half_NYC, w], " "))
            println(p1_file, join(U[1:NUM_X_CELLS+2, half_NYC, p1], " "))
            println(p2_file, join(U[1:NUM_X_CELLS+2, half_NYC, p2], " "))
            
            # println(eps11_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps11], " "))
            # println(eps12_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps12], " "))
            # println(eps22_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps22], " "))
            # println(eps33_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps22], " "))
            # println(eps_e11_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps_e11], " "))
            # println(eps_e12_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps_e12], " "))
            # println(eps_e22_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps_e22], " "))
            # println(eps_e33_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, eps_e22], " "))
            # println(sigma11_file, join(sigma11[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(sigma12_file, join(sigma12[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(sigma22_file, join(sigma22[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(sigma33_file, join(sigma22[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(v1_file, join(v1[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(v2_file, join(v2[half_NXC, 1:NUM_Y_CELLS+2], " "))
            # println(p_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, p], " "))
            # println(p_til_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, p_til], " "))
            # println(w_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, w], " "))
            # println(p1_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, p1], " "))
            # println(p2_file, join(U[half_NXC, 1:NUM_Y_CELLS+2, p2], " "))


            t_str = @sprintf("%04d", t)
            vtk_grid("$OUTPATH_VTK$VTK_NAME$t_str", xs, ys) do vtk # FIX THIS NOW!!!!!!!!!!
                vtk["eps11"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps11]
                vtk["eps12"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps12]
                vtk["eps22"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps22]
                vtk["eps33"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps33]
                vtk["eps_e11"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps_e11]
                vtk["eps_e12"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps_e12]
                vtk["eps_e22"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps_e22]
                vtk["eps_e33"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,eps_e33]
                vtk["sigma11"] = sigma11[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["sigma12"] = sigma12[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["sigma22"] = sigma22[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["sigma33"] = sigma33[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["v1"] = v1[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["v2"] = v2[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1]
                vtk["p"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,p]
                vtk["p_til"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,p_til]
                vtk["w"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,w]
                vtk["p1"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,p1]
                vtk["p2"] = U[2:NUM_X_CELLS+1,2:NUM_Y_CELLS+1,p2]
            end
        end
        compute_flux_riemann(U, Fx, Fy, Fx_shift, Fy_shift, v1, v2, sigma11, sigma12, sigma22, sigma33, eta, tr_eps, Mu, lambda, RHO, lc, beta, NUM_X_CELLS, NUM_Y_CELLS, p0, oneD)
        evolve_U_godunov(U, Fx_shift, Fy_shift, eta, cL, cR, RHO, NUM_X_CELLS, NUM_Y_CELLS, DELTA_T, DELTA_X, DELTA_Y)
        @time begin
        solve_source_term(U, Mu, lambda, RHO, sigma_y, kappa, k, beta, YOUNG_MOD_E, NUM_X_CELLS, NUM_Y_CELLS, DELTA_T, nsteps, oneD; solver=solver)
        print(t)
        end
        # Blowup test
        if abs(U[half_NXC, half_NYC, eps11]) >= 10 || isnan(U[half_NXC, half_NYC, eps11])
            print("blowup at timestep $t, with value e11 = ")
            print(U[half_NXC, half_NYC, eps11])
            print("\n")
            break
        end
    end
    end
end end end end end end end end end end end end end end end end end end end
   
    ### Plots stress strain curve
    plot(
        strain, stress,
        xlabel="Strain",
        ylabel="Stress",
        seriestype=:line,
        # title="Stress-Strain Curve",
        legend=false,
        color=:black,
        grid=true
    )
    savefig(OUTPATH_STRESS_STRAIN)

end

main()