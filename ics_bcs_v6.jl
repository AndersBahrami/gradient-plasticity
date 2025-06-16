using .Variables: eps_e11, eps_e12, eps_e22, eps_e33, eps11, eps12, eps22, eps33, u1, u2, p, p_til, w, p1, p2

function initial_conditions(U, xs, RHO, INIT_CONDITIONS)
    for i in range(1, length(xs))
        x = xs[i]
        U[i, :, u1] .= u0(x, xs, RHO, INIT_CONDITIONS)
    end
end

function u0(x, xs, RHO, INIT_CONDITIONS)
    return 100. * RHO * (xs[end]-x) / (xs[end])

    #----------- INITIAL CONDITIONS for U ----------#
    ### INITIAL CONDITION: X impact
    # for j in 1:NUM_Y_CELLS+2
    #     for i in 1:half_NXC+1
    #         U[i,j,u1] = 1. *RHO
    #         # U[i,j,p1] = 1. *RHO
    #     end
    #     for i in half_NXC+2:NUM_X_CELLS+2
    #         U[i,j,u1] = -1. *RHO
    #         # U[i,j,p1] = -1. *RHO
    #     end
    # end

    ### INITIAL CONDITION: Y impact
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:half_NYC+1
    #         U[i,j,u2] = 1. *RHO
    #         # U[i,j,eps_e11] = 0.01
    #     end
    #     for j in half_NYC+2:NUM_Y_CELLS+2
    #         U[i,j,u2] = -1. *RHO
    #         # U[i,j,eps_e11] = 0
    #     end
    # end


    ### INITIAL CONDITION: X pull test
    # for j in 1:NUM_Y_CELLS+2
    #     for i in 1:half_NXC+1
    #         U[i,j,u1] = -1. *RHO
    #     end
    #     for i in half_NXC+2:NUM_X_CELLS+2
    #         U[i,j,u1] = 1. *RHO
    #     end
    # end


    ### INITIAL CONDITION: Y pull test
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:half_NYC+1
    #         U[i,j,u2] = -1. *RHO
    #     end
    #     for j in half_NYC+2:NUM_Y_CELLS+2
    #         U[i,j,u2] = 1. *RHO
    #     end
    # end


    ### INITIAL CONDITION: X shear test
    # for j in 1:NUM_Y_CELLS+2
    #     for i in 1:half_NXC+1
    #         U[i,j,u2] = -1. *RHO
    #     end
    #     for i in half_NXC+2:NUM_X_CELLS+2
    #         U[i,j,u2] = 1. *RHO
    #     end
    # end


    ### INITIAL CONDITION: Y shear test
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:half_NYC+1
    #         U[i,j,u1] = 1. *RHO
    #     end
    #     for j in half_NYC+2:NUM_Y_CELLS+2
    #         U[i,j,u1] = -1. *RHO
    #     end
    # end


    ### INITIAL CONDITION: Left moving, right fixed wall, linear gradient velocity
    #=
        Use long bar, pay attention to the number of pieces broken into. 
        Try for all mesh size cases n=1,2,3,4,5 and compare least squares with n=8 (i.e. 100*2^8 cells).
        Try all of these for both 1D and 2D cases and also for both stress strain curves.
        Check notes for more specifications on the 2D case.
    =#
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:NUM_Y_CELLS+2
    #         U[i,j,u1] = 1. * RHO * (i-1) / (NUM_X_CELLS+1)
    #         U[i,j,u2] = 0
    #     end
    # end


    ### INITIAL CONDITION: Bomb
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:NUM_Y_CELLS+2
    #         if abs(i-half_NXC-1) <= 10 && abs(j-half_NYC-1) <= 10
    #             U[i,j,eps11] = 0.01
    #             U[i,j,eps_e11] = 0.01
    #             U[i,j,eps22] = 0.01
    #             U[i,j,eps_e22] = 0.01
    #             # U[i,j,w] = 0
    #         end
    #     end
    # end

    ### INITIAL CONDITION: LEFT/RIGHT STRAIN (like in J.A. Shaw's paper)
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:NUM_Y_CELLS+2
    #         if 105-5*i > j && 75-5*i < j
    #             U[i,j,eps11] = -0.001
    #             U[i,j,eps_e11] = -0.001
    #         end
    #     end
    # end
end








function boundary_conditions(U, BOUND_CONDITIONS)
    ### BOUNDARY CONDITIONS
    # Pass through
    U[1,:,:]=U[2,:,:]
    U[:,1,:]=U[:,2,:]
    U[end,:,:] = U[end-1,:,:]
    U[:,end,:] = U[:,end-1,:]

    ## Double pull
    # U[2,:,eps11].=-0.001 SHOULD BE 1, NOT 2,?
    # U[2,:,eps_e11].=-0.001
    # U[NUM_X_CELLS+1,:,eps11].= 0.001
    # U[NUM_X_CELLS+1,:,eps_e11].= 0.001

    ## Diagonal strip
    # for i in 1:NUM_X_CELLS+2
    #     for j in 1:NUM_Y_CELLS+2
    #         if 105-5*i > j && 75-5*i < j
    #             U[i,j,eps11] = -0.001
    #             U[i,j,eps_e11] = -0.001
    #         end
    #     end
    # end

    ## Reflective
    # U[1,:,u1]= -U[2,:,u1]
    # U[:,1,u2]= -U[:,2,u2]
    # U[NUM_X_CELLS+2,:,u1] = -U[NUM_X_CELLS+1,:,u1]
    # U[:,NUM_Y_CELLS+2,u2] = -U[:,NUM_Y_CELLS+1,u2]

    ## Other Reflective
    # U[1,:,eps11].= 0.001#-U[2,:,eps11]
    # U[1,:,eps_e11].= 0.001#-U[2,:,eps_e11]
    # U[:,1,eps22]= -U[:,2,eps22]
    # U[:,1,eps_e22]= -U[:,2,eps_e22]
    # U[NUM_X_CELLS+2,:,eps11] = -U[NUM_X_CELLS+1,:,eps11]
    # U[NUM_X_CELLS+2,:,eps_e11] = -U[NUM_X_CELLS+1,:,eps_e11]
    # U[:,NUM_Y_CELLS+2,eps22] = -U[:,NUM_Y_CELLS+1,eps22]
    # U[:,NUM_Y_CELLS+2,eps_e22] = -U[:,NUM_Y_CELLS+1,eps_e22]

    ## Left wall, right pull, weird boundary
    #=
        # NOTE: Change sigma12, sigma22 to 0 in the source term and maybe elsewhere?
        Use long bar, pay attention to the number of pieces broken into. 
        Try for all mesh size cases n=1,2,3,4,5 and compare least squares with n=8 (i.e. 100*2^8 cells).
        Try all of these for both 1D and 2D cases and also for both stress strain curves.
        Check notes for more specifications on the 2D case.
    =#
    U[end,:,u1].=0
    U[end,:,u2].=0
    # U[NUM_X_CELLS+2,:,eps11].= 0.001
    # U[NUM_X_CELLS+2,:,eps_e11].= 0.001
    # U[:,1,1:8].=0
    # U[:,NUM_Y_CELLS,1:eps22].=0
    # U[:,NUM_Y_CELLS,1:eps_e22].=0
    # U[:,NUM_Y_CELLS,1:eps12].=0
    # U[:,NUM_Y_CELLS,1:eps_e12].=0
end