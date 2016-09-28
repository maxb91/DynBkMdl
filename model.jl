using PyPlot

type ModelParams
    l_A
    l_B
end

function pacejka(a)
    B = 0.3#20
    C = 1.25
    mu = 0.234
    m = 1.98
    g = 9.81
    D = mu * m * g/2
    D = D*100

    C_alpha_f = D*sin(C*atan(B*a))
    return C_alpha_f
end


function simModel(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::ModelParams)

    # kinematic bicycle model
    # u[1] = acceleration
    # u[2] = steering angle
    l_A = modelParams.l_A
    l_B = modelParams.l_B

    bta = atan(l_A/(l_A+l_B)*tan(u[2]))

    zNext = z
    zNext[1] = z[1] + dt*(z[4]*cos(z[3] + bta))       # x
    zNext[2] = z[2] + dt*(z[4]*sin(z[3] + bta))     # y
    zNext[3] = z[3] + dt*(z[4]/l_B*sin(bta))        # psi
    zNext[4] = z[4] + dt*(u[1] - 0.63 * z[4]^2 * sign(z[4]))                     # v

    return zNext
end
function simDynModel_exact(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::ModelParams)
    dtn = dt/1000
    t = 0:dtn:dt
    z_final = z
    ang = zeros(2)
    for i=1:length(t)-1
        z_final, ang = simDynModel(z_final,u,dtn,modelParams)
    end
    return z_final, ang
end
function simDynModel(z::Array{Float64},u::Array{Float64},dt::Float64,modelParams::ModelParams)

    zNext::Array{Float64}
    L_f = modelParams.l_A
    L_r = modelParams.l_B
    m   = 1.98
    I_z = 0.24

    a_F = 0
    a_R = 0
    if z[3] != 0
        a_F     = atan((z[4] + L_f*z[6])/z[3]) - u[2]
        a_R     = atan((z[4] - L_r*z[6])/z[3])
    end

    C_alpha_f = pacejka(a_F)
    C_alpha_r = pacejka(a_R)

    FyF = -C_alpha_f# * a_F
    FyR = -C_alpha_r# * a_R

    zNext = z
    # compute next state
    zNext[1]        = zNext[1]       + dt * (cos(z[5])*z[3] - sin(z[5])*z[4])
    zNext[2]        = zNext[2]       + dt * (sin(z[5])*z[3] + cos(z[5])*z[4])
    zNext[3]        = zNext[3]       + dt * (u[1] + z[4]*z[6] - 0.63*z[3]^2*sign(z[3]))
    zNext[4]        = zNext[4]       + dt * (2/m*(FyF*cos(u[2]) + FyR) - z[6]*z[3])
    zNext[5]        = zNext[5]       + dt * (z[6])
    zNext[6]        = zNext[6]       + dt * (2/I_z*(L_f*FyF - L_r*FyR))

    return zNext, [a_F a_R]
end

function simDynModel_ODE(t::Float64,z::Array{Float64},u::Array{Float64},modelParams::ModelParams)

    zNext::Array{Float64}
    L_f = modelParams.l_A
    L_r = modelParams.l_B
    m   = 1.98
    I_z = 0.24

    a_F = 0
    a_R = 0
    if z[3] != 0
        a_F     = atan((z[4] + L_f*z[6])/z[3]) - u[2]
        a_R     = atan((z[4] - L_r*z[6])/z[3])
    end
    println("a_F = $a_F, a_R = $a_R")

    C_alpha_f = 10
    C_alpha_r = 10

    #C_alpha_f = pacejka(a_F)
    #C_alpha_r = pacejka(a_R)

    println(C_alpha_f)

    FyF = -C_alpha_f * a_F
    FyR = -C_alpha_r * a_R


    zNext = z
    # compute next state
    zNext[1]        = (cos(z[5])*z[3] - sin(z[5])*z[4])
    zNext[2]        = (sin(z[5])*z[3] + cos(z[5])*z[4])
    zNext[3]        = (u[1] + z[4]*z[6])
    zNext[4]        = (1/m*(FyF*cos(u[2]) + FyR) - z[6]*z[3])
    zNext[5]        = (z[6])
    zNext[6]        = (1/I_z*(L_f*FyF - L_r*FyR))

    #zNext = zNext + 0*randn(1,4)*0.001
    return zNext, [a_F a_R]
end

# modelParams = ModelParams(0.125,0.125)
# z0 = zeros(6)
# u = zeros(2)
# t,z = ode45((t,z)->simDynModel_ODE(t,z,u,modelParams),z0,1)


function startSim()
    modelParams = ModelParams(0.125,0.125)

    dt      = 0.01
    t       = 0:dt:3
    z       = zeros(size(t,1),6)
    z_kin   = zeros(size(t,1),6)
    #z_ode   = zeros(size(t,1),6)
    u       = zeros(size(t,1),2)
    ang     = zeros(size(t,1),2)
    #u[50:100,1] = 1
    #u[100:150,2] = 0.5

    for i=2:size(t,1)
        if t[i] < 1
            u[i,1] = 1
        elseif t[i] < 2
            u[i,1] = 0.0
            u[i,2] = 0.5
        #elseif t[i] < 4
            #u[i,1] = 1
            #u[i,2] = 0.5
        else
            u[i,:] = [0 0]
        end
        z[i,:],ang[i,:]      = simDynModel_exact(z[i-1,:],u[i-1,:],dt,modelParams)
        #t_sim,z_sim = ode45((t_sim,z_sim)->simDynModel_ODE(t_sim,z_sim,u,modelParams),z_ode[i-1,:],0:dt)
        #z_ode[i,:]  = z_sim[end]
        z_kin[i,:]  = simModel(z_kin[i-1,:],u[i-1,:],dt,modelParams)
    end

    plot(t,z,t,z_kin,"--")
    grid()
    legend(["x","y","v_x","v_y","psi","psi_dot","x","y","psi","v"])
    figure()
    plot(z[:,1],z[:,2],z_kin[:,1],z_kin[:,2])
    legend(["Dynamic","Kinematic"])
    grid()

    figure()
    plot(t,u)
    legend(["a","d_f"])
    grid()

    figure()
    plot(t,ang)
    legend(["a_F","a_R"])
    grid()
end

function showPacejka()
    a = -20:.01:20
    f = pacejka(a)
    plot(a,f)
    grid(1)
end
