include("../find_meteor.jl")
# using Plots

using JSON
using Meteor.Utility
using Meteor.ExactDiag


function bosonic_thermal_state(d::Int, n::Real)
    q = n / (n+1)
    r = [q^j for j in 0:d-1]
    r ./= sum(r)
    v = zeros(d, d)
    for i in 1:d
        v[i, i] = r[i]
    end
    return v
end

function four_level_eit_cooling(d, Delta_g, Delta_r, Delta_t, nu, eta, costheta, omega_L1, omega_L2, gamma_g, gamma_r, nbar, t)
    ps = nlevel(["g", "r", "e", "t"])
    pb = boson(d=d)    
    xop = pb["a"] + pb["adag"]
    
    model = build_model([ps, pb], isunitary=false)
    add_unitary!(model, (1,), ("e<-e",), coeff=-Delta_g)
    add_unitary!(model, (1,), ("r<-r",), coeff=-Delta_g+Delta_r)
    add_unitary!(model, (1,), ("t<-t",), coeff=-Delta_g+Delta_r-Delta_t)
    add_unitary!(model, (2,), ("n",), coeff=nu)
    
    exp_xop = exp((im*eta*costheta)*xop)
    add_unitary!(model, (1, 2), ("g<-e", exp_xop), coeff=omega_L1/2)
    add_unitary!(model, (1, 2), ("e<-g", exp_xop'), coeff=omega_L1/2)
    add_unitary!(model, (1, 2), ("r<-t", exp_xop), coeff=omega_L1/2)
    add_unitary!(model, (1, 2), ("t<-r", exp_xop'), coeff=omega_L1/2)    

    add_unitary!(model, (1, 2), ("r<-e", exp_xop'), coeff=omega_L2/2)
    add_unitary!(model, (1, 2), ("e<-r", exp_xop), coeff=omega_L2/2)
    
    dx = 0.015
    for x in -1:dx:1
        c = dx * 0.75 * (1 + x^2)
        jump = sqrt(c / 2) * exp(-im*eta*x*xop)
        add_dissipation!(model, (1, 2), ("g<-e", jump), coeff=gamma_g/2)
    end   
    for x in -1:dx:1
        c = dx * 0.75 * (1 + x^2)
        jump = sqrt(c / 2) * exp(-im*eta*x*xop)
        add_dissipation!(model, (1, 2), ("r<-t", jump), coeff=gamma_r/2)
    end 

    for x in -1:dx:1
        c = dx * 1.5 * (1 - x^2)
        jump = sqrt(c / 2) * exp(-im*eta*x*xop)
        add_dissipation!(model, (1, 2), ("r<-e", jump), coeff=gamma_r/2)
    end
    for x in -1:dx:1
        c = dx * 1.5 * (1 - x^2)
        jump = sqrt(c / 2) * exp(-im*eta*x*xop)
        add_dissipation!(model, (1, 2), ("g<-t", jump), coeff=gamma_g/2)
    end
    # add_dissipation!(model, (1,), ("g<-e",), coeff=gamma_g/2)
    # add_dissipation!(model, (1,), ("r<-t",), coeff=gamma_r/2)
    # add_dissipation!(model, (1,), ("r<-e",), coeff=gamma_r/2)
    # add_dissipation!(model, (1,), ("g<-t",), coeff=gamma_g/2)


    
    # set observers
    add_observer!(model, (1,), ("e<-e",), name="e")
    add_observer!(model, (1,), ("g<-g",), name="g")
    add_observer!(model, (1,), ("r<-r",), name="r")
    add_observer!(model, (1,), ("t<-t",), name="t")
    add_observer!(model, (2,), ("n",), name="n")
    
    s1 = zeros(ComplexF64, 4,4)
    s1[1,1] = 1
    s2 = bosonic_thermal_state(d, nbar)

    # model.state = tensor_product(s1, s2) 
    # results = evolve!(model, t=t, mt=0.5, driver="krylov")
    # obs = results.observables

    state = tensor_product(s1, s2) 
    ham, observer = model.h, model.observer

    L = size(state, 1)
    state = reshape(state, L*L)

    ham = matrix(ham)
    observer = immrep(observer)
    result = measure(observer, reshape(state, L, L))
    obs = Observables(result)

    mt = 0.5
    n = round(Int, t/mt)
    for tj in 1:n
        println("time step $tj....")
        @time state = evolve(ham, mt, state, ishermitian=false)
        tmp = measure(observer, reshape(state, L, L))
        append!(obs, tmp)
    end


    path_name = "result/ExactOmegag$(omega_L1)d$(d)nave$(nbar)t$(t).txt"

    println("save results to $path_name...")
    result = JSON.json(todict(obs))

    open(path_name, "w") do f
        write(f, result)
    end
    # return obs
end

nu = 2 * pi * 1.3
Delta_g = 2 * pi * 108 
Delta_r = 2 * pi * 108 
Delta_t = 2 * pi * 117.33
eta = 0.15

# Omega_L1 = 2 * pi * 1.
Omega_L2 = 2 * pi * 23.7 
costheta = sqrt(2)/2
gamma_g = 2 * pi * 6.89 
gamma_r = 2 * pi * 13.79  


d = 10
nave = 1.

omegas = [1.12, 1.41, 1.58, 1.93, 2.23, 2.5, 2.74, 5.0, 5.25, 5.48, 5.7, 5.92, 6.12]

# omegas = [1.12]

# omegas =  [2.236, 2.734, 3.16, 3.54, 3.87, 7.07, 7.42, 7.75, 8.06, 8.37, 8.66]

for i in omegas
    Omega_L1 = 2 * pi * i 
    @time four_level_eit_cooling(d, Delta_g, Delta_r, Delta_t, nu, eta, costheta, Omega_L1, Omega_L2, gamma_g, gamma_r, nave, 500.)
end



