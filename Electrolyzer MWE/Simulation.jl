function Simulation(np,par,Pnet)
    #Prepare the simulation
    SSEl =  Build_Model(par)
    nvars = length(all_variables(SSEl))
    
    #Pre-allocate arrays and vectors
    vars = Vector{VariableRef}(undef,nvars);
    vals = Array{Float64}(undef,nvars,np);
    ter = Vector{Any}(undef,np);
    
    for i in 1:np
        global vars, sol
        par[:Pnet] = Pnet[i]
        ## Setting up JuMP model 
        SSEl = Build_Model(par)
        optimize!(SSEl)
        #Save values of variables
        ter[i] = termination_status(SSEl)
        vars = name.(all_variables(SSEl))
        vals[:,i] = value.(all_variables(SSEl))
        if i == 1
            sol = Dict(vars[1] => vals[1,:])
        end 
        for j in 2:nvars
            push!(sol, vars[j] => vals[j,:])
        end 
    end 
    return sol, SSEl, ter
end 