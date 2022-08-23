#Use relevant packages
using JuMP, Ipopt, Plots, PGFPlotsX, LaTeXStrings

#Include necessary files and functions
include("Parameters.jl")
include("Model.jl")
include("Parameters.jl")
include("Simulation.jl")

#Define simulation settings
np   = 40;                      #Points to consider in the simulation
Pnet = LinRange(10000,2e6,np);  #Vector of power supply, in W

#Run simulation
sol, SSEl, ter = Simulation(np,par,Pnet)

#Plot results
plot(Pnet, sol["P"]/par[:maxP], ylims=(0,1.1), label = L"P", shape = :circle, xlabel = L"P_{net} [W]", ylabel = L"g/g_{max}", framestyle =:box, widen = false)
plot!(Pnet,sol["n_a_H2_o"]./sol["n_sepa_O2_o"]*100/2,ylims=(0,1.1), label = L"HTO", shape = :circle)
plot!(Pnet, sol["m_cw"]/26.6, label = L"\dot{m}_{cw}", shape = :circle)
plot!(Pnet,sol["TEl"]/(273+80), label = L"T^{El}", shape = :circle)
plot!(Pnet, sol["mel_lye_i"]/(par[:maxlye]), label = L"\dot{m}_{lye}", shape = :circle)