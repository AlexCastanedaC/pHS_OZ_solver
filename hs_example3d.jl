# Our most relevant tools come from the OrnsteinZernike package
using OrnsteinZernike

#= One must first define the parameters of the system, such as the dimensionality and the thermal energy. For this example, we use the hard-sphere potential =#
dims = 3; kBT = 1.0;

potential = HardSpheres(1.0)

using DataFrames, CSV

#= For the purpose of our research paper, we compare the OZ results with
simulation data (Brownian Dynamics). Therefore, we extract the packing 
fractions (in the liquid phase) from a csv file =#

df = DataFrame(CSV.File("eos_bd_3d.csv"; limit = 14))
ϕ = df[:, :phi]

#= the solver requires the numerical density (so we obtain the reduced density
from the packing fraction array=# 
ρ = (ϕ .* 6) ./ π

# So far, only the SimpleLiquid system has been implemented in the OZ package

system = SimpleLiquid(dims, ρ[end], kBT, potential)

#= Our favorite method for solving the OZ equation is the Ng Iteration scheme
which has been implemented already in the package =#

method = NgIteration()

#= However, for solving at high densities, it is convenient to use the DensityRamp solver, which takes the main method as one of its arguments. The idea is to iteratively solve systems of increasing densities. Thus, we can easily obtain 
the equation of state =#

method2 = DensityRamp(method, ρ[1:end-1])

#= Moreover, we need a closure to solve the OZ equation; as a solid reference, we use the Percus-Yevick closure =#

closure = PercusYevick()

sol = @time solve(system, closure, method2)

z_oz = Float64[]

for i in 1:length(ρ)
	pressure = compute_virial_pressure(sol[i], SimpleLiquid(dims, ρ[i], kBT, potential))
	zfactor = pressure / (ρ[i] * kBT)
	push!(z_oz, zfactor)
end

df_2 = DataFrame(ϕ = ϕ, Z_sim = df[:, :Z], Z_OZ = z_oz)

CSV.write("Comparison_Z.csv", df_2)
