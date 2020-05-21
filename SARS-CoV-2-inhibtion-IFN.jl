using DifferentialEquations
using Plots
gr(size=(500,500), show = true)  #use the gr backend for plotting
#Function for the lorenz equation
#u[1] = x; u[2] = y; u[3] = z
function P1(du,u,p,t)
   lambda = 0.003241
   d = 10e-3
   k = 0.095
   dau=0.625
   p=0.375
   c=3.0
  B=0.0
  tau=0.0
  ek=0.6
  ep=0.6
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end




function P6(du,u,p,t)
  lambda = 0.003241
d = 10e-3
k = 0.095
dau=0.625
p=0.375
c=3.0
 B=0.8
 tau=0.1
  ek=0.6
  ep=0.6
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end
u0 = [3.0;2.0;0.061]                          #intial conditions
tspan = (0.0,14.0)                     #start and end time

 prob1 = ODEProblem(P1,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol1= solve(prob1)

prob6 = ODEProblem(P6,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol6= solve(prob6)

plot(sol6,vars=(0,3),xaxis=("time (days)"),yaxis=("log10 Viral Loads"),label=["SARS-CoV-2 60% inhbition + IFN "])
plot!(sol1,vars=(0,3),label=["60% inhbition no IFN"])
