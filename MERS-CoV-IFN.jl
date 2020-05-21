using DifferentialEquations
using Plots
gr(size=(500,500), show = true)  #use the gr backend for plotting
#Function for the lorenz equation
#u[1] = x; u[2] = y; u[3] = z
function P1(du,u,p,t)
   lambda = 0.003241
   d = 10e-3
   k = 0.86
   dau=0.4
   p=0.513
   c=3.0
  B=0.8
  tau=0.1
  ek=0.0
  ep=0.0
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end

function P2(du,u,p,t)
   lambda = 0.003241
   d = 10e-3
   k = 0.86
   dau=0.4
   p=0.513
   c=3.0
  B=0.8
  tau=1.0
  ek=0.0
  ep=0.0
  du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
  du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
  du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end

function P3(du,u,p,t)
   lambda = 0.003241
   d = 10e-3
   k = 0.86
   dau=0.4
   p=0.513
   c=3.0
 B=0.8
 tau=2.0
  ek=0.0
  ep=0.0
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end



function P4(du,u,p,t)
   lambda = 0.003241
   d = 10e-3
   k = 0.86
   dau=0.4
   p=0.513
   c=3.0
 B=0.8
 tau=3.0
  ek=0.0
  ep=0.0
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt/dt
end



function P6(du,u,p,t)
  lambda = 0.003241
  d = 10e-3
k = 0.86
dau=0.4
p=0.513
c=3.0
 B=0.0
 tau=0.0
  ek=0.0
  ep=0.0
 du[1] = lambda-d*u[1] -(1-ek)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-ek)*k*u[3]*u[1]-(dau)*u[2];      #dy/dt
 du[3] = (1-ep)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end
u0 = [3.0;2.0;0.061]                          #intial conditions
tspan = (0.0,14.0)                     #start and end time

 prob1 = ODEProblem(P1,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol1= solve(prob1)
prob2 = ODEProblem(P2,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol2= solve(prob2)
 prob3 = ODEProblem(P3,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol3= solve(prob3)
 prob4 = ODEProblem(P4,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol4= solve(prob4)
 prob5 = ODEProblem(P5,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol5= solve(prob5)
prob6 = ODEProblem(P6,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol6= solve(prob6)

plot(sol6,vars=(0,3),xaxis=("time (days)"),yaxis=("log10 Viral Loads"),label=["MERS-CoV no IFN  "])
plot!(sol1,vars=(0,3),label=["tau=0.1"])
plot!(sol2,vars=(0,3),label=["tau=1.0"])
plot!(sol3,vars=(0,3),label=["tau=2.0"])
plot!(sol4,vars=(0,3),label=["tau=3.0"])
#plot!(sol5,vars=(0,3),label=["tau=3.0"])
