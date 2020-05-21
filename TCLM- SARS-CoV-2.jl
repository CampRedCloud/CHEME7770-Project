using DifferentialEquations
using Plots
gr(size=(500,500), show = true)  #use the gr backend for plotting

#Function for the lorenz equation
#u[1] = x; u[2] = y; u[3] = z
function lorenz!(du,u,p,t)
    lambda = 0.003241
   d = 10e-3
  k = 0.095
  dau=0.625
    p=0.375
  c=3.0
B=0.0
tau-0.0
epsk=0.0
epsp=0.0


 du[1] = lambda-d*u[1] -(1-epsk)*k*u[3]*u[1];               #dx/dt
 du[2] = (1-epsk)*k*u[3]*u[1]-dau*u[2];      #dy/dt
 du[3] = (1-epsp)*p*u[2]-c*u[3]-(t-tau)*B*u[3];     #dz/dt
end

u0 = [3.0;2.0;0.061]                      #intial conditions

tspan = (0.0,14.0)                     #start and end time
prob = ODEProblem(lorenz!,u0,tspan)     #Create an ODE problem for the Lorenz fxn


sol = solve(prob)                      #Solve the system


#Plot the results; the vars=(0,1) argument specifies to plot X (column 1 of sol)
#vs t (column 0 of sol)
plt1 = plot(sol,vars=(0,1), xaxis="time", yaxis = "log10 Uninfected cells ",label="SARS-CoV-2" )

display(plt1)


#Plot the results; the vars=(0,1) argument specifies to plot Y (column 1 of sol)
#vs t (column 0 of sol)
plt2 = plot(sol,vars=(0,2), xaxis="time", yaxis = "log10 Infected cells",label="SARS-CoV-2" )

display(plt2)


#Plot the results; the vars=(0,3) argument specifies to plot Z (column 3 of sol)
#vs t
plt3 = plot(sol,vars=(0,3), xaxis="time", yaxis = "log10 Viral loads",label="SARS-CoV-2")

display(plt3)
