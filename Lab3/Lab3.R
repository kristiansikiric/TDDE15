####### Innit ######
set.seed(123)
z1 = runif(1,0,100)
T = 100
M = 100
sd_emission = 1
####### Functions ######
transition = function(z,sd = 1){
  comp = sample(0:2,1)
  return (rnorm(1,z+comp,sd))
}

emission = function(z,sd = sd_emission){
  comp = sample(-1:1,1)
  return (rnorm(1,z+comp,sd))
}

density.emission = function(z,x,sd = sd_emission){
  p = dnorm(x,z,sd) + dnorm(x,z-1,sd) + dnorm(x,z+1,sd)
  return(mean(p))
}

particle_filter = function(weight = TRUE, obs){
  w = rep(1,100)
  particles = runif(M,0,100)
  expected = rep(0,100)
  X = matrix(NA,nrow = T, ncol = M)
  for (t in 1:T){
    if(weight){
      for (m in 1:M){
        w[m] = density.emission(z = particles[m],x = obs[t])
      }
    }
    particles = sample(particles, size = 100, replace = TRUE, prob = w)
    X[t,] = particles
    
    w = w / sum(w)
    tmp = NULL
    for(m in 1:M){
      tmp[m] = w[m]*particles[m]
    }
    
    expected[t]= sum(tmp)
    for(m in 1:M){
      particles[m] = transition(particles[m])
    }
  }
  return (list("Expected" = expected,"Particles" = X))
}

SSM = function(){
  states = rep(0,T)
  observations = rep(0,T)
  states[1] = z1
  
  observations[1] = emission(states[1])
  for (t in 2:T){
    states[t] = transition(states[t-1])
    observations[t] = emission(states[t])
  }
  return (list("States" = states, "Observations" = observations))
}

####### Assignment 1 ######
#SSM
ssm = SSM()
states = ssm$States
observations = ssm$Observations

kalman = particle_filter(obs = observations)
expected = kalman$Expected
X = kalman$Particles

plot(1:T,expected,type = "l", col = "red")
lines(1:T,states, col = "black")
points(rep(1,100), X[1,],pch = 1)
points(rep(30,100), X[30,],pch = 1)
points(rep(50,100), X[50,],pch = 1)
points(rep(100,100), X[100,],pch = 1)

mean((states-expected)^2)


###### Assignment 2 ######
sd_emission = 5
ssm.5 = SSM()
states.5 = ssm.5$States
observations.5 = ssm.5$Observations

kalman = particle_filter(obs = observations.5)
expected.5 = kalman$Expected
X.5 = kalman$Particles

plot(1:T,states.5, col = "black", type = "l")
lines(1:T,expected.5,col = "green")
points(rep(1,100), X.5[1,],pch = 1)
points(rep(30,100), X.5[30,],pch = 1)
points(rep(50,100), X.5[50,],pch = 1)
points(rep(100,100), X.5[100,],pch = 1)
mean((states.5-expected.5)^2)

sd_emission = 50
ssm.50 = SSM()
states.50 = ssm.50$States
observations.50 = ssm.50$Observations

kalman = particle_filter(obs = observations.50)
expected.50 = kalman$Expected
X.50 = kalman$Particles

plot(1:T,expected.50, col = "blue", type = "l")
lines(1:T,states.50,col = "black")
points(rep(1,100), X.50[1,],pch = 1)
points(rep(30,100), X.50[30,],pch = 1)
points(rep(50,100), X.50[50,],pch = 1)
points(rep(100,100), X.50[100,],pch = 1)
mean((states.50-expected.50)^2)

###### Assignment 3 ######
set.seed(123)
sd_emission = 1
ssm.w = SSM()
states.w = ssm.w$States
observations.w = ssm.w$Observations

kalman = particle_filter(weight = FALSE,obs = observations.w)
expected.w = kalman$Expected
X.w = kalman$Particles

plot(1:T,states.w, col = "black", type = "l", ylim = c(0,max(X.w)))
lines(1:T,expected.w,col = "yellow")
points(rep(1,100), X.w[1,],pch = 1)
points(rep(30,100), X.w[30,],pch = 1)
points(rep(50,100), X.w[50,],pch = 1)
points(rep(100,100), X.w[100,],pch = 1)
mean((states.w-expected.w)^2)

