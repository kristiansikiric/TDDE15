library(HMM)
###### Assignment 1 ######
transprob = matrix(1:100,ncol = 10)*0
for(i in 1:10)
  for(j in 1:10)
    if(i == j | j == i +1)
    transprob[i,j] = 0.5
#Special case to make the path a ring
transprob[10,1] = 0.5

emissionprob = matrix(1:100,ncol = 10)*0
for(i in 1:10)
  for(j in 1:10)
    if(i == j | j == i +1 | j == i+2 | j == i-1 | j == i -2)
      emissionprob[i,j] = 0.2
#Some special cases
emissionprob[10,1] = 0.2
emissionprob[10,2] = 0.2
emissionprob[9,1] = 0.2
emissionprob[1,10] = 0.2
emissionprob[1,9] = 0.2
emissionprob[2,10] = 0.2

#Make sure the robot is in each sector at least two time steps
#states = c("S1","s1","S2","s2","S3","s3","S4","s4","S5","s5")
#symbols = c("S1","s1","S2","s2","S3","s3","S4","s4","S5","s5")
#startprob = c(0.2,0,0.2,0,0.2,0,0.2,0,0.2,0)
states = 1:10
symbols = 1:10


hmm = initHMM(States = states, 
              Symbols = symbols, 
              transProbs = transprob,
              startProbs = startprob,
              emissionProbs = emissionprob)

###### Assignment 2 ######
sim = simHMM(hmm,100)
sim$states
sim$observation

###### Assignment 3 ######
f = exp(forward(hmm,sim$observation))
b = exp(backward(hmm,sim$observation))

print(f)
filtering = f/sum(f)
smoothing = (f*b)/sum(f*b)

#Normalized probability distributions
prob.dist.filter = prop.table(filtering,2)
prob.dist.smoothing = prop.table(smoothing,2)

path = viterbi(hmm, sim$observation)

###### Assignment 4 ######
path.tab = table(Actual = sim$states,Predicted = path)
path.acc = sum(diag(path.tab))/sum(path.tab)

predict.filter = apply(prob.dist.filter,2,which.max)
filter.bool = states[predict.filter] == sim$states
filter.acc = sum(filter.bool == TRUE)/length(filter.bool)

predict.smoothing = apply(prob.dist.smoothing,2,which.max)
smoothing.bool = states[predict.smoothing] == sim$states
smoothing.acc = sum(smoothing.bool == TRUE)/length(smoothing.bool)

###### Assignment 5 ######
#Rerun 1-4 with different sample sizes...
###### Assignment 6 ######
library(entropy)
sim = simHMM(hmm,200)
f.200 = exp(forward(hmm,sim$observation))
filtering.200 = f.200/sum(f.200)
prob.dist.filter.200 = prop.table(filtering.200,2)


entropy.empirical(prob.dist.filter)
entropy.empirical(prob.dist.filter.200)
###### Assignment 7 ######
res = transprob %*% prob.dist.filter[,100]


