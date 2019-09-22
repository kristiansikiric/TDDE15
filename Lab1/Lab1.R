###### Init ######
library(bnlearn)
library(gRain)
data('asia')

###### Functions ######
predict.bn = function(junction.tree,data,target_nodes, evidence_nodes)
{
  predictions = rep(0,dim(data)[1])
  for (i in 1:dim(data)[1]) {
    states = NULL
    for (node in evidence_nodes) {
      states[node] = ifelse(data[i,node]=="yes","yes","no")
    }
    evidence = setEvidence(junction.tree,
                           nodes = evidence_nodes,
                           states = states)
    prob = querygrain(evidence,
                      nodes = target_nodes)
    predictions[i] = ifelse(prob$S["yes"] >= prob$S["no"],"yes","no")
  }
  return(predictions)
}

###### Exercise 1 ######
set.seed(1234)
bn = hc(x = asia)
bn_with_restart = hc(x = asia, restart=10, score = 'aic')
plot(bn)
plot(bn_with_restart)
all.equal(bn,bn_with_restart)

# Answer: This is due to ...


###### Exercise 2 ######
smp_size = dim(asia)[1]
set.seed(123)
id = sample(1:smp_size,floor(smp_size*0.8))
train = asia[id,]
test = asia[-id,]

bn = hc(x=train)
fitted.bn = bn.fit(bn,data=train)
junction.tree = compile(as.grain(fitted.bn)) #RIP

true.bn = model2network("[A][S][T|A][L|S][B|S][D|B:E][E|T:L][X|E]")
true.fitted.bn = bn.fit(true.bn,data=train)
true.junction.tree=compile(as.grain(true.fitted.bn))

nodes = colnames(asia)
evidence_nodes = nodes[nodes != "S"]
target_nodes = nodes[nodes == "S"]

prediction = predict.bn(junction.tree = junction.tree,
                        data = test,
                        target_nodes = target_nodes,
                        evidence_nodes = evidence_nodes)
true.prediction = predict.bn(junction.tree =  true.junction.tree,
                             data = test,
                             target_nodes = target_nodes,
                             evidence_nodes = evidence_nodes)
print(table(prediction,test$S))
print(table(true.prediction,test$S))

###### Exercise 3 ######
markov.blanket = mb(fitted.bn,"S")
true.markov.blanket = mb(true.fitted.bn,"S")
prediction.mb = predict.bn(junction.tree = junction.tree,
                     data = test,
                     target_nodes = target_nodes,
                     evidence_nodes = markov.blanket)
true.prediction.mb = predict.bn(junction.tree = true.junction.tree,
                     data = test,
                     target_nodes = target_nodes,
                     evidence_nodes = true.markov.blanket)
print(table(prediction.mb,test$S))
print(table(true.prediction.mb,test$S))

###### Exercise 4 ######
naive.bayes.dag = empty.graph(nodes = nodes)
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "A")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "D")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "X")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "E")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "B")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "L")
naive.bayes.dag = set.arc(naive.bayes.dag, from = "S", to = "T")
plot(naive.bayes.dag)
nbd.fit = bn.fit(naive.bayes.dag,data=train)
junc.tree = compile(as.grain(nbd.fit))

predict.nbd = predict.bn(junc.tree, test, target_nodes, evidence_nodes)
print(table(predict.nbd,test$S))
print(table(true.prediction,test$S))
