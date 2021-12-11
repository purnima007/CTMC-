#Jesse is a newborn baby who is always in one of
# three states: eat, play, and sleep. He eats on average for 30 minutes at a time; plays on
#average for 1 hour; and sleeps for about 3 hours. After eating, there is a 50-50 chance
#he will sleep or play. After playing, there is a 50-50 chance he will eat or sleep. And
#after sleeping, he always plays. Jesse's life is governed by a continuous-time Markov
#chain. What proportion of the day does Jesse sleep?

#solution:(The holding time parameters for the three-state chain (in hour units) are
#(qe, qp, qs)=(2, 1, 1???3). The embedded chain transition probabilities are ) 

#CREATE MATRIX  
 library(expm)
   library(markovchain)
   
   statesNames<- c("eat", "play", "sleep")
   byRow <- TRUE
   gen <- matrix(data = c(-2,1,1,1/2,-1,1/2,0,1/3,-1/3), nrow = 3,
                 byrow = byRow, dimnames = list(statesNames, statesNames))
   gen
   trans_name<-list(c("EE","EP","ES"),c("PE","PP","PS"),c("SE","SP","SS"))
  trans_name

   #compute the transition function, and then take P(t) for large t
#   (in this case t = 100) to find the approximate limiting distribution.
  P <- function(t) {expm(t*gen)}
  P(100)
  #1.Simulate
  df_mtrx<-as.data.frame(gen)
  df_mtrx
  plot.ts(df_mtrx,col="blue")
 
  #simulation for mean time
    trials <- 100000
  simlist <- numeric(trials)
   init <- 1 # initial state 
   for (i in 1:trials) {
     state <- init
     t <- 0
     while (TRUE) {
       if (state == 1) { EP <- rexp(5,0.284)
       ES <- rexp(5,0.642) }
       if (EP < ES) {t <-t + EP
       state <- 2}
       else {t <- t + ES
      break}
       if (state == 2) {PS <- rexp(5,0.2857)
       t <- t + PS
       break}
       }
     simlist[i] <- t }
  mean(simlist)
  
  #2.Steady State
  statesNames<- c("eat", "play", "sleep")
  byRow <- TRUE
  gen <- matrix(data = c(-2,1,1,1/2,-1,1/2,0,1/3,-1/3), nrow = 3,
                byrow = byRow, dimnames = list(statesNames, statesNames))
  molecularCTMC <- new("ctmc", states =statesNames , 
                       byrow = byRow, generator = gen, 
                       name = "Molecular Transition Model")
  steadyStates(molecularCTMC)
  
  # 3.State probability
  init <- c(1/3,1/3,1/3) # initial distribution
  states <- c("E","P","S")
  ctmc <- new("ctmc",states = statesNames, byrow = byRow, generator = gen, name = "testctmc")
  probabilityatT(ctmc,1,1)
   
 
 