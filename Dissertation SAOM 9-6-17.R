# Stochastic Actor Oriented Modeling 7/6/17
##############################################
# Create Network Slices and put them in a list
# Differnt list for each year -- each year has 12 time points
# Will model the network each year and see if from year to year they are different
# QUESTION -- INCLUDE PRIOR COVARIATES IN YEAR?
# QUESTION -- CAN I MODEL ALL 108 TIME PERIODS AS 1 NETOWRK?
# QUESTION -- SHOULD I CONVERT TO NETWORK FOR EACH YEAR?
# RSiena requires that each timepoint be a matrix and that all timepoints are an array of these matrices.
# Each matrix in the array must have the same dimensions.

setwd("e:/profiles/desktop/Dissertation")

# Open data
all.data = read.csv("interaction data for RSiena 8-16-17.csv", na.strings = NA)
#n=40954
#vars=4
names(all.data)
names(all.data)[1] <- "org"
head(all.data, n=25)
table(all.data$t.month)
hist(all.data$t.month)

# create year variable
all.data$year <- 0
all.data$year[all.data$t.month > 1 & all.data$t.month < 14] <- 1
all.data$year[all.data$t.month > 13 & all.data$t.month < 26] <- 2
all.data$year[all.data$t.month > 25 & all.data$t.month < 38] <- 3
all.data$year[all.data$t.month > 37 & all.data$t.month < 50] <- 4
all.data$year[all.data$t.month > 49 & all.data$t.month < 62] <- 5
all.data$year[all.data$t.month > 61 & all.data$t.month < 74] <- 6
table(all.data$t.month, all.data$year)

###########################################################################################################
# Create the network matrices -- network and dyadic 
# meeting is all is binary (0/1) based on committee attendance
# TA is binary (0/1) based on TA attendance
# contract is binary based on a contract for 12 months
library(Matrix)
library(igraph)

# create list by year
all <- split(all.data, all.data$year)

# create adjacency matrices from each dataframe (to find the zero links in the time period)
# First create a sparse matrix of orgs by meeting/contract/TA
# Then create an adjacency matrix of orgs by orgs IS THIS WEIGHTED? IF SO SAVE FOR DYADIC VARIABLE
# Then create a igraph object and simplify it to remove loops and duplicates
# Then convert back to adjacency matrix
# Function Adapted from Solomon Messing
# https://solomonmessing.wordpress.com/2012/09/30/working-with-bipartiteaffiliation-network-data-in-r/
matrices <- function(df) {
  test = df
  test.m <- spMatrix(nrow=length(unique(test$org)),
                     ncol=length(unique(test$interactionID)),
                     i = as.numeric(factor(test$org)),
                     j = as.numeric(factor(test$interactionID)),
                     x = rep(1, length(as.numeric(test$org))) )
  row.names(test.m) <- levels(factor(test$org))
  colnames(test.m) <- levels(factor(test$interactionID))
  test.cp <- tcrossprod(test.m)
  test.g <- graph_from_adjacency_matrix(test.cp, weighted = NULL, mode = 'undirected', diag = FALSE)
  test.g <- simplify(test.g, remove.multiple = TRUE, remove.loops = TRUE)
  test.g.a <- as_adjacency_matrix(test.g, type = "both", names = TRUE, sparse = FALSE)
  df = test.g.a
}

matrices.all <- lapply(all, matrices)


# Take a look
head(matrices.all$`2`)
sum(matrices.all$`2`)

# remove what I don't need
rm(all)

###################################################################################################
# Problem is that graphs have different dimensions. Need to be the same -- and to include all active
# orgs in the year

names.all <- lapply(matrices.all, rownames)
names.all.all <- unlist(names.all)
names.all.unique <- unique(names.all.all)
names.all.list <- list(names.all.unique, names.all.unique)

# remove what I don't need
rm(names.all)
rm(names.all.all)
#rm(names.all1.list)

#rm(names.all6.list)

# create empty matrices
empty.all <- matrix(data = NA, nrow = 402, ncol = 402, dimnames =  names.all.list)

# convert empty matrix to dataframe
library(reshape2)
df.empty.all <- setNames(melt(empty.all), c("Var1", "Var2", "value"))
df.empty.all$value <- as.numeric(df.empty.all$value)
df.empty.all$Var1 <- as.character(df.empty.all$Var1)
df.empty.all$Var2 <- as.character(df.empty.all$Var2)

head(df.empty.all)

# convert time slice matrices to dataframes
df.all.slices <- lapply(matrices.all, melt)

dim(df.all.slices$`2`)
head(df.all.slices$`2`)
tail(df.all.slices$`2`)

# join full and empty dataframes -- actually edgelists
# here I make no link with an org that did not attend any all, or have a contract, or attend TA
# this month a zero -- meaning the org could join a meeting in this network if wanted to.

library(dplyr)
join <- function(df) {
  temp1 <- df
  temp1$Var1 <- as.character(temp1$Var1)
  temp1$Var2 <- as.character(temp1$Var2)
  temp <- left_join(df.empty, temp1, by=c("Var1", "Var2"))
  temp$value.y[is.na(temp$value.y)] <- 0 #making this 0 means no link, making it 10 means structural zero
  temp$value.x <- NULL
  temp$value <- temp$value.y
  temp$value.y <-NULL
  df <- temp
}

df.empty <- df.empty.all
join.all <- lapply(df.all.slices, join)

class(join.all)
class(join.all$`2`)
dim(join.all$`2`)

# create matrices
matrices2 <- function(df) {
  test <- empty
  for (i in 1:NROW(df)) test[ df[i,1], df[i,2] ] <- df[i,3]
  df <- test
}

empty <- empty.all
mat.all <- lapply(join.all, matrices2)

dim(mat.all$`2`)

# Convert them all to dgTMatrix which is the object type required by RSiena
m.all <- lapply(mat.all, as, "dgTMatrix")

class(m.all)
class(m.all$`5`)
m.all$`5`

# lost the row and column names -- why?
m.all <- lapply(m.all, "rownames<-", names.all.unique); m.all <- lapply(m.all, "colnames<-", names.all.unique)
class(m.all)
class(m.all$`1`)
m.all$`6`

# remove what I don't need now
rm(empty.all)
rm(df.all.slices)
rm(join.all)
rm(mat.all)
rm(matrices.all)
rm(names.all.list)
rm(empty)

#########################################################################################
# Create the portfolio and resolution DVs and covariate files
library(tidyr)

org.data <- read.csv("Data 9-4-17.csv")
#n=2412
#vars=27
names(org.data)
names(org.data)[1] <- "org"
names(org.data)
head(org.data)

# resolution (binary 0/1 -- once turned on does not turn off)
resolution <- as.data.frame(org.data[, c(1:2,23)])
resolution <- spread(resolution, Year, resolution_cum)
resolution <- as.matrix(resolution[, 2:7])

# portfolio
portfolio <- as.data.frame(org.data[, c(1:2,25)])
portfolio <- spread(portfolio, Year, portfolio_cum)
portfolio <- as.matrix(portfolio[, 2:7])

# TA
TA <- as.data.frame(org.data[, c(1:2,27)])
TA <- spread(TA, Year, TA_cum)
TA <- as.matrix(TA[, 2:7])

# contract
contract <- as.data.frame(org.data[, c(1:2,13)])
contract <- spread(contract, Year, contract_amount)
contract <- as.matrix(contract[, 2:7])


# provider
provider <- as.data.frame(org.data[, c(1:2,4)])
provider <- spread(provider, Year, sector)
provider <- provider[, 2]
provider <- as.numeric(provider)

#########################################################################################################
# Move on now to modeling
# Create 3 dependent networks: all, contracts, TA (3 different forms of interactions)
##########################################################################################################
library(RSiena)

# Create the network objects -- the basic network structure
interactions <- sienaDependent(m.all)
interactions

# Create the dependent variable -- cumulative portfolios
# allowOnly constrains the movement to up (or down) only -- TRUE is the default but I include it here to
# remind myself that I am accounting for 
port <- sienaDependent(portfolio, type = "behavior", allowOnly = TRUE)

# Create non-varying monadic variable
prov <- coCovar(val = provider)

# Create varying monadic variables
cont <- varCovar(contract)
res <- varCovar(resolution)
TA <- varCovar(TA)

# Create systems
# Leave out TA until the 4th year becuase it is always 0 and model with not run "computationally singular"
system <- sienaDataCreate(interactions, res, port, cont, TA, prov)
system

# Print a desctiptives report-- this goes to the working directory as a text file
print01Report(system, modelname="System")

# Now, create the basic effects objects
effects <- getEffects(system)

# All the effects that are available given the structure of this data set can be seen from
# 890 possible effects
effectsDocumentation(effects)

##########################################################################################################
# 4 models per Steglich, Snijders, & Pearson 2012 Dynamic Networsk and Behavior: Separating Selection from Influence
# They call these 4 degress of complexity for the coevolution (Trend, Control, Full, Straw)
# Now create the algorithm
# Can only specify network modelType = OR behModelType = (NOT BOTH!)
# REMEMBER -- when I find a solution I like, then I need to run it with lessMem = FALSE to be able to do 
# the model fit statistics
# test it with small n3 -- then increase to 3000 for publication
# cond can only be true if there is only 1 dependent variable. When there is a behavioral DV, 
# then it must be FALSE becuase the network itself is the 1st DV
# lessMem TRUE will be faster. But when have model I like need to be FALSE to run sienaGOF on the model

# TREND MODEL (NULL MODEL)
a.trend <- sienaAlgorithmCreate(projname = "dissertation", cond = FALSE, n3 = 100, dolby = FALSE,
                          condname = "port1", behModelType = 1, lessMem = TRUE)
model.trend <- siena07(a.trend, data = system, effects = effects, useCluster = TRUE, nbrNodes = 30)


# Now add effects iteratively and based on theory
# these are my understandings of these concepts now 8/17/17
# First add effects that describe the network -- it is useful to simply describe how the network works
# outdegree is the overall tendency to have ties, it is the density of the network
effects1 <- includeEffects(effects1, density, type = "eval", name = "interactions1")
# transitive ties -- the network tends to display triadic closure -- both immediately and at distance-2
effects1 <- includeEffects(effects1, transTriads, type = "eval", name = "interactions1")
effects1 <- includeEffects(effects1, nbrDist2, type = "eval", name = "interactions1")
# balance -- tendency to have ties to others with similar structure (in terms of degrees)
effects1 <- includeEffects(effects1, balance, type = "eval", name = "interactions1")
# betweeness: tendencey to occupy an intermediary position between unrelated others (broker position)
effects1 <- includeEffects(effects1, between, type = "eval", name = "interactions1")
# homofilia (similiarity) there is similarity in ties by provider_bin
effects1 <- includeEffects(effects1, simX, type = "eval", name = "interactions1", interaction1 = "prov1")

# Then add effects on the DV of interest -- the resolution and the portfolio
# org having a portfolio is assimilated to connected orgs (an org influences other orgs
effects1 <- includeEffects(effects1, avSim, type = "eval", name = "port1", interaction1 = "interactions1")
# effects of being a provider, having a resolution, amount of contracts (and amount of TA) on having a portfolio
effects1 <- includeEffects(effects1, effFrom, type = "eval", name = "port1", interaction1 = "prov1")
effects1 <- includeEffects(effects1, effFrom, type = "eval", name = "port1", interaction1 = "res1")
effects1 <- includeEffects(effects1, effFrom, type = "eval", name = "port1", interaction1 = "cont1")
#effects1 <- includeEffects(effects1, effFrom, type = "eval", name = "port1", interaction1 = "ta1")
# Add quadtratic shape -- needed for behavioral DVs
effects1 <- includeEffects(effects1, quad, type = "eval", name = "port1")
effects1

#####################################################################################################
# MODEL 6
# Now add effects iteratively and based on theory
# these are my understandings of these concepts now 8/17/17
# First add effects that describe the network -- it is useful to simply describe how the network works
# outdegree is the overall tendency to have ties, it is the density of the network
effects6 <- includeEffects(effects6, density, type = "eval", name = "interactions6")
# transitive ties -- the network tends to display triadic closure -- both immediately and at distance-2
effects6 <- includeEffects(effects6, transTriads, type = "eval", name = "interactions6")
effects6 <- includeEffects(effects6, nbrDist2, type = "eval", name = "interactions6")
# balance -- tendency to have ties to others with similar structure (in terms of degrees)
effects6 <- includeEffects(effects6, balance, type = "eval", name = "interactions6")
# betweeness: tendencey to occupy an intermediary position between unrelated others (broker position)
effects6 <- includeEffects(effects6, between, type = "eval", name = "interactions6")
# homofilia (similiarity) there is similarity in ties by provider_bin
effects6 <- includeEffects(effects6, simX, type = "eval", name = "interactions6", interaction1 = "prov6")
# Then add effects on the DV of interest -- the resolution and the portfolio
# org having a portfolio is assimilated to connected orgs (an org influences other orgs
effects6 <- includeEffects(effects6, avSim, type = "eval", name = "port6", interaction1 = "interactions6")
# effects of being a provider, having a resolution, amount of contracts (and amount of TA) on having a portfolio
effects6 <- includeEffects(effects6, effFrom, type = "eval", name = "port6", interaction1 = "prov6")
effects6 <- includeEffects(effects6, effFrom, type = "eval", name = "port6", interaction1 = "res6")
effects6 <- includeEffects(effects6, effFrom, type = "eval", name = "port6", interaction1 = "cont6")
effects6 <- includeEffects(effects6, effFrom, type = "eval", name = "port6", interaction1 = "ta6")
# Add quadtratic shape -- it represents the basic shape of the DV (similar to the intercept)
effects6 <- includeEffects(effects6, quad, type = "eval", name = "port6")
effects6




