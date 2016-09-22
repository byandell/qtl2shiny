library(networkD3)

# Create fake data

src <- c("A", "A", "A", "A",
         "B", "B", "C", "C", "D")
target <- c("B", "C", "D", "J",
            "E", "F", "G", "H", "I")
networkData <- data.frame(src, target)

networkData <- read.csv("inst/doc/ShinyLogicFlow.csv")
networkData <- read.csv("inst/doc/ShinyLogicServer.csv")

# Plot
simpleNetwork(networkData, textColour = "black", fontSize=9)

library(igraph)
ig <- graph_from_edgelist(as.matrix(networkData[,1:2]))
