library(sna)
library(tsna)
library(ndtv)
library(networkDynamic)
library(network)
PHStaticEdges <- read.csv("TNAWR_StaticEdgelist.csv")

PHVertexAttributes <- read.csv(
  "TNAWR_VertexAttributes.csv",
  stringsAsFactors = FALSE
)

thenetwork <- network(
  PHStaticEdges,
  vertex.attr = PHVertexAttributes,
  vertex.attrnames = c("vertex.id", "name", "region"),
  directed = FALSE,
  bipartite = FALSE
)
#plot(thenetwork)


# Import Temporal Network Data
PHDynamicNodes <- read.csv("TNAWR_DynamicNodes.csv")
PHDynamicEdges <- read.csv("TNAWR_DynamicEdges.csv")


dynamicCollabs <- networkDynamic(
  thenetwork,
  edge.spells = PHDynamicEdges,
  vertex.spells = PHDynamicNodes
)

# Calculate how to plot an animated version of the dynamic network
compute.animation(
  dynamicCollabs,
  animation.mode = "kamadakawai",
  slice.par = list(
    start = 1260,
    end = 1300,
    interval = 1,
    aggregate.dur = 20,
    rule = "any"
  )
)


# Render the animation and open it in a web brower
render.d3movie(
  dynamicCollabs,
  displaylabels = TRUE,
  # This slice function makes the labels work
  vertex.tooltip = function(slice) {
    paste(
      "<b>Name:</b>", (slice %v% "name"),
      "<br>",
      "<b>Region:</b>", (slice %v% "region")
    )
  }
)


#########################xx


nodes <- read.csv("Dataset1-Media-Example-NODES.csv", header=T, as.is=T)
links <- read.csv("Dataset1-Media-Example-EDGES.csv", header=T, as.is=T)
net3 <- network(links,  vertex.attr=nodes, matrix.type="edgelist", 
                loops=F, multiple=F, ignore.eval = F)

#net3[,]
net3 %n% "net.name" <- "Media Network" #  network attribute
net3 %v% "media"    # Node attribute
net3 %e% "type"     # Edge attribute


net3 %v% "col" <- c("gray70", "tomato", "gold")[net3 %v% "media.type"]
plot(net3, vertex.cex=(net3 %v% "audience.size")/7, vertex.col="col")

par(mar=c(0,0,0,0))

render.d3movie(net3, usearrows = T, displaylabels = T, bg="#111111", 
               vertex.border="#ffffff", vertex.col =  net3 %v% "col",
               vertex.cex = (net3 %v% "audience.size")/8, 
               edge.lwd = (net3 %e% "weight")/3, edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (net3 %v% 'media') , "<br>",
                                      "<b>Type:</b>", (net3 %v% 'type.label')),
               edge.tooltip = paste("<b>Edge type:</b>", (net3 %e% 'type'), "<br>", 
                                    "<b>Edge weight:</b>", (net3 %e% 'weight') ),
               launchBrowser=F, filename="Media-Network.html" )  


# dynamic

vs <- data.frame(onset=0, terminus=50, vertex.id=1:17)
es <- data.frame(onset=1:49, terminus=50, 
                 head=as.matrix(net3, matrix.type="edgelist")[,1],
                 tail=as.matrix(net3, matrix.type="edgelist")[,2])

net3.dyn <- networkDynamic(base.net=net3, edge.spells=es, vertex.spells=vs)


compute.animation(net3.dyn, animation.mode = "kamadakawai",
                  slice.par=list(start=0, end=50, interval=1, 
                                 aggregate.dur=1, rule='any'))

render.d3movie(net3.dyn, usearrows = T, 
               displaylabels = T, label=net3 %v% "media",
               bg="#ffffff", vertex.border="#333333",
               vertex.cex = degree(net3)/2,  
               vertex.col = net3.dyn %v% "col",
               edge.lwd = (net3.dyn %e% "weight")/3, 
               edge.col = '#55555599',
               vertex.tooltip = paste("<b>Name:</b>", (net3.dyn %v% "media") , "<br>",
                                      "<b>Type:</b>", (net3.dyn %v% "type.label")),
               edge.tooltip = paste("<b>Edge type:</b>", (net3.dyn %e% "type"), "<br>", 
                                    "<b>Edge weight:</b>", (net3.dyn %e% "weight" ) ),
               launchBrowser=T, filename="Media-Network-Dynamic.html",
               render.par=list(tween.frames = 30, show.time = F),
               plot.par=list(mar=c(0,0,0,0)), output.mode = "HTML" )

plot(tEdgeFormation(dynamicCollabs, time.interval = .25), 
     xlab = "Time (years)",
     ylab = "Edge formation",
     main = "Edge formation rate")


dynamicBetweenness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1260,
  end = 1320,
  time.interval = 1,
  aggregate.dur = 20,
  FUN = "betweenness"
)
plot(dynamicBetweenness,xlab = "Time (years)",
     ylab = "Betweenness",
     main = "Betweenness centralization")


dynamicCloseness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1260,
  end = 1320,
  time.interval = 1,
  aggregate.dur = 20,
  FUN = "closeness"
)
plot(dynamicCloseness,xlab = "Time (years)",
     ylab = "Closeness",
     main = "Closeness centrality")


dynamicDegree <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1260,
  end = 1320,
  time.interval = 1,
  aggregate.dur = 20,
  FUN = "closeness"
)
windows()
plot(dynamicDegree,xlab = "Time (years)",
     ylab = "Degree",
     main = "Degree centrality")

windows()
filmstrip(dynamicCollabs, displaylabels = FALSE)
