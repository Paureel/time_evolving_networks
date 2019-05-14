library(sna)
library(tsna)
library(ndtv)
library(networkDynamic)
library(network)
library(plotly)
library(ggplot2)
library(tidyr)
library(gplots)
library(leaflet)
library(geojsonR)
library(dplyr)
require(graphics)
library(TTR)
library(anytime)
library(forecast)
library(imputeTS)
library(reshape2)
library(matrixStats)
library(leaflet.extras)

df1 <- read.table("output-file-new.gz", header = F, 
                  stringsAsFactors = F, sep = "\t")
nodes_selected <- read.table("rnodes.txt", header = F, 
                             stringsAsFactors = F, sep = "\t")
file_js = FROM_GeoJson(url_file_string = "milano-grid.geojson")


attr <- file_js$features

idvec <- vector()
geomvec1 <- vector()
geomvec2 <- vector()
for (index in c(1:10000)){
  
  
  idvec <- c(idvec, attr[[index]][["id"]])
}

for (index in c(1:10000)){
  
  
  geomvec1 <- c(geomvec1, attr[[index]][["geometry"]][["coordinates"]][1,1])
  geomvec2 <- c(geomvec2, attr[[index]][["geometry"]][["coordinates"]][1,2])
}

names(df1) <- c("onset","tail", "head", "weight")

df1 <- df1[(df1$tail %in% nodes_selected$V1) & (df1$head %in% nodes_selected$V1),]
nodelist <- unique(c(unique(df1$tail), unique(df1$head)))

dfcord <- data.frame(id = idvec, g1 = geomvec1, g2 = geomvec2)
dfcord$id <- dfcord$id + 1
dfcord.f <- dfcord[dfcord$id %in% nodelist,]

names(df1) <- c("onset","tail", "head", "weight")

m=leaflet(data = dfcord.f) %>% addTiles() %>% addMarkers(~g1, ~g2, popup = ~as.character(id))
m

leaflet(data = dfcord.f) %>% addTiles() %>%
  addWebGLHeatmap(lng=~g1, lat=~g2, intensity = ~id/1000, size=1000)



# Begin analysis

df <- df1
rm(df1)

df$terminus <- df$onset + 600000

df$duration <- 600000
df$onset.censored <- FALSE
df$terminus.censored <- FALSE

df <- df[with(df, order(onset)), ]

df$onset2 <- df$onset
df$onset <- as.POSIXct(df$onset/1000, origin="1970-01-01")
df$edgeid <- paste(df$tail, df$head, sep = "+")

nodevec <- vector()
for (index1 in nodelist){
  for (index2 in nodelist){
    if(index1 == index2){
      
      nodevec <- c(nodevec,paste(index1, index2, sep = "+"))
    }
    
  }
}

node_df <- data.frame(tail_b = nodelist, head_b = nodelist, numerized = 1:length(nodelist))



inv_only_pairs <- df[!df$edgeid %in% nodevec,]

inv_only_pairs <- merge(inv_only_pairs, node_df, by.x = "tail", by.y = "tail_b")
inv_only_pairs$tail <- inv_only_pairs$numerized

inv_only_pairs$head_b <- NULL
inv_only_pairs$numerized <- NULL

inv_only_pairs <- merge(inv_only_pairs, node_df, by.x = "head", by.y = "head_b")
inv_only_pairs$head <- inv_only_pairs$numerized

inv_only_pairs$tail_b <- NULL
inv_only_pairs$numerized <- NULL

inv_only_pairs$edgeid <- paste(inv_only_pairs$tail, inv_only_pairs$head, sep = "+")


nodelist_new <- nodelist
df_vertexattr <- data.frame(vertex.id = nodelist_new)
inv_only_pairs$onset2 <- (as.numeric(inv_only_pairs$onset2) / 600000)*10
inv_only_pairs$onset2 <- inv_only_pairs$onset2  - min(inv_only_pairs$onset2)
inv_only_pairs$duration <- 10
inv_only_pairs$terminus <-inv_only_pairs$onset2 + inv_only_pairs$duration
inv_only_pairs <- inv_only_pairs[with(inv_only_pairs, order(onset2)), ]

#inv_only_pairs_sampled <- inv_only_pairs[inv_only_pairs$onset2 < 92000,]
inv_only_pairs_sampled <- inv_only_pairs

#thenetwork <- network(
#  inv_only_pairs_sampled,
#  vertex.attr = df_vertexattr,
#  vertex.attrnames = c("vertex.id"),
#  directed = F,
#  bipartite = F
#)

df.edges <- data.frame(onset = inv_only_pairs_sampled$onset2,terminus = inv_only_pairs_sampled$terminus,
                       tail.vertex.id = inv_only_pairs_sampled$tail, head.vertex.id = inv_only_pairs_sampled$head)
df.edges$onset <- c(1:length(df.edges$onset))
df.edges$terminus <- df.edges$onset + 10
dynamicCollabs <- networkDynamic(
  # thenetwork,
  edge.spells = df.edges
  #  vertex.spells = df
)

#compute.animation(
#  dynamicCollabs,
#  animation.mode = "kamadakawai",
#  slice.par = list(
#    start = 0,
#    end = length(inv_only_pairs_sampled$terminus),
#    interval = 10,
#    aggregate.dur = 10,
#    rule = "any"
#  )
#)

#render.d3movie(
 # dynamicCollabs,
  #displaylabels = TRUE,
  #filename="milan.html"
  # This slice function makes the labels work
  #vertex.tooltip = function(slice) {
  # paste(
  #  "<b>Name:</b>", (slice %v% "name"),
  # "<br>",
  #"<b>Region:</b>", (slice %v% "region")
  #)
  #}
#)

plot(tEdgeFormation(dynamicCollabs, time.interval = 1
), 
xlab = "Time (minutes)",
ylab = "Edge formation",
main = "Edge formation rate")

dynamicBetweenness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df.edges$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "betweenness"
)
plot(dynamicBetweenness,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

test_dfma <- SMA(dynamicBetweenness,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(dynamicBetweenness,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


dynamicDeg <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df.edges$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "degree"
)
plot(dynamicCent,xlab = "Time (minutes)",
     ylab = "Degree centrality",
     main = "Degree centralization")


dynamicEvcent <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df.edges$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "evcent"
)

dynamicCloseness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df.edges$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "closeness"
)
plot(dynamicCent,xlab = "Time (minutes)",
     ylab = "Degree centrality",
     main = "Degree centralization")
plot(dynamicCloseness,xlab = "Time (minutes)",
     ylab = "Degree centrality",
     main = "Degree centralization")

test_dfma <- SMA(dynamicBetweenness,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(dynamicBetweenness,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(dynamicEvcent,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)






test_dfma <- SMA(inv_only_pairs_sampled$weight,n=3000)
plot(test_dfma, type = 'l')
#test_dfma <- ts(test_dfma, frequency = 8800)
test_dfma <- ts(test_dfma, frequency = 8800)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)

# Begin randomization

df_random <- inv_only_pairs_sampled
df_rndnodes <- unique(c(unique(inv_only_pairs_sampled$tail), unique(inv_only_pairs_sampled$head)))

for (index in c(1:length(inv_only_pairs_sampled$onset2))){
  sampled_A <- sample(df_rndnodes, 1)
  df_random[index,'tail'] <-  sampled_A
  df_random[index,'head'] <-  sample(df_rndnodes[df_rndnodes != sampled_A], 1)
  df_random$edgeid[index] <- paste(df_random$tail[index], df_random$head[index], sep = "+")
  
}

#thenetwork <- network(
#  df_random,
#  vertex.attr = df_vertexattr,
#  vertex.attrnames = c("vertex.id"),
#  directed = F,
#  bipartite = F
#)
#plot(thenetwork)



df.edges_random <- data.frame(onset = df_random$onset2,terminus = df_random$terminus,
                       tail.vertex.id = df_random$tail, head.vertex.id = df_random$head)
dynamicCollabs_random <- networkDynamic(
  # thenetwork,
  edge.spells = df.edges_random
  #  vertex.spells = df
)

#compute.animation(
#  dynamicCollabs_random,
#  animation.mode = "kamadakawai",
#  slice.par = list(
#    start = 0,
#    end = max(df_random$terminus),
#    interval = 10,
#    aggregate.dur = 10,
#    rule = "any"
#  )
#)

dynamicBetweenness_random <- tSnaStats(
  dynamicCollabs_random,
  snafun = "centralization",
  start = 0,
  end = max(df_random$terminus),
  time.interval = 10,
  aggregate.dur = 10,
  FUN = "betweenness"
)

test_dfma_random <- SMA(dynamicBetweenness_random,n=20)
plot(test_dfma_random, type = 'l')
test_dfma_random <- ts(test_dfma_random, frequency = 1000)

test_seriessma_components <- decompose(test_dfma_random)
plot(test_seriessma_components)

plot(dynamicBetweenness_random,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")
