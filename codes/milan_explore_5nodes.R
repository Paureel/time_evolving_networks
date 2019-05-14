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

set.seed(42)

df1 <- read.table("output-file_5", header = F, 
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

dfcord <- data.frame(id = idvec, g1 = geomvec1, g2 = geomvec2)



#df2 <- read.table("MItoMI-2013-12-02.txt", header = F, 
 #                stringsAsFactors = F, sep = "\t")

#df3 <- read.table("MItoMI-2013-12-03.txt", header = F, 
#                   stringsAsFactors = F, sep = "\t")
  
df <- df1
#datevec <- as.POSIXct(df1$onset, origin="1970-01-01")
names(df) <- c("onset","tail", "head", "weight")
names(df1) <- c("onset","tail", "head", "weight")
nodelist <- unique(c(unique(df$tail), unique(df$head)))

#rnodes <- sample(df$tail, 50)
rnodes <- nodes_selected$V1
df <- df[(df$tail %in% rnodes) & (df$head %in% rnodes),]

#df <- df[(df$tail %in% nodelist[1:3]) & (df$head %in% nodelist[1:3]),]

df2 <- df
df2$edgeid <- paste(df2$tail, df2$head, sep = "+")

edgelist <- unique(df2$edgeid)
#df2 <- df2[df2$edgeid %in% edgelist[1:50],]

df2$tail <- NULL
df2$head <- NULL

ggplot(df2, aes(x = onset, y = weight)) + 
  geom_line(aes(color = edgeid), size = 1)  +
  theme_minimal()
df2s <- spread(df2, key = edgeid, value = weight)
rownames(df2s) <- df2s$onset
df2s$onset <- NULL

cormat <- cor(df2s, use = "pairwise.complete.obs")
dev.new()
heatmap.2(cormat)

df$edgeid <- paste(df$tail, df$head, sep = "+")

#df <- df[df$onset,]

df.meanagg <- aggregate(weight ~ edgeid, data=df, mean)
p <- plot_ly(alpha = 0.6, data = df.meanagg) %>% add_histogram(x = ~weight)
p



nodelist <- unique(c(unique(df$tail), unique(df$head)))

#df$onset <- df$onset / 60000
#df$onset <- df$onset  - min(df$onset)



df$terminus <- df$onset + 600000

df$duration <- 600000
df$onset.censored <- FALSE
df$terminus.censored <- FALSE

df <- df[with(df, order(onset)), ]
dfcord$id <- dfcord$id + 1

dfcord.f <- dfcord[dfcord$id %in% nodelist,]

m=leaflet(data = dfcord.f) %>% addTiles() %>% addMarkers(~g1, ~g2, popup = ~as.character(id))
m

leaflet(data = dfcord.f) %>% addProviderTiles(providers$CartoDB.DarkMatter) %>%
  addWebGLHeatmap(lng=~g1, lat=~g2, intensity = ~id/100, size=1000)


i <- 1

#for (index in unique(df$edgeid)){
 # df.temp <- df[df$edgeid == index,]
  
#}
df$onset2 <- df$onset
df$onset <- as.POSIXct(df$onset/1000, origin="1970-01-01")


only_pairs <- df[df$edgeid %in% c("1008+1008", "3785+3785", "4507+4507", 
                                  "5229+5229", "8069+8069"),]

inv_only_pairs <- df[!df$edgeid %in% c("1008+1008", "3785+3785", "4507+4507", 
                                      "5229+5229", "8069+8069"),]

df <- df %>% group_by(edgeid) %>% mutate(Nor = (weight-mean(weight))/sd(weight)) 
only_pairs <- only_pairs %>% group_by(edgeid) %>% mutate(Nor = (weight-mean(weight))/sd(weight)) 
inv_only_pairs <- inv_only_pairs %>% group_by(edgeid) %>% mutate(Nor = (weight-mean(weight))/sd(weight)) 



ggplotly(ggplot(df[1:2500,], aes(x=onset, y=Nor, color=edgeid)) +
  geom_line())

ggplotly(ggplot(only_pairs, aes(x=onset, y=Nor, color=edgeid)) +
           geom_line() +labs(
                             x = "Date",
                             y = "Normalized call intensity",
                             title ="Normalized call intensity of Milan provinces"))

ggplotly(ggplot(inv_only_pairs, aes(x=onset, y=Nor, color=edgeid)) +
           geom_line() +labs(
             x = "Date",
             y = "Normalized call intensity",
             title ="Normalized call intensity between Milan provinces"))

#test_series <- df[1:2500,]
test_series <- df
edgsel <- "3785+3785"
test_series <- test_series[test_series$edgeid == edgsel,]



test_seriessma <- SMA(test_series$weight,n=20)
test_seriessma<-ts(test_seriessma, frequency=60)
plot.ts(test_seriessma)
test_seriessma_components <- decompose(test_seriessma)
plot(test_seriessma_components)

acfplot <- Acf(test_seriessma)
plot(acfplot, main = edgsel)

#milano.model <- Arima(test_seriessma,order=c(100,1,1),lambda=0)
#plot(forecast(milano.model,h=200))
error

#df <- df[df$weight > median(df$weight),]

df <- inv_only_pairs
#df <- df[df$edgeid %in% inv_only_pairs,]

df$tail[df$tail == 1008] <- 1
df$tail[df$tail == 3785] <- 2
df$tail[df$tail == 4507] <- 3
df$tail[df$tail == 5229] <- 4
df$tail[df$tail == 8069] <- 5

df$head[df$head == 1008] <- 1
df$head[df$head == 3785] <- 2
df$head[df$head == 4507] <- 3
df$head[df$head == 5229] <- 4
df$head[df$head == 8069] <- 5

nodelist_new <- nodelist

nodelist_new[nodelist_new == 1008] <- 1
nodelist_new[nodelist_new == 3785] <- 2
nodelist_new[nodelist_new== 4507] <- 3
nodelist_new[nodelist_new == 5229] <- 4
nodelist_new[nodelist_new == 8069] <- 5


df_vertexattr <- data.frame(vertex.id = nodelist_new)

df$onset <- (as.numeric(df$onset) / 600000)*10
df$onset <- df$onset  - min(df$onset)
df$duration <- 10
df$terminus <-df$onset + df$duration

df$onset <- 1:(length(df$onset))
df$terminus <-df$onset + df$duration
####

test_dfma <- SMA(df$weight,n=20)
#plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 23)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


####





df.edges <- data.frame(onset = df$onset,terminus = df$terminus,
                       tail.vertex.id = df$tail, head.vertex.id = df$head)

thenetwork <- network(
  df.edges[,c("tail.vertex.id", "head.vertex.id")],
  vertex.attr = df_vertexattr,
  vertex.attrnames = c("vertex.id"),
  directed = F,
  bipartite = F
)
#plot(thenetwork)




dynamicCollabs <- networkDynamic(
  # thenetwork,
  edge.spells = df.edges,
  directed = F
  #  vertex.spells = df
)

compute.animation(
  dynamicCollabs,
  animation.mode = "kamadakawai",
  slice.par = list(
    start = 1,
    end = length(df$onset),
    interval = 1,
    aggregate.dur = 1,
    rule = "any"
  )
)

render.d3movie(
  dynamicCollabs,
  displaylabels = TRUE,
  filename="milan.html"
  # This slice function makes the labels work
  #vertex.tooltip = function(slice) {
  # paste(
  #  "<b>Name:</b>", (slice %v% "name"),
  # "<br>",
  #"<b>Region:</b>", (slice %v% "region")
  #)
  #}
)

#plot(tEdgeFormation(dynamicCollabs, time.interval = 1), 
    # xlab = "Time (minutes)",
    # ylab = "Edge formation",
   #  main = "Edge formation rate")

plot(tEdgeFormation(dynamicCollabs, time.interval = 1
                    ), 
      xlab = "Time (minutes)",
      ylab = "Edge formation",
       main = "Edge formation rate")

dynamicBetweenness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "betweenness", mode = "graph"
)

dynamicCloseness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "closeness", mode = "graph"
)

dynamicDegree <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "degree"
)
ggplotly(qplot(c(1:length(dynamicDegree)),dynamicDegree))

dynamicEvcent <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 1,
  end = max(df$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "evcent", mode = "graph"
)



plot(dynamicBetweenness,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicCloseness,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicDegree,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicEvcent,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")




test_dfma <- SMA(dynamicBetweenness,n=20)
plot(test_dfma, type = 'l')
dynbet_backup <- test_dfma
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)

acfplot <- Acf(test_dfma,lag.max = 120)


test_dfma <- SMA(dynamicCloseness,n=20)
plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(dynamicDegree,n=20)
plot(test_dfma, type = 'l')
dyndeg_backup <- test_dfma
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(na.remove(dynamicEvcent),n=20)
plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)

proximity.timeline(dynamicCollabs,mode='sammon',default.dist=10,start=0,end=200,
                   labels.at=c(1,16,25),label.cex=0.7,vertex.col=c(rep('gray',14),'green','blue') )



randoms = data.frame(matrix(nrow = length(dynamicBetweenness)))
randoms_deg <- data.frame(matrix(nrow = length(dynamicBetweenness)))

for (i in c(1:50)){
cat(i)
df_random <- df
df_rndnodes <- unique(c(unique(df$tail), unique(df$head)))

for (index in c(1:length(df$onset))){
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



df.edges <- data.frame(onset = df_random$onset,terminus = df_random$terminus,
                       tail.vertex.id = df_random$tail, head.vertex.id = df_random$head)
dynamicCollabs_random <- networkDynamic(
  # thenetwork,
  edge.spells = df.edges,
  directed = F
  #  vertex.spells = df
)

#compute.animation(
#  dynamicCollabs_random,
#  animation.mode = "kamadakawai",
#  slice.par = list(
#    start = 0,
#    end = 900,
#    interval = 1,
#    aggregate.dur = 1,
#    rule = "any"
#  )
#)

#render.d3movie(
 # dynamicCollabs_random,
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


dynamicBetweenness_random <- tSnaStats(
  dynamicCollabs_random,
  snafun = "centralization",
  start = 1,
  end = max(df_random$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "betweenness", mode = "graph"
)

dynamicCloseness_random <- tSnaStats(
  dynamicCollabs_random,
  snafun = "centralization",
  start = 1,
  end = max(df_random$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "closeness", mode = "graph"
)

dynamicDegree_random <- tSnaStats(
  dynamicCollabs_random,
  snafun = "centralization",
  start = 1,
  end = max(df_random$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "degree", mode = "graph"
)

dynamicEvcent_random <- tSnaStats(
  dynamicCollabs_random,
  snafun = "centralization",
  start = 1,
  end = max(df_random$onset),
  time.interval = 1,
  aggregate.dur = 1,
  FUN = "evcent", mode = "graph"
)



plot(dynamicBetweenness_random,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicCloseness_random,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicDegree_random,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

plot(dynamicEvcent_random,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")




test_dfma <- SMA(dynamicBetweenness_random,n=20)
randoms <- cbind(randoms, as.numeric(test_dfma))

plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)



acfplot <- Acf(test_dfma,lag.max = 120)

test_dfma <- SMA(dynamicCloseness_random,n=20)
plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(dynamicDegree_random,n=20)
randoms_deg <- cbind(randoms_deg, as.numeric(test_dfma))
plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)


test_dfma <- SMA(na.remove(dynamicEvcent_random),n=20)
plot(test_dfma, type = 'l')
test_dfma <- ts(test_dfma, frequency = 15)

test_seriessma_components <- decompose(test_dfma)
plot(test_seriessma_components)
}
## 1##
rand_back <- randoms
rand_back_deg <- randoms_deg

names(randoms) <- 1:length(names(randoms))
randoms$`1` <- 1:length(randoms[,2])
meansrand <- rowMeans(randoms[,-1])
meansmed <- rowMedians(as.matrix(randoms[,-1]))
sds <- rowSds(as.matrix(randoms[,-1]))
firstq <- rowQuantiles(as.matrix(randoms[,-1]), probs = 0.25)
secondq <- rowQuantiles(as.matrix(randoms[,-1]), probs = 0.75)
randoms$meansrandcol  <- meansrand
randoms$meansmedcol  <- meansmed
randoms$sds <- sds
randoms$firstq <- firstq
randoms$secondq <- secondq

randoms$plussds_betC <- randoms$meansrandcol + randoms$sds
randoms$minussds_betC <- randoms$meansrandcol - randoms$sds

randoms$firstq_betC <- randoms$meansmedcol + randoms$firstq
randoms$secondq_betC  <- randoms$meansmedcol - randoms$secondq

#todo randoms <- cbind(randoms,as.numeric(test_dfma))
randoms <-  cbind(randoms,as.numeric(dynbet_backup))
#names(randoms)[length(names(randoms))] <- "True_TS_betC"
melted_randoms <- melt(randoms, id.vars = c('1', 'plussds_betC','minussds_betC'))
melted_randoms <- melted_randoms[complete.cases(melted_randoms),]
melted_randoms$variable <- as.character(melted_randoms$variable)
names(melted_randoms)[1] <- "Time"
melted_randoms$variable[melted_randoms$variable == "as.numeric(dynbet_backup)"] <- "True_TS_betC"
melted_randoms$variable[melted_randoms$variable == "meansrandcol"] <- "Mean_of_randoms_betC"

melted_randoms_q <- melt(randoms, id.vars = c('1', 'firstq_betC','secondq_betC'))
melted_randoms_q <- melted_randoms_q[complete.cases(melted_randoms_q),]
melted_randoms_q$variable <- as.character(melted_randoms_q$variable)
names(melted_randoms_q)[1] <- "Time"
melted_randoms_q$variable[melted_randoms_q$variable == "as.numeric(dynbet_backup)"] <- "True_TS_betC"
melted_randoms_q$variable[melted_randoms_q$variable == "meansmedcol"] <- "Median_of_randoms_betC"

## 2 ##
#rand_back <- randoms_deg
names(randoms_deg) <- 1:length(names(randoms_deg))
randoms_deg$`1` <- 1:length(randoms_deg[,2])
meansrand <- rowMeans(randoms_deg[,-1])
meansmed <- rowMedians(as.matrix(randoms_deg[,-1]))
sds <- rowSds(as.matrix(randoms_deg[,-1]))
firstq <- rowQuantiles(as.matrix(randoms_deg[,-1]), probs = 0.25)
secondq <- rowQuantiles(as.matrix(randoms_deg[,-1]), probs = 0.75)
randoms_deg$firstq <- firstq
randoms_deg$secondq <- secondq

randoms_deg$meansrandcol  <- meansrand
randoms_deg$meansmedcol  <- meansmed
randoms_deg$sds <- sds
randoms_deg$plussds_degreeC <- randoms_deg$meansrandcol + randoms_deg$sds
randoms_deg$minussds_degreeC <- randoms_deg$meansrandcol - randoms_deg$sds
randoms_deg$firstq_degreeC <- randoms_deg$meansmedcol + randoms_deg$firstq
randoms_deg$secondq_degreeC  <- randoms_deg$meansmedcol - randoms_deg$secondq

#todo randoms <- cbind(randoms,as.numeric(test_dfma))
randoms_deg <-  cbind(randoms_deg,as.numeric(dyndeg_backup))
#names(randoms_deg)[length(names(randoms_deg))] <- "True_TS_betC"
melted_randoms_deg <- melt(randoms_deg, id.vars = c('1', 'plussds_degreeC','minussds_degreeC'))
melted_randoms_deg <- melted_randoms_deg[complete.cases(melted_randoms_deg),]
melted_randoms_deg$variable <- as.character(melted_randoms_deg$variable)
names(melted_randoms_deg)[1] <- "Time"
melted_randoms_deg$variable[melted_randoms_deg$variable == "as.numeric(dyndeg_backup)"] <- "True_TS_degreeC"
melted_randoms_deg$variable[melted_randoms_deg$variable == "meansrandcol"] <- "Mean_of_randoms_degreeC"

melted_randoms_deg_q <- melt(randoms_deg, id.vars = c('1', 'firstq_degreeC','secondq_degreeC'))
melted_randoms_deg_q <- melted_randoms_deg_q[complete.cases(melted_randoms_deg_q),]
melted_randoms_deg_q$variable <- as.character(melted_randoms_deg_q$variable)
names(melted_randoms_deg_q)[1] <- "Time"
melted_randoms_deg_q$variable[melted_randoms_deg_q$variable == "as.numeric(dyndeg_backup)"] <- "True_TS_degreeC"
melted_randoms_deg_q$variable[melted_randoms_deg_q$variable == "meansmedcol"] <- "Median_of_randoms_degreeC"



error2


g_bet <- ggplotly(ggplot(melted_randoms[melted_randoms$variable == 'True_TS_betC' | melted_randoms$variable == 'Mean_of_randoms_betC',], aes(x=Time, y=value, color=variable, 
                                    ymin=minussds_betC, ymax=plussds_betC)) +
  geom_line() + geom_ribbon(alpha=0.5) + ylab("Betweenness centrality") + guides(fill=guide_legend(title="Legend"))
  )
g_deg <- ggplotly(ggplot(melted_randoms_deg[melted_randoms_deg$variable == 'True_TS_degreeC' | melted_randoms_deg$variable == 'Mean_of_randoms_degreeC',], aes(x=Time, y=value, color=variable, 
                                                                                                                                   ymin=minussds_degreeC, ymax=plussds_degreeC)) +
                    geom_line() + geom_ribbon(alpha=0.5) + ylab("Degreee centrality") + guides(fill=guide_legend(title="Legend"))
)

g_bet_q <- ggplotly(ggplot(melted_randoms_q[melted_randoms_q$variable == 'True_TS_betC' | melted_randoms_q$variable == 'Median_of_randoms_betC',], aes(x=Time, y=value, color=variable, 
                                                                                                                                             ymin=firstq_betC, ymax=secondq_betC)) +
                    geom_line() + geom_ribbon(alpha=0.5) + ylab("Betweenness centrality") + guides(fill=guide_legend(title="Legend"))
)

g_deg_q <- ggplotly(ggplot(melted_randoms_deg_q[melted_randoms_deg_q$variable == 'True_TS_degreeC' | melted_randoms_deg_q$variable == 'Median_of_randoms_degreeC',], aes(x=Time, y=value, color=variable, 
                                                                                                                                                               ymin=firstq_degreeC, ymax=secondq_degreeC)) +
                    geom_line() + geom_ribbon(alpha=0.5) + ylab("Degreee centrality") + guides(fill=guide_legend(title="Legend"))
)

p <- subplot(g_bet, g_deg,g_bet_q, g_deg_q)

htmlwidgets::saveWidget(as.widget(p), "index_combined.html")
melted_randoms_backup <- melted_randoms


newframe <- data.frame(Time = c(1:length(randoms$`1`)), true_betweenness = randoms$`as.numeric(dynbet_backup)`,
                       true_degree = randoms_deg$`as.numeric(dyndeg_backup)`,
                       random_mean_betweenness = randoms$meansrandcol, 
                       random_mean_degree = randoms_deg$meansrandcol)
newframe_melted <- melt(newframe, id.vars = "Time")


f <- ggplotly(ggplot(newframe_melted, aes(x=Time, y=value, color=variable, 
                        )) +
                    geom_line()  + ylab("Centrality") + guides(fill=guide_legend(title="Legend"))
)
f
htmlwidgets::saveWidget(as.widget(f), "centralities.html")

#histograms
hist(as.numeric(randoms[100,c(2:51)]), main = "Histogram of the random values of the betweenness centrality at timestep = 100", xlab = "Betweenness centrality")

#randoms_deg$sds <- sds
#randoms$meansrandcol  <- meansrand
#randoms$meansmedcol  <- meansmed
ggplotly(qplot(randoms$meansrandcol, randoms$meansmedcol) + 
           xlab("Mean of the betweenness centrality")+
           ylab("Median of the betweenness centrality")
           )
ggplotly(qplot(randoms_deg$sds, randoms$firstq) + 
           xlab("SD of the betweenness centrality")+
           ylab("0.25 quantile of the betweenness centrality")
)
ggplotly(qplot(randoms_deg$sds, randoms$secondq) + 
           xlab("SD of the betweenness centrality")+
           ylab("0.75 quantile of the betweenness centrality")
)
ggplotly(qplot(randoms_deg$sds, randoms$secondq-randoms$firstq) + 
           xlab("SD of the betweenness centrality")+
           ylab("0.75 - 0.25 quantile of the betweenness centrality")
)

meantrues <- mean(melted_randoms[melted_randoms$variable == 'True_TS_degreeC',]$value)
meanrandoms <- mean(melted_randoms[melted_randoms$variable == 'Mean_of_randoms_degreeC',]$value)

c2 <- ggplotly(qplot(c(1:length(dynamicBetweenness)),dynamicBetweenness)+
                 geom_line() )
c1 <- ggplotly(qplot(c(1:length(dynamicDegree)),dynamicDegree)+
                 geom_line())
c4 <- ggplotly(qplot(c(1:length(dynamicDegree)),(dynamicDegree-dynamicBetweenness))+
                 geom_line())
c3 <- subplot(c1, c2, c4,titleX = TRUE, titleY = TRUE) 
c3 

#fig['layout'].update(height=600, width=800, title='i <3 annotations and subplots')
proximity.timeline(dynamicCollabs_random,mode='sammon',default.dist=10,start=0,end=200,
                   labels.at=c(1,16,25),label.cex=0.7,vertex.col=c(rep('gray',14),'green','blue') )



data(toy_epi_sim)
# set up layout to draw plots under timeline
#layout(matrix(c(1,1,1,2,3,4),nrow=2,ncol=3,byrow=TRUE))
# plot a proximity.timeline illustrating infection spread
proximity.timeline(dynamicCollabs,vertex.col = 'ndtvcol',spline.style='color.attribute',mode = 'sammon',default.dist=100,chain.direction='reverse')
# plot 3 static cross-sectional networks
# (beginning, middle and end) underneath for comparison
plot(network.collapse(dynamicCollabs,at=1),vertex.col='ndtvcol',vertex.cex=2,main='toy_epi_sim network at t=1')
plot(network.collapse(dynamicCollabs,at=17),vertex.col='ndtvcol',vertex.cex=2,main='toy_epi_sim network at=17')
plot(network.collapse(dynamicCollabs,at=25),vertex.col='ndtvcol',
vertex.cex=2,main='toy_epi_sim network at t=25')
layout(1)


proximity.timeline(dynamicCollabs,start=0,end=100,mode='sammon',spline.style='inactive.ignore')

times <- df1$V1
times_conv <- as.POSIXct(times/1000, origin="1970-01-01")
