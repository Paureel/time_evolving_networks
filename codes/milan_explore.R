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

set.seed(42)

df1 <- read.table("MItoMI-2013-12-01.txt", header = F, 
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

rnodes <- sample(df$tail, 50)

df <- df[(df$tail %in% rnodes) & (df$head %in% rnodes),]

#df <- df[(df$tail %in% nodelist[1:3]) & (df$head %in% nodelist[1:3]),]

df2 <- df1
df2$edgeid <- paste(df2$tail, df2$head, sep = "+")
edgelist <- unique(df2$edgeid)
df2 <- df2[df2$edgeid %in% edgelist[1:50],]

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

df$onset <- df$onset / 60000
df$onset <- df$onset  - min(df$onset)



df$terminus <- df$onset + 10

df$duration <- 10
df$onset.censored <- FALSE
df$terminus.censored <- FALSE

df <- df[with(df, order(onset)), ]
dfcord$id <- dfcord$id + 1

dfcord.f <- dfcord[dfcord$id %in% nodelist,]

m=leaflet(data = dfcord.f) %>% addTiles() %>% addMarkers(~g1, ~g2, popup = ~as.character(id))
m

i <- 1

#for (index in unique(df$edgeid)){
 # df.temp <- df[df$edgeid == index,]
  
#}

df.temp1 <- df[df$edgeid == unique(df$edgeid)[1],]
df.temp2 <- df[df$edgeid == unique(df$edgeid)[2],]

plot(df.temp1$onset, df.temp1$weight, type = 'l')
plot(df.temp2$onset, df.temp2$weight, type = 'l')



df <- df[df$weight > median(df$weight),]

df_vertexattr <- data.frame(vertex.id = nodelist)

thenetwork <- network(
  df,
  vertex.attr = df_vertexattr,
  vertex.attrnames = c("vertex.id"),
  directed = F,
  bipartite = T
)
#plot(thenetwork)


df.edges <- data.frame(onset = df$onset,terminus = df$terminus,
                       tail.vertex.id = df$tail, head.vertex.id = df$head)
dynamicCollabs <- networkDynamic(
  # thenetwork,
  edge.spells = df.edges
  #  vertex.spells = df
)

compute.animation(
  dynamicCollabs,
  animation.mode = "kamadakawai",
  slice.par = list(
    start = 0,
    end = max(df$terminus),
    interval = 50,
    aggregate.dur = 50,
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

plot(tEdgeFormation(dynamicCollabs, time.interval = 50), 
     xlab = "Time (minutes)",
     ylab = "Edge formation",
     main = "Edge formation rate")


dynamicBetweenness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 0,
  end = max(df$terminus),
  time.interval = 200,
  aggregate.dur = 200,
  FUN = "betweenness"
)
plot(dynamicBetweenness,xlab = "Time (minutes)",
     ylab = "Betweenness",
     main = "Betweenness centralization")

proximity.timeline(dynamicCollabs,mode='sammon',default.dist=10,start=0,end=100,
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


proximity.timeline(dynamicCollabs,start=0,end=1000,mode='sammon',spline.style='inactive.ignore')
