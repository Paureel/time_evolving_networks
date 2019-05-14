library(sna)
library(tsna)
library(ndtv)
library(networkDynamic)
library(network)

df <- read.table("ia-contacts_hypertext2009.txt", header = F, stringsAsFactors = F, sep = ",")
names(df) <- c("tail", "head", "onset")
df$terminus <- df$onset + 100
df$duration <- 100
df$onset.censored <- FALSE
df$terminus.censored <- FALSE
#df$edge.id <- 1:length(df$tail)

nodelist <- unique(c(unique(df$tail), unique(df$head)))

df_vertexattr <- data.frame(vertex.id = nodelist)
  
thenetwork <- network(
  df,
  vertex.attr = df_vertexattr,
  vertex.attrnames = c("vertex.id"),
  directed = FALSE,
  bipartite = FALSE
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
    start = 20,
    end = 50000,
    interval = 2000,
    aggregate.dur = 2000,
    rule = "any"
  )
)

render.d3movie(
  dynamicCollabs,
  displaylabels = TRUE,
  filename="conference.html"
  # This slice function makes the labels work
  #vertex.tooltip = function(slice) {
   # paste(
    #  "<b>Name:</b>", (slice %v% "name"),
     # "<br>",
      #"<b>Region:</b>", (slice %v% "region")
    #)
  #}
)

plot(tEdgeFormation(dynamicCollabs, time.interval = 1000), 
     xlab = "Time (seconds)",
     ylab = "Edge formation",
     main = "Edge formation rate")

dynamicBetweenness <- tSnaStats(
  dynamicCollabs,
  snafun = "centralization",
  start = 20,
  end = 212360,
  time.interval = 1000,
  aggregate.dur = 2000,
  FUN = "betweenness"
)
plot(dynamicBetweenness,xlab = "Time (seconds)",
     ylab = "Betweenness",
     main = "Betweenness centralization")
