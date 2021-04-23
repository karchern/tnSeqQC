library(tidyverse)
#library(ggembl)
library(patchwork)
library(ggrepel)

argv <- commandArgs(trailingOnly = TRUE)

data <- read_tsv(argv[1], col_names = F)
colnames(data) <- c("sampleID", "repeatInfo", "MapPos")
data$mapped <- ifelse(data$MapPos != 0, " and read mapped", " and read not mapped")
data$output <- str_c(data$repeatInfo, data$mapped)
data$output <- as.factor(data$output)
data2 <- data %>% group_by(sampleID) %>% do(., {
  tmp <- .
  tbl <- table(tmp$output)
  as.data.frame(tbl/sum(tbl))
  #as.data.frame(tbl)
})

data3 <- data2 %>% filter(Var1 == "repeat not found and barcode replacement sequence not found and read not mapped") %>% arrange(Freq) %>% pull(sampleID)
l <- data3
data2$sampleID <- factor(as.vector(data2$sampleID), levels  = data3)
p <- ggplot(data2, aes(x = sampleID, y = Freq, fill = Var1)) + geom_bar(stat = 'identity')
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + ylab("Proportion")

data4 <- data %>% group_by(sampleID) %>% summarize(n = length(repeatInfo))
data4$sampleID <- factor(as.vector(data4$sampleID), levels = data3)
p2 <- ggplot(data4, aes(x = sampleID, y = n)) + geom_bar(stat = 'identity')
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- p2 + theme(axis.text.x = element_blank())
p2 <- p2 + ylab("Total number of forward reads")

pfin <- p2/p

#ggsave(filename = argv[4], width = 16, height = 9, dpi = 300)

#########
### 2 ###
#########

data <- read_tsv(argv[3], col_names = F)
#data$X2 <- NULL
colnames(data) <- c("sampleID", "position", "depth")
#data <- data %>% group_by(sampleID) %>% arrange(positions)

#print(head(data))
data2 <- data %>% filter(depth != 0) %>% group_by(sampleID) %>% nest() %>% mutate(hists = map(data, function(x){
  p <- ggplot(x %>% arrange(position), aes(x = position, y = depth, group = 1)) + geom_line()
  #p <- p + theme_embl()
  p <- p + xlim(-40000, 4660999)
  return(p)
}))

data <- read_tsv(argv[2], col_names = F)
colnames(data) <- c("sampleID", "readID", "readSequence", "position", "barcodeSequence")

# percentage0Positions
data <- data %>% group_by(sampleID) %>% nest() %>% mutate(percentage0Positions = map_dbl(data, function(x) {
  dims <- dim(x)[1]
  sum0s <- x %>% pull(position)
  sum0s <- sum(sum0s == 0)
  return(sum0s/dims)
}))

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

ModeFreq <- function(x) {
  ta <- table(x)
  return(ta[which.max(ta)]/sum(ta))
}

data3 <- data %>% mutate(rollingWindowBarcodeData = map(data, function(x){
  minM <- min(x$position) - 11000
  maxM <- max(x$position) + 10000
  seqq <- seq(minM, maxM, 5000)
  x$cutt <- cut(x$position, seqq)
  x <- x %>% group_by(cutt) %>% summarize(numReads = length(cutt),
                                          modeBarCode = Mode(barcodeSequence),
                                          modeFreq = ModeFreq(barcodeSequence))
  return(x)
}))

barcodeInfo <- data3 %>% select(sampleID, rollingWindowBarcodeData) %>% unnest()
barcodeInfo <- barcodeInfo %>% group_by(sampleID) %>% nest()
barcodeInfo$sampleID <- map_chr(barcodeInfo$sampleID, function(x) str_split(string = x, pattern = "_1_sequence")[[1]][1])
#print(head(data2$sampleID))
#print(head(barcodeInfo$sampleID))
#saddss
data2 <- left_join(data2, barcodeInfo, by = 'sampleID')
#print(head(data2))
data2 <- data2 %>% mutate(barcodeBarplots = map(data.y, function(x){
  if (is.null(x)){
    return(NULL)
  }
  #print(x)
  x$midPoint <- map_dbl(x$cutt, function(y){
    #print(y)
    a <- str_split(string = y, pattern = ",")[[1]][[1]]
    #print(a)
    a <- as.numeric(str_split(string = a, pattern = '\\(')[[1]][2])
    #print(a)
    b <- str_split(string = y, pattern = ",")[[1]][[2]]
    #print(b)
    b <- as.numeric(gsub(pattern = "\\]", replacement = "", x = b))
    #print(b)
    #print(c(a, b))
    return(mean(c(a, b)))
  })
  #print('what')
  #print(x)
  #return(x)
  x <- x %>% select(midPoint, numReads, modeBarCode, modeFreq)
  p <- ggplot() + geom_bar(data = x, aes(x = midPoint, y = numReads, fill = modeFreq), stat='identity', width = 15000)
  p <- p + geom_text_repel(data = x %>% filter(numReads > 10), aes(x = midPoint, y = numReads, label = modeBarCode), min.segment.length = 0)
  #p <- p + theme_embl()
  p <- p + xlim(-40000, 4660999)
  return(p)
}))

# Important: bars below 0 on the x-axis correspond to UNMAPPED reads with a barcode detected.
data2 <- data2 %>% mutate(finPlots = map2(hists, barcodeBarplots, function(x,y) x/y))
map2(data2$finPlots, data2$sampleID, function(x, y) ggsave(plot = x, filename = str_c("testPlot", y, ".png")))

###
allDepths <- data2 %>% summarize(depth = map_dbl(data.x, function(x){
  return(sum(x$depth))
}))
meanDepth <- mean(allDepths$depth)
data2 <- data2 %>% mutate(data.z = map(data.y, function(x){
  if (is.null(x)){
    return(data.frame(NULL))
  }
  if (sum(x$numReads) < 500){
    return(data.frame(NULL))
  }
  # Keep only those peaks with at least 0.5% coverage of the highest peak
  x <- x %>% mutate(numReadsFrac = numReads/max(numReads)) %>% filter(numReadsFrac > 0.005) %>% mutate(modeBarCode = as.factor(as.numeric(as.factor(modeBarCode))))
  #print(x)
  #asdsad
}))
tmp <- data2 %>% select(sampleID, data.z) %>% unnest()
print(tmp)
print(l)
l <- map_chr(l, function(x) str_split(string = x, pattern = "_1_sequence")[[1]][1])
tmp$sampleID <- factor(as.vector(tmp$sampleID), levels = l)
# print(tmp)
# sadasdasd
p3 <- ggplot(tmp, aes(x = sampleID, y = numReadsFrac, color = modeBarCode)) + geom_point()
p3 <- p3 + theme(axis.text.x = element_blank(),
                 axis.title.x = element_blank())
p3 <- p3 + scale_color_discrete(drop=FALSE) + scale_x_discrete(drop=F)
p2 <- p2 + theme(axis.title.x = element_blank())
pfin <- p3/p2/p
#pfin <- p3
# 
ggsave(filename = argv[4], plot = pfin, width = 16, height = 9, dpi = 300)