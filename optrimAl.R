# optrimAl.R - run this second

cutoff_trim <- readLines('cutoff_trim.txt')

# create one multiple tables for each threshold value to store AMAS results

#### Maybe this is the first cut off trim i.e. zero which also sets up the  table 
### YES - THIS IS CORRETC cutoff must include 0 and 1
amas_table <- read.table('summary_0.txt', header = TRUE)
sites <- data.frame(row.names = amas_table$Alignment_name)
pct <- data.frame(row.names = amas_table$Alignment_name)
filled <- data.frame(row.names = amas_table$Alignment_name)
lost <- data.frame(row.names = amas_table$Alignment_name)

for(i in 1:length(cutoff_trim)){
  amas_table <- read.table(paste('summary_', cutoff_trim[i], '.txt', sep = ''), header = TRUE)
  for(j in amas_table$Alignment_name){
    sites[rownames(sites) == j,i] <- amas_table$Parsimony_informative_sites[amas_table$Alignment_name == j]
    pct[rownames(pct) == j,i] <- as.numeric(amas_table$Proportion_parsimony_informative[amas_table$Alignment_name == j])
    filled[rownames(filled) == j,i] <- amas_table$Total_matrix_cells[amas_table$Alignment_name == j] * (1 - amas_table$Missing_percent[amas_table$Alignment_name == j] / 100)
  }
}

# calculate data loss for each trimming threshold

sites[is.na(sites)] <- 0
pct[is.na(pct)] <- 0

for(i in 1:ncol(filled)){
  lost[,i] <- 1 - filled[,i] / filled[,1]
}

lost[is.na(lost)] <- 1

colnames(sites) <- cutoff_trim
colnames(pct) <- cutoff_trim
colnames(filled) <- cutoff_trim
colnames(lost) <- cutoff_trim

# select optimal trimming threshold
# current criterion is maximum proportion of parsimony informative sites where data loss is no more than one median absolute deviation above the median

optrim <- numeric()
optrim_loss <- numeric()

for(i in rownames(pct)){
  lost_i <- unlist(lost[rownames(lost) == i, ])
  pct_i <- unlist(pct[rownames(pct) == i, ])
  dldp <- data.frame(pct_i, lost_i, row.names = cutoff_trim)
  write.csv(dldp, paste('dldp_', i, '.csv', sep = ''))
  real_loss <- dldp$lost_i[dldp$lost_i < 1]
  diff_loss <- real_loss[2:length(real_loss)] - real_loss[1:(length(real_loss) - 1)]
  median_loss <- median(diff_loss[diff_loss != 0])
  dldp <- subset(dldp, dldp$lost_i <= (median(real_loss) + median_loss))
  if(length(dldp$pct_i) > 0){
    optrim[i] <- rownames(dldp)[dldp$pct_i == max(dldp$pct_i)][[1]]
    optrim_loss[i] <- dldp$lost_i[rownames(dldp) == optrim[i][[1]]]
  } else {
    optrim[i] <- 0
    optrim_loss[i] <- 0
  }
}

# generate graphs to show effect of trimming on informativeness and data loss

for(i in rownames(pct)){
  dldp <- read.csv(paste('dldp_', i, '.csv', sep = ''))
  png(paste('dldp_', i, '.png', sep = ''))
  par(mar = c(5,5,2,5))
  plot(main = i, dldp$lost_i ~ cutoff_trim, ylim = c(0,1), ylab = 'proportion of data lost', xlab = 'strictness of trimming (trimAl gap threshold)', pch = 18, col = 'red')
  par(new = T)
  plot(dldp$pct_i ~ cutoff_trim, xlab = NA, ylab = NA, ylim = c(0,1), axes = F, pch = 16, col = 'blue')
  axis(side = 4)
  mtext(side = 4, line = 3, 'proportion parsimony informative')
  legend(x = 0, y = 1, legend = c('proportion of data lost', 'proportion of parsimony informative sites', 'selected trimming threshold'), pch = c(18, 16, NA), lty = c(NA, NA, 2), col = c('red', 'blue', 'black'), cex = 0.9, bty = 'n')
  if(is.na(optrim[i]) == FALSE){
    lines(c(-0.5, optrim[i]), c(optrim_loss[i], optrim_loss[i]), lty = 2)
    lines(c(-0.5, optrim[i]), c(dldp$pct_i[dldp$X == optrim[i]], dldp$pct_i[dldp$X == optrim[i]]), lty = 2)
    lines(c(optrim[i], optrim[i]), c(-0.5, max(optrim_loss[i], dldp$pct_i[dldp$X == optrim[i]])), lty = 2)
  }
  dev.off()
}

overlost <- names(optrim_loss[optrim_loss > 0.3])

write(overlost, 'overlost.txt', sep = '\n')

# Paul B. - copies the chosen alignment file to the pwd optrimal dir:
#file.copy(paste(optrim, '/', names(optrim), sep = ''), getwd()) - Paul B. changed to overwrite file for pipeline
file.copy(paste(optrim, '/', names(optrim), sep = ''), getwd(), overwrite = TRUE)

file.remove(paste(overlost, sep = '')) 