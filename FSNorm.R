
## read comandline args
args = commandArgs(trailingOnly=TRUE)

coverage_track_path = args[1]
capture_track_path = args[2]
FS_track_path = args[3]
chr_sizes_path = args[4]
binsize = args[5]
binsize = as.numeric(binsize)
capture_coefficient_path = args[6]
outdir_path = args[7]

## read capture track
binned_capture_track <- function(capture_track, genome_size, binsize, ic){
  n = 0
  for (i in 1:nrow(capture_track)) {
    if(capture_track[i,2]%/%binsize!=capture_track[i,3]%/%binsize){
      capture_track = rbind(capture_track, c(capture_track[i,1], (capture_track[i,2]%/%binsize + 1)*binsize, capture_track[i,3]))
      capture_track[i,3] = (capture_track[i,2]%/%binsize + 1)*binsize - 1
      n = n + 1
    }
  }
  Z <- matrix(0 , nrow = genome_size, ncol = 1)
  x = capture_track[,2]%/%(binsize) + 1 + ic[capture_track[,1],1]
  y = capture_track[,3] - capture_track[,2]
  for (i in 1:length(x)) {
    if(i%%100000==0){ print(i) }
    Z[x[i],1] = Z[x[i],1] + y[i]
  }
  Z = (binsize - Z)/binsize
  return(Z)
}
## read chr sizes

## genome indexing
chr_s <- read.table(chr_sizes_path, stringsAsFactors = F, head=F)
get_genome_size <- function(chr_s, binsize){
  genome_size = 0
  for (i in 1:24) {
    genome_size = genome_size + chr_s[i,2]%/%binsize + 1
  }
  return(genome_size)
}
genome_size100kb = get_genome_size(chr_s, 100000)
genome_size = get_genome_size(chr_s, binsize)

get_chr_index_correction <- function(chr_s, binsize){
  ic = matrix(0, nrow = 24, ncol = 1)
  ic[1,1] = 0
  for (i in 2:24) {
    ic[i,1] = sum(chr_s[1:(i - 1),2]%/%binsize + 1)
  }
  return(ic)
}
ic100kb = get_chr_index_correction(chr_s, 100000)
ic = get_chr_index_correction(chr_s, binsize)

## get sample coverage
get_coverage <- function(s, genome_size, binsize, ic){
  Z <- matrix(0 , nrow = genome_size, ncol = 1)
  x = s[,2]%/%(binsize) + 1 + ic[s[,1],1]
  for (i in 1:length(x)) {
    if(i%%100000==0){ print(i) }
    Z[x[i],1] = Z[x[i],1] + s[i,4]
  }
  return(Z)
}
## read capture track
capture_track <- read.table(capture_track_path, stringsAsFactors = F, head=F)
capture_track = capture_track[,1:3]
capture_track[which(capture_track[,1]=="chrX"),1] = "chr23"
capture_track[which(capture_track[,1]=="chrY"),1] = "chr24"
capture_track[,1] = as.numeric(substr(capture_track[,1],4,5))
capture_track100kb = binned_capture_track(capture_track, genome_size100kb, 100000, ic100kb)
capture_track = binned_capture_track(capture_track, genome_size, binsize, ic)

capture_coefficient = read.table(capture_coefficient_path, stringsAsFactors = F, head=F)
capture_coefficient = capture_coefficient[1,1]

## coverage
coverage_pre <- read.table(coverage_track_path, stringsAsFactors = F, head=F)
coverage100kb = get_coverage(coverage_pre, genome_size100kb, 100000, ic100kb)
coverage = get_coverage(coverage_pre, genome_size, binsize, ic)
N100kb = length(coverage100kb)
N = length(coverage)


## FS track
# read normalization track
FS_track_pre <- scan(FS_track_path, what="", sep="\n")
FS_track <- strsplit(FS_track_pre, "[[:space:]]+")
names(FS_track) <- sapply(FS_track, `[[`, 1)
FS_track <- lapply(FS_track, `[`, -1)


# capture normalization of coverage100kb
for (i in 1:N100kb) {
  coverage100kb[i] = coverage100kb[i]/(capture_track100kb[i] + capture_coefficient*(1 - capture_track100kb[i]))
}

# eval Norm_coef_100kb
Norm_coef_100kb = rep(1,N100kb)
for (i in 1:N100kb) {
  w = as.numeric(unlist(FS_track[[i]]))
  # w[1] == 0 for indefined FS value
  if(length(w)!=0){
    Norm_coef_100kb[i] = median(coverage100kb[w])
  }
}

# eval Norm_coef_binsize
Norm_coef_binsize = 1:N

bin = 0
bin100kb = 0
n = 100000/binsize
for (i in 1:24) {
  chr_i_size_100kb = chr_s[i,2]%/%100000
  for (bin_local1 in 1:chr_i_size_100kb) {
    bin100kb = bin100kb + 1
    for (bin_local2 in 1:n) {
      bin = bin + 1
      Norm_coef_binsize[bin] = Norm_coef_100kb[bin100kb]/n
    }
  }
  chr_i_size_100kb_rest = (chr_s[i,2] - chr_i_size_100kb*100000)%/%binsize + 1
  bin100kb = bin100kb + 1
  for (j in 1:chr_i_size_100kb_rest) {
    bin = bin + 1
    Norm_coef_binsize[bin] = Norm_coef_100kb[bin100kb]/n
  }
}




## capture normalization
for (i in 1:N) {
  coverage[i] = coverage[i]/(capture_track[i] + capture_coefficient*(1 - capture_track[i]))
}


## FS normalization
# find low or very high coverage bins
#qup = as.numeric(quantile(coverage, 0.9995))
#qdown = as.numeric(quantile(coverage[which(coverage!=0)], 0.05))
wq = which(coverage==0)
#wq = which(coverage<qdown | coverage>qup)
#wq0 = which(coverage==0)
#wq01 = wq0 - 1
#wq01 = wq01[which(wq01>0)]
#wq02 = wq0 + 1
#wq02 = wq02[which(wq02<=N)]
#wq0 = unique(c(wq0,wq01,wq02))
#wq = setdiff(wq,wq0)
# normalization
normalized_coverage = rep(0,N)
for (i in 1:N) {
  if(Norm_coef_binsize[i]!=0){
    normalized_coverage[i] = coverage[i]/Norm_coef_binsize[i]
  }
}
normalized_coverage[wq] = 0.001

## write result to files
ind_i_st = 1
for (i in 1:24) {
  chr_i_bin_count = chr_s[i,2]%/%binsize + 1
  chr_i_coverage = matrix(0, nrow = chr_i_bin_count, ncol = 3)
  
  ind_i_en = ind_i_st + chr_i_bin_count - 1
  chr_i_coverage[,1] = 0:(chr_i_bin_count - 1)*binsize
  chr_i_coverage[,2] = 1:chr_i_bin_count*binsize
  chr_i_coverage[,3] = normalized_coverage[ind_i_st:ind_i_en]
  ind_i_st = ind_i_st + chr_i_bin_count
  write.table(chr_i_coverage, paste(outdir_path, "/chr", as.character(i),"_normalized_coverage", sep = ""), row.names =F, col.names = F, append = F)
}

write.table(normalized_coverage, paste(outdir_path,"/normalized_coverage", sep = ""), row.names =F, col.names = F, append = F)




