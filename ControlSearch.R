args = commandArgs(trailingOnly=TRUE)

outdir_path = args[1]
DIR = args[2]
binsize = args[3]
binsize = as.numeric(binsize)

control_samples_list <- scan(paste(outdir_path,"/control_samples_list", sep = ""), what="", sep="\n")
case_normalized_coverage = read.table(paste(outdir_path,"/normalized_coverage", sep = ""), stringsAsFactors = F, head=F)
chr_sizes = read.table(paste(DIR,"/chr_sizes_hg19", sep = ""), stringsAsFactors = F, head=F)
chr_sizes[,2] = as.numeric(chr_sizes[,2])

N = length(control_samples_list)
vars = matrix(0, nrow = N, ncol = 1)
for (control_sample in 1:N) {
  control_sample_i = read.table(paste(DIR,"/control_samples/",control_samples_list[control_sample], sep = ""), stringsAsFactors = F, head=F)
  l = log2(case_normalized_coverage[,1]/control_sample_i[,1])
  vars[control_sample,1] = var(l[which(!is.na(l) & l!=Inf & l!=-Inf)])
}


control0 = read.table(paste(DIR,"/control_samples/",control_samples_list[1], sep = ""), stringsAsFactors = F, head=F)
control0[which(control0[,1]==0.001),1] = 0
for (i in 1:N) {
  ord = order(vars[,1])
  control = control0
  for (j in ord[2:5]) {
    control_sub = read.table(paste(DIR,"/control_samples/",control_samples_list[ord[j]], sep = ""), stringsAsFactors = F, head=F)
    control_sub[which(control_sub[,1]==0.001),1] = 0
    control = control + control_sub
  }
  normalized_coverage = control - control0
  qlow = as.numeric(quantile(normalized_coverage[which(normalized_coverage[,1]!=0),1], 0.05))
  qup = as.numeric(quantile(normalized_coverage[which(normalized_coverage[,1]!=0),1], 0.9998))
  wq = which(normalized_coverage[,1]<qlow | normalized_coverage[,1]>qup)
  
  normalized_coverage = as.numeric(normalized_coverage[,1])/4
  normalized_coverage[wq] = -1
  ## write result to files
  ind_i_st = 1
  for (i in 1:24) {
    chr_i_bin_count = chr_sizes[i,2]%/%binsize + 1
    chr_i_coverage = matrix(0, nrow = chr_i_bin_count, ncol = 3)
    
    ind_i_en = ind_i_st + chr_i_bin_count - 1
    chr_i_coverage[,1] = 0:(chr_i_bin_count - 1)*binsize
    chr_i_coverage[,2] = 1:chr_i_bin_count*binsize
    chr_i_coverage[,3] = normalized_coverage[ind_i_st:ind_i_en]
    ind_i_st = ind_i_st + chr_i_bin_count
    write.table(chr_i_coverage, paste(outdir_path, "/chr", as.character(i),"_normalized_coverage_control", sep = ""), row.names =F, col.names = F, append = F)
  }
}
