#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"



#####################################################################################################
### PATHS TO INPUT FILES
#####################################################################################################
## default
capture_track_path=${DIR}/MedExome_hg19_capture_targets.bed
FS_track_path=${DIR}/fs_track_100kb.txt
chr_sizes_path=${DIR}/chr_sizes_hg19
control_coverage_track_path=${DIR}/control_samples
binsize=10000

## edit here
valid_pairs_file_path=
outdir_folder_path=





#####################################################################################################
## preparations
#####################################################################################################

mkdir ${outdir_folder_path}
cd ${outdir_folder_path}

rm coverage
rm capture_coef_pre
echo ${binsize} > binsize

cat ${capture_track_path} | awk '{
print $1 " " $2 " " $3
}' > capture_exome
time sed -i -r 's/X/23/g' capture_exome
time sed -i -r 's/Y/24/g' capture_exome

cat capture_exome | awk -F'chr' '{
print $2
}' > capture_exome1



time sort -k1,1n -k2,2n capture_exome1 > capture_exome2




#####################################################################################################
## bin valid_pairs and evaluate capture coefficient
#####################################################################################################

## read valid pairs file
cat ${valid_pairs_file_path} | awk '{
print $0
}' > allValidPairs_all
time sed -i -r 's/X/23/g' allValidPairs_all
time sed -i -r 's/Y/24/g' allValidPairs_all


for chr in {1..24}
do
echo ${chr}
echo ${chr} > chr


# capture
cat capture_exome2 | awk 'BEGIN{getline chr < "chr"}{
if($1==chr){
  print $0
}
}' > capture_exome3
echo "25 0 0" >> capture_exome3
cat capture_exome3 | awk '{
print $1
}' > capture_exome_chr

cat capture_exome3 | awk '{
print $2
}' > capture_exome_st

cat capture_exome3 | awk '{
print $3
}' > capture_exome_en



## generate bins
cat ${chr_sizes_path} | awk 'BEGIN{getline chr < "chr"; getline binsize < "binsize";}{
  if(NR==chr){
  max = $2/binsize
  max = max - max%1
  max = max + 1
  for (i = 1; i <= max; ++i){
    bin_st=(i-1)*binsize
    bin_en=bin_st+binsize
    print bin_st " " bin_en
    }
  }
}' > chr${chr}_bins

cat chr${chr}_bins | awk '{
  print $1
}' > bin_st

cat chr${chr}_bins | awk '{
  print $2
}' > bin_en

rm chr${chr}_bins

## get coverage track
cat allValidPairs_all | awk 'BEGIN{
getline chr_j < "chr"
}{
chr1=$1
pos1=$2
chr2=$3
pos2=$4
if(chr1==chr_j){
  print pos1
}
if(chr2==chr_j){
  print pos2
}
}' > allValidPairs_cnv_coord_pre

time sort -k1,1n allValidPairs_cnv_coord_pre > allValidPairs_cnv_coord
rm allValidPairs_cnv_coord_pre

cat allValidPairs_cnv_coord | awk 'BEGIN{
getline bin_st < "bin_st"
getline bin_en < "bin_en"
getline chr < "chr"
val = 0
}{
coord = $1
while (coord>bin_en){
  print chr " " bin_st " " bin_en " " val
  getline bin_st < "bin_st"
  getline bin_en < "bin_en"
  val = 0
}
val = val + 1

}' >> coverage


## capture coefficient
cat allValidPairs_cnv_coord | awk 'BEGIN{ 
getline chr_exome < "capture_exome_chr"
getline chr < "chr"
getline st < "capture_exome_st"
getline en < "capture_exome_en"
n=0
k=0
}{
if(chr==chr_exome){
n=n+1
while($1>en && chr_exome==chr){
  getline chr_exome < "capture_exome_chr"
  getline st < "capture_exome_st"
  getline en < "capture_exome_en"
}
if($1<st) {
  d=1
} else {
  k=k+1
}
}

}END{ print k " " n }' >> capture_coef_pre
done

cat capture_coef_pre | awk 'BEGIN{k=0;n=0;}{
k = k + $1
n = n + $2
}END{
k = k/1.5
n = n/98.5
print k/n
}' > capture_coef

# rm byproducts
rm chr
rm allValidPairs_all
rm allValidPairs_cnv_coord_pre
rm allValidPairs_cnv_coord
rm bin_st
rm bin_en
rm capture_exome_st
rm capture_exome_en
rm capture_coef_pre
rm capture_exome_chr
rm capture_exome
rm capture_exome1
rm capture_exome2
rm capture_exome3


#####################################################################################################
## run fs_normalization Rscript
#####################################################################################################
capture_coefficient_path=${outdir_folder_path}/capture_coef
coverage_track_path=${outdir_folder_path}/coverage
Rscript ${DIR}/FSNorm.R ${coverage_track_path} ${capture_track_path} ${FS_track_path} ${chr_sizes_path} ${binsize} ${capture_coefficient_path} ${outdir_folder_path}

#####################################################################################################
## search for the most suitable control sample
#####################################################################################################
ls ${control_coverage_track_path} > control_samples_list
Rscript ${DIR}/ControlSearch.R ${outdir_folder_path} ${DIR} ${binsize}

rm control_samples_list

#####################################################################################################
## dividing a sample into a control 
#####################################################################################################
for chr in {1..24}
do
sed -i -r 's/,/./g' chr${chr}_normalized_coverage_control
ed -i -r 's/,/./g' chr${chr}_normalized_coverage 
cat chr${chr}_normalized_coverage_control | awk '{
print $3
}' > control_values
cat chr${chr}_normalized_coverage | awk '{
getline control_value < "control_values" 
val = $3/control_value
val = val*1000
val = val - val%1
val = val/1000
print $1+0. " " $2+0. " " $3 " " control_value " " val
}' > chr${chr}_normalized_coverage_case_vs_control
rm control_values
rm chr${chr}_normalized_coverage
rm chr${chr}_normalized_coverage_control
done

exit





