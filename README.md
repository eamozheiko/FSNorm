# FSNorm

This script help to get normalized coverage from DNAse Hi-C libraries. It results in 24 files named as chr${chr}_normalized_coverage_case_vs_control. The format of this file is:
start_position end_position normalized_case_value normalized_control_value normalized_case_value/normalized_control_value

For every case script search for most suitable four control samples from control semples set. Control semples are in the folder /control_sample.
Negative value means undefined value in the bin (too low value in the control).

## running FSNorm

Enter paths to input files such as valid_pairs_file_path and outdir_folder_path in FSNorm.sh file.

to run program from command line:
```
./FSNorm.sh
```
## valid_pairs_file format:

chr1 pos1 chr1 pos2

examples:
```
1 10000 X 20000
```

## requirements 

R, Shell
