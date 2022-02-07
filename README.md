
This script produces peptide and MHC moleculs binded by a predefind set of amino acids.
The outputs are two separate files. 

`--sample_num` number of samples (pairs of peptide and MHC molecule) to produce

`--decoys` number of decoys per sample

`--peptide_min` minimum length of the peptides

`--peptide_max` maximum length of the peptides

`--mhc_min` minimum length of the MHC molucles

`--mhc_max` maximum length of the MHC molucles

```
python .\data_creator.py --sample_num 100 --decoys 99 --peptide_min 7 --peptide_max 15 --mhc_min
 15 --mhc_max 35
```
