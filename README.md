# crispr_donor_pipe
## Purpose:
To generate primers using refseq IDs, a table with sequence flanking stop codons of trancripts and primer3 software.  Requires that primer3 software already be installd/compiled
## Usage
### Mac/Linux
```
python crispr_donor_pipe.py (<config> <list> <LR_len> <RF_len> <min> <max>) [options]
```
### Windows
```
path_to_python.exe crispr_donor_pipe.py (<config> <list> <LR_len> <RF_len> <min> <max>) [options]
```

Arguments:  
    &lt;config&gt;    json formatted config file with refs and tool locations  
    &lt;list&gt;      list of refSeq IDs, one per line  
    &lt;LR_len&gt;    length of left reverse sequence primer to try  
    &lt;RF_len&gt;    length of right forward sequence primer to try  
    &lt;min&gt;       min length to use when trying to find an alternative primer  
    &lt;max&gt;       max length to use when trying to find an alternative primer  

Options:  
    -h --help  

## Inputs
#### config_example:
{  
  "primer3": "/Users/Miguel/Documents/Miguel/2016-Nov-2_moran_req/primer3-2.3.7/src/primer3_core", # location of compiled primer3 software  
  "master": "/Users/Miguel/Documents/Miguel/2016-Nov-2_moran_req/crispr_donor_pipe/REFS/hg19_stop_region_master.txt.gz", # location of rerference stop flank sequences(see below)   
  "Lsettings": "/Users/Miguel/Documents/Miguel/2016-Nov-2_moran_req/crispr_donor_pipe/REFS/L_primer3_settings.txt", # primer3-specific settings file to use for left flank  
  "Rsettings": "/Users/Miguel/Documents/Miguel/2016-Nov-2_moran_req/crispr_donor_pipe/REFS/R_primer3_settings.txt", # primer3-specific settings file to use for right flank   
  "lf_gibson": "tgaattcgagctcggtaccc", #gibson sequences to prepend to primer results  
  "lr_gibson": "cttgcagtagttggaatccc",  
  "rf_gibson": "ttgacgagttcttctgagtc",  
  "rr_gibson": "agcttgcatgcctgcaggtc",  
  "GC_trigger": "0.5" # GC content midpoint to use when detmining whether to try smaller or larger sequences when wasn't found with defaults  
}  
#### Note on primer3 settings file:
The first time you set up, you must add or edit the path to the primer3 melting temp config.  For example, add the following adjusting the path:
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/Miguel/Documents/Miguel/2016-Nov-2_moran_req/primer3-2.3.7/src/primer3_config/
#### input_ids_example:
NM_001135564
NM_001242333
NM_024303
NM_006492

#### master table example:
refseq_ID       Gene_symbol     L_seq_region    R_seq_region  
NM_017734       PALMD   CAGAGCACATTCTCTCAGTTCATCCTAAATTCTAGGGAATATTCTTGGCCCTGCCACTAGAACTGATCATCCCTGTTAAAGTTTCAAGGGAAATAATGCTCTTTAGCTTCATTGAAATTAAAATATTTCATATAAGTGGAGTGACAATTAAAATAATCATGGTCTTGTTACAGAATTAGATCCCATCCCTTCTCATATTTGTGCAAAAAGGTTTCTTAGGTATGGGGTCATTTCCATTGTCAAAAGGTAACCAATATATGGTATACAAATGAGAAAACATGTCTTCCATTTTACTTTCCATGAGACAAAATTGTGCAATCAGTTACATAGCTTGTTCTTGGGCCATATTACATGTGAATGTAATATTAGATCATTAACATAATCTTTAGATCAGTGCCTTGAAAGCATCTATCAACTTTTGGAACATGACTAAGGCAGCTCTACACTGAAATCCTAAATAATAAGAATAGAATACCAGAATACTCATCCAGAATACCAGAAAAGTGTTATACATTCTAAATCTACATTTTAAAAAATGCTTTCTACATTCAATATATTATTAGCAGAGGAATATATCCAAAAAGATATATTTGCCATTTTAAAACATACTTATATCTTGACAAAGACATTTTATAATCACTCCAGCTTTTAATAACTTGTTCCTGTTGTTTCATTCGTAAAGAAATTTGTTTTTGTTATTGTATATCCAAAGGCCAGAGCTCCTTTGGGAAAGGAAGTTGAATGAGCTGAAAGCAAACACGATCATGAAATAGAGTTTAAGTGTAAACACCATTTTATTTTGGCAGTACATAACTCAGAAATGTATGTAAAACAGTCATCAACCTATAGGAATCTGTATGTATGGCTAATAGGGAAACTGTGCATTAATCCTGTATTTGTTCCACTTGTACTGTTTGCTCTACTAATTTTTCTAATCTTTATAGAAATAGGTTCCTTCTAGTTCTGAACCTGCTTCTGCCCTAAGTTAACCATAGCTTATTATTTCTACTATCTTTTTGTGCTACTAAATCAATGCTCACACTGTTGGATAAAACCATCTGAGGTATGTTCTAAACTTTCCTGCATTCTATGAATGGCCATTTAGAGTAAAAATTGAGGATTTTTTTTATAACTCTCTTTTTTTTTCTTTAATTTGCAGCTTTAAGGATGAGAATGGCAAAGCTGGGAAAAAAGGTGATC        GAGTTGTACCACCTATATAAACATCCTTTGAAGAAGAAACTAAGAAGCATTTGCAAATTTCTCTTCTGGATATTTTGTTTATTTTTTCTGAAGTCCAAAAAATTATCATTACAGTGTACCATATTAAGCCATGTGAATAAGTAGTAGTCATTATTTGTGAAAAATTCCCAAAAAGCTGGGGAAAACAAATGTGTAACTTTTCCAGTTACTTGACACGATTCAGTGGGGGAAAACCAGCATTTTTTATTCTATTGATACCAAAGCATTTCTAATAAGAGCTTGTTAAATTTAAGAATAAAGTTATTTAAAATATTCTGAGTATAGTATATTAACTGGCATTGTAATTTTGATGATACAAAGATTGAAAGATCATAGGAAAGCATTGCCCTTCATCACAGAAGTATTCAACTCTGACAAATAAATATGTCATCCTGAATTAATAATGCCTTAATAAAAGTACATCCTCCTGCTAACTATCACTTCAGTTGGTTTGGTAGAGACTTCAAACCCTTGCACATATCCAAAAGGTGCATCCCGTAATGATCCTCAGATTTGGGGACTGAGTTTGGATAGTCCGGCCAGTGCACAGATGATAACCGGCACATTCATTCTTCTTAAATTATCATCCTTTCTTTCACTCAGTTCAATGTGAACCCTGTGGTGCAATGACAGACTATGACCTCAGGGTGTTGGATAGAGGATGTCCCATACCCTCAAAGAGCACAGATCTGGTCTAATATTTCACTGATTCTTGTGAACAAGTTAGACAGGGGTGGCCCATCTTTTGGCTTCCCTGAGCCACATTGGAAGAATTGTCTGGGCCACACATAAAATACACGAACACTAACAGTAGCTGATGATCTAAAAAAAAACGCAAAATAATCTCATACTGTTTTAATAAAGTTTATGAATTTGTGTTGGGATGCATTCAAAGCCATCCTGGGCCACATGCGGCCCGCGGGCGGCAGGTTGGGCAAGCTTGAGTTCGAATCTTAATATTAAAACTAATAGGAATGTTTACAGCATATTACTTCTAAGTATTCTTGCCAAGATATTCAATCCTAATTTGAAGCAGAGTACTCCCATGACAATAAATTCTTTCCTACCCTGCCTTATTTTACGTCCCTCCTTCAAACCCCCTCGAGATCCACTCCTACACACCTTCCCCTGAATCATGCCAGTATATATTATCCAGCTCAT  

## Outputs
Note - all output files prepended today's date, time in 24 hour format as an integer, and lengths of the desired left and right primers
#### Temp folder
2016-11-10_1308_20_20_TEMP # contains all primer3 input SEQUENCE files and output results
#### Results file
2016-11-10_1308_20_20_results.xls # contains parsed results with primers found and their melting temps
#### Warning log
2016-11-10_1308_20_20_warnings.txt # outputs whether a refseq was not found in the reference and results of attempts to try different lengths of primers when one could not be found for the default
