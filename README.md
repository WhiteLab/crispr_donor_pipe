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
#### input_ids example:
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


# create_fa_from_ref
## Purpose:
Creates fasta file using list of gene symbols or refseq IDs from reference file.  Useful for creating input batch fasta file for site to work with crispr_guide_select.  Outputs to stdout and also prints an error if some ids not found.
## Usage:
### Mac/Linux:
```
python ./create_fa_from_ref.py (<id_list> <ref>)  [options] > output.fa
```
### Windows
```
path_to_python.exe create_fa_from_ref.py (<id_list> <ref>)  [options] > output.fa
```

Arguments:  
    &lt;id_list&gt;    new line-seperated list of gene symbols or refseq IDs    
    &lt;ref&gt;      tab-separated reference file with header - ids, symbols, sequence  

Options:  
    -h --help  

#### id_list example:
NM_001135564  
ZSCAN16  
NM_001242333  
CLOCK  
#### ref example:
refseq_ID	Gene_symbol	Seq  
NM_017734	PALMD	TTTAGAGTAAAAATTGAGGATTTTTTTTATAACTCTCTTTTTTTTTCTTTAATTTGCAGCTTTAAGGATGAGAATGGCAAAGCTGGGAAAAAAGGTGATCTAAGAGTTGTACCACCTATATAAACATCCTTTGAAGAAGAAACTAAGAAGCATTTGCAAATTTCTCTTCTGGATATTTTGTTTATTTTTTCTGAAGTCCAAAA  


# crispr_guide_select
## Purpose:
Parses pasted output from http://crispr.mit.edu/ and selects prioritized guide to fit our needs
## Usage:
### Mac/Linux:
```
python ./crispr_guide_select.py (<flank> <pasted>) [options] > output.fa
```
### Windows
```
path_to_python.exe crispr_guide_select.py (<flank> <pasted>) [options] > output.fa
```

Arguments:  
    &lt;flank&gt;    100 bp flanked stop codon reference file    
    &lt;pasted&gt;      pasted output from browser  

Options:  
    -h --help  

#### flank example:
refseq_ID	Gene_symbol	Seq  
NM_017734	PALMD	TTTAGAGTAAAAATTGAGGATTTTTTTTATAACTCTCTTTTTTTTTCTTTAATTTGCAGCTTTAAGGATGAGAATGGCAAAGCTGGGAAAAAAGGTGATCTAAGAGTTGTACCACCTATATAAACATCCTTTGAAGAAGAAACTAAGAAGCATTTGCAAATTTCTCTTCTGGATATTTTGTTTATTTTTTCTGAAGTCCAAAA  

#### pasted example:
CLOCK  
Guide #1	79	AGGAAGCACGTGTGCTACTGTGG  
Guide #2	77	GTGCTACTGTGGTTGAACCTTGG  
Guide #3	75	GGTTGAACCTTGGAAGGGTCGGG  
Guide #4	72	AGCAGCAACTCAGCCGGCACAGG  
Guide #5	72	TACTGTGGTTGAACCTTGGAAGG  
Guide #6	69	TCCTCCCTTGATGTCAAGAGAGG  
Guide #7	68	ACTGTGGTTGAACCTTGGAAGGG  
ZSCAN16  
Guide #1	86	AGAGCTGAACTTACTCGGAAAGG  
Guide #2	68	CAAAAGTTCTTTGGACACTCAGG  
Guide #3	64	TAATAAGAGCTGAACTTACTCGG  
Guide #4	63	AACCTTGAGTTTCCTCAATGTGG  
Guide #5	60	AGATTCTTTGATAACTAGTTAGG  

# Outputs to stdout example:
Gene	Category	Distance	Sequence	Sense	Name  
ZSCAN16	UTR	15	CAAAAGTTCTTTGGACACTCAGG	+	Guide #2  
CLOCK	STOP	0	TACTGTGGTTGAACCTTGGAAGG	-	Guide #5  