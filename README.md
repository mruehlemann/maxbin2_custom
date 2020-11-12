# maxbin2_custom script
Decluttered and enhanced scripts for [MaxBin2](https://sourceforge.net/projects/maxbin/)

This repository holds custom scripts for using MaxBin2 on previously identified aminoacid sequences using [prodigal](https://github.com/hyattpd/Prodigal). It skips the gene identification step and directly starts at the identification of marker genes using [hmmsearch (v3.3.1)](http://hmmer.org/).

```
MaxBin - a metagenomics binning software [Edited by: Malte Ruehlemann (m.ruehlemann\@ikmb.uni-kiel.de)].
Usage:
  run_MaxBin_prodigal.pl
    -contig (contig file)
    -prodigal (prodigal amino-acid output file)
    -out (output file)

   (Input reads and abundance information)
    [-abund (abundance file) -abund2 (abundfile) -abund3 (abundfile) -abund4 ... ]

   (You can also input lists consisting of reads and abundance files)
    [-abund_list (list of abundance files)]

   (Other parameters)
    [-min_contig_length (minimum contig length. Default 1000)]
    [-max_iteration (maximum Expectation-Maximization algorithm iteration number. Default 50)]
    [-thread (thread num; default 1)]
    [-prob_threshold (probability threshold for EM final classification. Default 0.9)]
    [-plotmarker]
    [-markerset (marker gene sets, 107 (default), 40 (Bacteria & Archaea),
      or 120 (recommended; GTDB Bacteria Marker Genes), and 122 (GTDB Archaea Marker Genes)  See README for more information.)]

  (for debug purpose)
    [-version] [-v] (print version number)
    [-verbose]
    [-preserve_intermediate]

  Please specify -abund information.
  You can input multiple abundance files at the same time.
  Please read README file for more details.\n
  ```
  
 The files gtdb_bac120.hmm and gtdb_ar122.hmm have been generated using the script in the newscripts folder. They hold HMM profiles for the marker genes defined by the [Genome Taxonomy Database](https://gtdb.ecogenomic.org/) and used in the [GTDB toolkit](https://github.com/Ecogenomics/GTDBTk). 
 
 To use these profiles, specify the option ```-markerset 120``` (for Bacteria) or ```-markerset 122``` (for Archaea).
 
 All files should simply be placed in the unpacked MaxBin-2.2.7 folder (this is the version it was tested with).
 The headers in the prodigal amino-acid ouptut file should have the prodigal standard format: ```>[CONTIG_NAME]\_[###]```
 
