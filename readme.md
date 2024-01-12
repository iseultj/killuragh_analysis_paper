# Analysis pipeline for "Ancient genomes from Bronze Age remains reveal deep diversity and recent adaptive episodes for human oral pathobionts"

This file outlines all steps and commands needed to reproduce the analysis in the paper "Ancient genomes from Bronze Age remains reveal deep diversity and recent adaptive episodes for human oral pathobionts".

Where appropriate, scripts are available in ./scripts directory; otherwise one-line commands are outlined below. 





## Part 0: Data Download

Sequencing data generated as part of this project is archived on the ENA (Accession: PRJEB64128).

Intermediate data files are archived on Zenodo (doi: https://doi.org/10.5281/zenodo.10024810)

This code repository is organised assuming that the ./data directory is structured identically to the archive on Zenodo. 

The initial sequence data processing scripts are independent of this structure, and assume that you've put the raw FASTQ files in a directory called raw_data.

**NB: Make sure you put paired end FASTQs in a different directory to the single end FASTQs**



## Part 1: Raw sequence data processing and reference alignment

This step is to be carried out on the submitted fastq files on the ENA; Paired and single end files must be treated separately.

### a. FASTQs aligned to human reference genome

LCLab alignment script (seqPro_v7.py. [https://github.com/LCLabTCD/aDNA_data_processing])
Parameters used in paper:

`python2.7 /Programs/seqPro_v7.py -RGCN Macrogen_Korea -RGPM HiSeq2500 -end SE -un Y # for HiSeq data`

`python2.7 /Programs/seqPro_v7.py -end PE -un Y # for Novaseq`

`python2.7 /Programs/seqPro_v7.py -end SE -un Y # for MiSeq`

The same command is used for reference alignment to bacterial reference genomes, with the additional option `-R /path/to/reference/`
​	See **Table S7** for accessions to use.

Additionally, for bacterial alignments, reads were filtered for minimum length 34, mapping quality 25 and were softclipped for 2bp at the start and end of reads (UDG treated; 5bp for non-treated data).

See [https://github.com/LCLabTCD/aDNA_data_processing] for alignment and post-alignment post-processing and quality control scripts.

### b. Unaligned BAM processing & Metataxonomic profiling

Use `./scripts/processing/metascreen.py` in unaligned directory.

Bracken reports are combined to OTU tables using kraken-biom

`ls *report  |  xargs kraken-biom --fmt tsv -o killuragh_teeth.dental_data_with_sources.tsv --max S`  

Change --min and --max for different taxonomic levels

## Part 2: Metataxonomic analysis

This analysis is carried out on the OTU table generated in Step 1b.

NOTE: for sourcetracker analysis, you must include source data as well. The OTU table archived in ./data/ already includes this.

### a. Sourcetracker analysis:

For output files, see: `./data/results/metagenomics/sourcetracker_species_rarefy_2000_17-06-22/mixing_proportions.txt`

To run yourself:

`sourcetracker2 gibbs -i ./data/results/metagenomics/OTU_tables/killuragh_teeth.dental_data_with_sources.25-04-22.species.tsv -m ./data/results/metagenomics/metadata/Map_17-06-22.txt  -o sourcetracker_species_rarefy_2000_${date} --jobs $N --source_rarefaction_depth 2000 --sink_rarefaction_depth 2000` 
	

Plotting output: 

see script in `./scripts/metagenomic_analysis/sourcetracker_results_plot.R`

### b. Relative abundance plotting & calculations from OTU tables.

Input: species/genus/order level OTU tables & metadata from` ./data/results/metagenomics/`

Paths are hard-coded into the script- if you are running from a different directory to the one you have the data directory in, you'll need to edit them. 

## Part 3: Metagenomic Assembly

This analysis is carried out on trimmed reads downloaded from the ENA project: PRJEB63463.

Commands used for each step are outlined in `./scripts/metagenomic_assembly/`
YAML files for necessary conda environments are also in this directory.

1. `metagenomic_assembly_steps.pre_concoct.sh` : prepare data for input to CONCOCT for binning
2. Binning: used docker for concoct:



```
docker run -v /path/to/metagenomic_assembly/:/opt/Data -i -t binpro/concoct_latest bash
```

Within docker container:

```
mkdir -p Concoct
cd Concoct/
concoct --coverage_file ../Coverage.tsv --composition_file ../Assembly/killuragh.contigs_c10k.fa
```

3. `metagenomic_assembly_steps.post_concoct.txt` : split bins and assess using checkm and gtdbtk

## Part 4 : *S. mutans* analysis

This analysis is carried out on the assembly from Part 3 and assemblies downloaded from NCBI: see SI table for accessions. 

YAML files for conda environments used in these analyses are available here:

`./scripts/mutans_analysis/conda_envs`

These are:

-  PROKKA
- roary
- diamond
- gubbins
- snp-sites
- IQTree
- BEAST
- gargammel

Metadata for sequences downloaded from NCBI were obtained from the associated GBFF files using the script `./scripts/mutans_analysis/metadata_ext_gbff.py`

### a. Pangenome analysis

- Annotation using PROKKA:

  Refseq genes from the UA159 assembly were used as a training set. 

`prokka --outdir ${id}_prokka_${date} --kingdom Bacteria --proteins GCF_000007465.2_ASM746v2_cds_from_genomic.fna ${id}.fa`

- Pangenome inference using the pangenome pipeline roary

Command used:

`roary -f roary_${date} -p ${threads} -e --mafft -r -v  ./path/to/gffs/*gff`

- Additional annotation for functional analysis using diamond blastx:

  Database used: NCBI nr, downloaded 13th January 2022. 

  `diamond blastx -d /path/to/diamond_ncbi-nr_13-01-22.dmnd --threads ${threads} --log -q pan_genome_reference.fa --fast -o S_mutans_pangenome.tsv`

Subsequent functional analysis scripts (analysed in tandem with phylogenetic data) are described below



### b. Multiple sequence alignment and phylogenetic analysis.

- Multiple sequence alignment: SKA (requires gubbins conda environment)

  Generates multiple sequence alignment to reference.

  `generate_ska_alignment.py --reference S_mutans_UA159.ASM746v2_genomic.fna  --fasta list_of_fastas.txt --out S_mutans_aln.ska.${date}.fa`

- Gubbins: mask recombinant regions

Input: multiple sequence alignment from SKA

Output: recombination predictions as GFF file - this is then used to create masked alignment 



` run_gubbins.py --prefix S_mutans_aln.ska.${date}.gubbins.100_bootstraps.rm_identical S_mutans_aln.ska.${date}.fa  -c ${cores} --remove-identical-sequences  --bootstrap 100`

Mask alignment:

`mask_gubbins_aln.py --aln S_mutans_aln.ska.${date}.gubbins.100_bootstraps.rm_identical.filtered_polymorphic_sites.fasta --gff S_mutans_aln.ska.${date}.gubbins.100_
bootstraps.rm_identical.recombination_predictions.gff --out S_mutans_aln.ska.${date}.gubbins.100_bootstraps.rm_identical.masked.fa`

- Filter polymorphic sites from alignment for input to phylogenetic analysis
	
- Step 1: Convert multiple sequence alignment to VCF using snp-sites (use conda environment snp-sites)

	`snp-sites -v -o S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.snp-sites.vcf S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.fa`
	
	Monomorphic sites (required for beast) also obtained:
	
	`snp-sites -v -m -o S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.monomorphic.vcf S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.fa`
	
	- Step 2: Remove multiallelic sites and ensure missing data is encoded correctly
	
	Script: `./scripts/mutans_analysis/fix_gubbins_snp-site_vcf.py`
	
	- Step 3: Filter for maximum 5% missingness using plink1.9
	
	`plink1.9 --vcf S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.snp-sites.recode_multiallelic.vcf --geno 0.05 --recode vcf --out S_mutans_aln.ska.25pc_missingness.mask
	_alignment_${date}.snp-sites.recode_multiallelic.geno05`
	
	**NOTE: Make sure that if any IDs that are delimited with "_", they are converted correctly.** 
	
	- Step 4: convert VCF back to multiple sequence alignment:
	
	First, convert VCF to tab
	
	`vcf-to-tab < ${input} > merged_tabbed.vcf`
	
	Then `./scripts/mutans_analysis/vcf2fasta.py ${n_samples}` 

Concatenate resulting fastas to create multi-fasta file to be used as input for phylogenetic analysis. 

- IQtree:

Use conda environment iqtree.

Use dates file extracted from gbff metadata files (see above).

Command used:

```
seed=$(echo $RANDOM)
iqtree -s S_mutans_aln.ska.25pc_missingness.mask_alignment_${date}.snp-sites.recode_multiallelic.geno05.fna -pre S_mutans.beast_sites.1000boots.seed${seed} -m MFP --date mutans_dates.txt --date-ci 100 --seed ${seed} -B 1000 -mem ${mem} -nt ${threads}
```



- BEAST analysis:
  - xml file included here: `./scripts/mutans_analysis/beast/S_mutans_aln.snp-sites.rm_multi.geno05.bmod.bcoalsky.500mchain.1millburnin.xml`
  - Output logs , tree files & final maximum clade credibility tree here: `./data/results/beast/`

### c. Functional and phylogenetic data

Plot BEAST tree in conjunction with mutacin data:

`./scripts/mutans_analysis/beast/tree_and_mutacin_analysis.R`

IQtree plotting:

`./scripts/mutans_analysis/iqtree_beast_comparisons.R`

### d. Mutacin Simulations

- data: input for simulations in `./data/intermediate/mutans/mutacin_simulation/`
- scripts: run simulation script `./scripts/mutans_analysis/gargammel_1000_sims_for_all.sh` in mutacin simulation directory
- Align output fastqs to pangenome reference as in reference genome alignment section (**Part 1**)

## Part 5: *T. forsythia* analysis

### a. Data preparation: Alignment

Raw FASTQ files aligned to reference as in **Part 1**

Modern FASTA files were converted to pseudo-reads using `./scripts/forsythia_analysis/chop_fasta_file_gzip.py` and aligned to the same reference as the ancient FASTQ files, and were subject to the same post-alignment filtering as the ancient files, but were not softclipped (see **Part 1**).

### b. SNP Ascertainment and Calling

SNPs were ascertained in modern and ancient alignments with >= 2X average genomic coverage.

`bash ./scripts/forsythia_analysis/SNP_ascertainment.forsythia.sh list_of_bams.txt`

These sites were then called in the ancient and modern BAM files with ` ./scripts/forsythia_analysis/HapSNP_v10.py`, using  `-consensus T -mincov 2` to call consensus SNPs for these samples, with at least 2 reads covering called sites. 



Heteroplasmy at the ascertained sites was assessed using `./scripts/forsythia_analysis/bacterial_heteroplasmy.py` 



Distribution of singletons using different calling filters can be assessed using data from`./data/intermediate/forsythia/singleton_analysis` and the script 

`./scripts/forsythia_analysis/singleton_check.R`



### c. Phylogenetic Analysis



#### IQtree:

Maximum per-individual SNP missingness of 40%; this was relaxed to include the African sequence TAF008, the neanderthal-derived sequence GOY005 and KGH2-F.

Maximum per-site missingness 5%

SNPs converted to multi-sequence fasta using

 `./scripts/mutans_analysis/vcf2fasta.py ${n_samples}`

Initial phylogenetic analysis with dog outgroup:

```
for i in {1..4}; do seed=$(echo $RANDOM); iqtree -s T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.
mind40.keep_special.fna -pre run${i}/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.seed_${seed}.B1000.iqtree --seed ${seed} -B 1000 -m MFP --date t_forsythia_dates.t
xt --date-ci 1000 -nt 5 -o OH2617-COT023; done

```



IQTree with Neanderthal outgroup:

```
for i in {1..4}; do seed=$(echo $RANDOM); iqtree -s T_
forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.geno05.fna -seed ${seed} -pre run${i}/T_forsythia.ind40.keep_spe
cial.geno05.seed_${seed}.B1000.iqtree -B 1000 --date ../t_forsythia_dates.txt --date-ci 1000 -m MFP -nt 3 -o GOY005; done

```



#### BEAST:

Temporal signal in medieval data assessed with IQTree without date information



```
for file in $(ls *cod*/*fna); do out=$(echo $file | cut -d/ -f1); seed=$(echo $RANDOM); ec
ho "iqtree -s $file -pre tempest_tests_partitions/${out}.seed_${seed} -m MFP -o OH2617-COT023 --seed $seed -B 1000"; done > iqtree_partitions_no_dates.sh
```



Invariant sites for BEAST ascertainment correction:

```
plink --bfile T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.homozyg --keep keep_medieval_onwards.txt --max-maf 0.0
01 --recode rlist --out medieval.invariant_sites.rlist --allow-extra-chr
## For partitioned BEAST, extract codon positions & non-coding sites for ascertainment correction for each partition

```



BEAST XML file:



`./scripts/forsythia_analysis/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.medieval_on_no_dog.geno05.rm_invariant.partitioned.xml`

#### Phylogenetic Data Plotting 

IQTree plots:

`forsythia_analysis/phylogenetic_analysis/plot_iqtree_output.R`



BEAST and Bayesian skyline plots:



`forsythia_analysis/phylogenetic_analysis/medieval_beast_data_plot.R`

### d. Functional analysis

Coverage across virulence factors:

```
for bam in $(ls /path/to/*bam); do id=$(basename $bam | cut -d- -f1); bedtools coverage -a ./data/results/forsythia/functional_analysis/92A2_virulence_factors.bed -b $bam > ${id}.92A2_virulence_factors.txt; done

```

Similarly for KLIKK Protease BED in same directory.



Coverage data concatenated to `virulence_factors_long.txt`

Script for processing & plotting this output:

`./scripts/forsythia_analysis/functional_analysis/virulence_analysis.R`



## Part 6: *T. denticola* analysis

Ancient samples were aligned to reference genomes and SNPs were ascertained as for *T. forsythia* (**Part 1, Part 5 b**).



Heteroplasmy levels for *T. denticola* alignments are in `./data/results/denticola/heteroplasmy`

Heteroplasmy aggregation accross genomes:

`./scripts/denticola_analysis/aggregate_heterozygosity_data.R`

Calculate rolling heterozygosity accross genome:

`./scripts/denticola_analysis/T_denticola_htz_estimates.R`

