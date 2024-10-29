
# MNase-seq analysis pipeline

## description

An analysis pipeline for paired-end MNase-seq data with the following major steps:

- alignment with [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) (.bam)
- indexing and sorting alignment files with [`samtools`](http://www.htslib.org/) (.bam and .bai)
- summary of fragment sizes with [`deeptools PEFragmentSize`](https://deeptools.readthedocs.io/en/develop/content/tools/bamPEFragmentSize.html)
- counting reads aligning to experimental (*S. cerevisiae*) and spike-in (*S.pombe*) genomes with [`samtools view`](http://www.htslib.org/doc/samtools-view.html)
- calculation of per-library spike-in normalization scaling factors using a custom python script
- generation of coverage tracks (.bw) scaled by spike-in normalization using [`deeptools bamCoverage`](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)
- generation of matrices for calculating correlation between datasets using [`multiBigwigSummary`](https://deeptools.readthedocs.io/en/develop/content/tools/multiBigwigSummary.html)
- calculation and plotting of correlation between datasets using [`plotCorrelation`](https://deeptools.readthedocs.io/en/stable/content/tools/plotCorrelation.html)
- generation of matrices for plotting scaled coverage data using deeptools (.gz) or directly in python (.tab) using [`deeptools computeMatrix`](https://deeptools.readthedocs.io/en/develop/content/tools/computeMatrix.html#reference-point)
- data visualization ([heatmaps](https://deeptools.readthedocs.io/en/develop/content/tools/plotHeatmap.html) and [metagenes](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html)) using `deeptools` plotting functions and custom python scripts

## requirements

- Designed to be run on O2, the high-performance computing cluster at HMS
- Uses slurm job scheduler to batch submit jobs
- Paired-end FASTQ files from MNase-seq libaries
	- FASTQ files should be demultiplexed and can (should) be compressed (.gz)
	- FASTQ filenames are used by the scripts throughout this pipeline, and should easily identify the sample. For example, my filenames usually take the format: strain_treatment_replicate (e.g. 425_93_D_rep1). After demultiplexing, there is usually trailing info added to the filename (e.g. _S3_R1_001.fastq.gz). The S indicates the index number on your sample sheet, R1 and R2 indicate the paired reads from the sequencer, and 001 is a trailing number added for reasons beyond my comprehension.
- Bowtie2 index files (.bt2) for genome alignment, included for *S. cerevisiae* and *S. pombe* in the `genomes/bowtie2_index/` directory in this repository.
- BED files to tell deeptools which portions of the genome you would like to plot. I have included the standard non-overlapping ORF BED file from James Chuang in the `genomes/annotations/` directory of this repository.


## instructions

**1. Clone this repository into a clean directory on O2.**

```bash
# clone the repostiory
git clone https://github.com/jamwarner/mnase_seq.git

#navigate to the newly created directory
cd mnase_seq
```

Run all commands from the base `mnase_seq/` directory.


**2. Download or transfers your .fastq.gz files into the `fastq/` directory.**


**3. Align your libraries to the experimental and spike-in genomes.**

We will submit two alignment jobs for each library: one to the experimental genome and the othe to the spike-in genome. Submitting the jobs separately allows all of the alignments to run in parallel.

```bash
# use a for loop to submit alignments for each set of paired reads separately
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/batch_aligner.sh $name; done
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/spike_batch_aligner.sh $name; done
```

This will generate two sets (experimental and spike-in) of three files for each library:
- in `bam/`
	- filename_unsorted.bam
	- filename_sorted.bam
	- filename_sorted.bam.bai
- in `bam/spike-in/`
	- filename_spikein_unsorted.bam
	- filename_spikein_sorted.bam
	- filename_spikein_sorted.bam.bai

Also generated is a summary of each alignment in the `logs/` directory:
- filename_bowtie2.txt

We will use the 'sorted.bam' and 'sorted.bam.bai' files in subsequent steps. I don't think that the 'unsorted.bam' files need to be saved, but I have not made a habit of deleting them.


**4. Determine the distribution of insert sizes in your MNase samples.**

Since the reads are paired, we can determine the size of each fragment that was sequenced from its two ends.

```bash
# use deeptools to generate summary statistics for each sorted.bam file
sbatch scripts/PEFragmentSize.sh
```

This script will look at all of the 'sorted.bam' files and will generate two new files in the `fragment_sizes/` directory:
- a histogram showing the distribution of the insert sizes for each library
- a CSV file with this information in tabular format. This is usually easier to interpret.

These files will be generated for the experimental and spike_in alignments separately.


**5. Count aligned reads for experimental and spike-in genomes.**

This step uses [samtools view](http://www.htslib.org/doc/samtools-view.html) to count the number of reads in each 'sorted.bam' file. 
- `-c` makes `samtools` output only the count and not the reads themselves
- `-F` filters the BAM files to EXCLUDE reads that fit the [flag condition](https://broadinstitute.github.io/picard/explain-flags.html)
	- the flag here is `388` which = unmapped OR not primary alignment OR second in pair (which confirms read pairs are only counted once and not twice)

To run:

```bash
sbatch scripts/read_counter.sh
```

This scripts creates two files in the `logs/` directory: 'experimental_counts.log' and 'spikein_counts.log'. These files contain pairs of lines:
- the first line is the BAM filename with '_sorted.bam' stripped
- the second line is the number of paired reads that are mapped to the respective genome in each 'sorted.bam' file
The lines continue to alternate for each BAM file processed.


**6. Do spike-in normalization math.**

To be completed.


**7. Generate coverage tracks for each library scaled by spike-in normalization.**

```bash
# use a for loop to submit each coverage job in parallel
for name in bam/*_sorted.bam; do sbatch scripts bamCoverage_spikein.sh $name; done
```


**8. Calculate and plot correlation scores**

```bash
sbatch scripts/multiBigwigSummary.sh
```

This script calculates coverage scores in non-overlapping 75 bp bins genome-wide (skipping the mitochondrial genome).


```bash
sbatch scripts/plotCorrelation.sh
```

This script takes the above binned coverage scores, calculates the Pearson's and Spearman's correlation coeeficients between all libraries, and plots the results as two heatmaps in the `correlation/` directory. It also outputs the results in tabular (.tab) format in the same directory if you would like to plot it yourself.


**9. Generate matrices to plot data.**

Using the coverage (.bw) files, generate matrices with which to plot the data. The default BED file lines all non-overlapping protein-coding genes up according to their +1 nucleosome. Plotting data is generated for 500 bp upstream and 1500 bp downstream of this position.

```bash
# use a for loop to submit each replicate in parallel
for rep in rep1 rep2 rep3; do sbatch scripts/computeMatrix_reference.sh $rep; done
```

There should now be new .gz files in the `deeptools/` directory that will be used in the next step to plot the results. The uncompressed matrices are also saved in `deeptools/tab` in tab delimeted format so that you can import them into python and plot the data manually.

**10. Plot data.**

```bash
# use a for loop to submit each replicate in parallel
for rep in rep1 rep2 rep3; do sbatch scripts/makeplots.sh $rep; done
```

Plots are found in the `plots/` directory.  Two types of plots are generated: profile (metagenes) and heatmaps.