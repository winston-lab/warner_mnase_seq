
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

**2. Prepare virtual environments**

The code written for these analyses relies on two virtual environments `deeptools` and `spike_in` that you will have to set up using `virtualenv` and `pip`. The following commands should be run line-by-line.
(More on virtual environments on O2 [here](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1588662166/Personal+Python+Packages).)

```bash
# start an interactive session
srun --pty -p interactive -t 0-1:00 --mem=1G bash
```

To create the `deeptools` environment:
```bash
# load modules
module load gcc/9.2.0

module load python/3.10.11

# create a local virtual environment
virtualenv env/deeptools --system-site-packages

# activate the virtual environment
source env/deeptools/bin/activate

# install packages to the virtual environment
pip3 install deepTools

#deactivate the environment
source deactivate
```

To create the `spike_in` environment:
```bash
# create a local virtual environment
virtualenv env/spike_in --system-site-packages

# activate the virtual environment
source env/spike_in/bin/activate

# install packages to the virtual environment
pip3 install numpy

pip3 install pandas

pip3 install matplotlib

#deactivate the environment
source deactivate
```

There should now be two folders, one for each environment, in the `env/` directory.

**3. Download or transfer your .fastq.gz files into the `fastq/` directory.**


**4. Align your libraries to the experimental and spike-in genomes.**

We will submit two alignment jobs for each library: one to the experimental genome and the othe to the spike-in genome. Submitting the jobs separately allows all of the alignments to run in parallel.

```bash
# use a for loop to submit alignments for each set of paired reads separately
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/batch_aligner_new.sh $name; done
for name in fastq/*R1_001.fastq.gz; do sbatch scripts/spike_batch_aligner_new.sh $name; done
```

Alignment is done using [`bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) using the following options:
- `--sensitive`
- `--no-unal`
- `--no-mixed`
- `--no-discordant`

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

We will use the 'sorted.bam' and 'sorted.bam.bai' files in subsequent steps. The 'unsorted.bam' files can be safely deleted.
```bash
rm bam/*unsorted*
rm bam/spike_in/*unsorted*
```


**5. Determine the distribution of insert sizes in your MNase samples.**

Since the reads are paired, we can determine the size of each fragment that was sequenced from its two ends.

```bash
# use deeptools to generate summary statistics for each sorted.bam file
sbatch scripts/PEFragmentSize.sh
```

This script will look at all of the 'sorted.bam' files and will generate files in the `fragment_sizes/` directory:
- a histogram showing the distribution of the insert sizes for each library
- a CSV file with this information in tabular format. This is usually easier to interpret.

These files will be generated for the experimental and spike_in alignments separately.


**6. Count aligned reads for experimental and spike-in genomes.**

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


**7. Do spike-in normalization math.**

For spike-in normalization, we will use the *S. pombe* reads to calculate a normalization factor.  Essentially, the normalization factor is a per-library number so that, when scaled by this factor, each library would have the same number of reads aligning to the *S. pombe* genome.

```bash
# start an interactive session (if not already interactive)
srun --pty -p interactive -t 0-1:00 --mem=1G bash

# load modules
module load gcc/9.2.0

module load python/3.10.11

# activate the virtual environment
source env/spike_in/bin/activate

# run the python script
python scripts/mnase_spikein_norm.py

# deactivate the environment
deactivate
```

The output of all of this is a file called 'normalization_table.csv' in `logs/` that consists of two columns:
- column 1 is the library name
- column 2 is the scaling factor that will be used for normalization

This step also generates a plot of the proportion of each library that aligned to either the *S. pombe* or *S. cerevisiae* genome. It's called 'proportion_reads_mapped.png' and can also be found in `logs/`.


**8. Generate coverage tracks for each library scaled by spike-in normalization.**

```bash
# use a for loop to submit each coverage job in parallel
for name in bam/*_sorted.bam; do sbatch scripts/bamCoverage_spikein.sh $name; done
```


**9. Calculate and plot correlation scores.**

```bash
sbatch scripts/multiBigwigSummary_Correlation.sh
```

This script calculates coverage scores in non-overlapping 75 bp bins genome-wide (skipping the mitochondrial genome).

The binned coverage scores are used to calculate Pearson correlation coeeficients between all libraries. The results are plotted as a heatmap in the `correlation/plots` directory. It also outputs the results in tabular (.tab) format in the `correlation/tab` directory if you would like to plot it yourself. Sample code to plot for yourself in python is provided as `scripts/plot_correlation.py`.


**10. Average biological replicates to generate coverage tracks.**

If the correlation between biological replicates looks good, you can average your coverage tracks before plotting your results.

```bash
# calculate average coverage across all replicates in 10 bp windows
sbatch scripts/average_bigwigAverage.sh
```

The averaged files (.bw) can be found in `deeptools/averaged/`. These files can be loaded by IGV to visualize nucleosome dyad coverage at individual loci.

**11. Generate matrices to plot data.**

Using the coverage (.bw) files, generate matrices with which to plot the data. The BED files line 3,087 non-overlapping protein-coding genes up according to their +1 nucleosome.

```bash
# Sorted by gene length. Plotting data is generated for 500 bp upstream and up to 4500 bp downstream. Positions after the end of the gene are filled with nan.
sbatch scripts/computeMatrix_reference.sh

# Sorted by Rpb1 ChIP-seq occupancy in wild-type. Plotting data is generated for 500 bp upstream and 1500 bp downstream of this position.
sbatch scripts/computeMatrix_reference_Rpb1sorted.sh
```

There should now be new .tab files in the `deeptools/averaged/tab` directory that will be used in the next step to plot the results. The compressed matrices are also saved in `deeptools/gz` if you would like to use deeptools to generate plots.

**12. Import and plot data.**

The .tab files can be used to plot the results in python. I do this by moving the files to my local machine, though you could certainly also try to do the plotting on O2. Sample code to plot for yourself in python is provided as `scripts/plot_mnase.py`.

To transfer .tab files to your machine:
```bash
# from a local session
scp -r <YOUR_USER_ID>@transfer.rc.hms.harvard.edu <PATH/TO/DIR/>deeptools/averaged/tab <./PATH/TO/LOCAL/DIR>
```
To fix the first few lines of the files and make them python-readable, run the following code (provided as `scripts/tab_converter.sh` on the directory containing the .tab files:
```bash
for name in <FILE/PTAH/>*.tab; do sed '3s/genes:3087\t/#genes:\n/g' $name > lala.tab && mv lala.tab $name; done
```

