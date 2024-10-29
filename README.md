
# MNase-seq analysis pipeline

## description



## requirements



## instructions

**1. Clones this repository into a clean directory on O2.**

```bash
# clone the repostiory
git clone https://github.com/jamwarner/mnase_seq.git

#navigate to the newly created directory
cd mnase_seq
```

Run all commands from the base `mnase_seq` directory.

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


