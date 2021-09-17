## Experimental and Bioinformatic upgrade for genome-wide mapping of nucleosomes in Trypanosoma cruzi

<a href="https://doi.org/10.1101/2021.07.02.450927" target="_blank">bioRxiv pre-print</a>

These scripts require already installed software listed <a href="https://github.com/paulati/nucleosome/blob/master/tools/tools.txt" target="_blank">here</a>

Custom configuration values should be set <a href="https://github.com/paulati/nucleosome/blob/master/source/common/config.yaml" target="_blank">here</a>

### Repository outline

The <code>source</code> directory contains
<ul>
<li><code>main.R</code> - Entry point</li>

<li><code>pipeline</code> - scripts for each step involved in nucleosome data analysis</li>
	
<li><code>common</code> - configuration and helpers scripts</li>
</ul>
The <code>data</code> directory contains
<ul>
	<li><code>genomes</code> - genomes fasta files</li>	
	<li>overrepresented sequences in analyzed samples in fasta format</li>	
	<li>UTRme output defining positions for value of plot2DO --sites</li>
</ul>  
The <code>output</code> directory contains
<ul>
	<li><code>qc_output</code> - fastqc output for each sample</li>	
	<li><code>alignment</code> - summary of alignments for each strategy (bowtie2, hisat2) and each sample</li>	
	<li><code>alignment_summary</code> - csv files for consolidated <code>alignment</code> data</li>	
	<li><code>plots</code> - plot2DO output for each strategy (bowtie2, hisat2) and each sample</li>	
</ul>  

