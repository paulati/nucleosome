## Improving genome-wide mapping of nucleosomes in <i>Trypanosome cruzi</i>

<a href="https://doi.org/10.1101/2021.07.02.450927" target="_blank">bioRxiv pre-print</a>
<a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0293809" target="_blank">paper</a>


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

### Pipeline

Available steps: 
<ul>
   <li>preprocess</li>
   <li>align</li>
   <li>convert_to_bam</li>
   <li>plot</li>
   <li>summary_tables</li>
   <li>convert_to_bigwig</li>
</ul>

#### Suggested Pipeline

1. <code>main('preprocess')</code>
2. <code>main('align')</code>
3. <code>main('summary_tables')</code>
4. <code>main('convert_to_bam')</code>
5. <code>main('plot')</code>
6. <code>main('convert_to_bigwig')</code>

