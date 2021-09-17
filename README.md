## Experimental and Bioinformatic upgrade for genome-wide mapping of nucleosomes in Trypanosoma cruzi

These scripts require already installed software listed <a href="https://github.com/paulati/nucleosome/blob/master/tools/tools.txt">here</a>

Custom configuration values should be set <a href="https://github.com/paulati/nucleosome/blob/master/source/common/config.yaml">here</a>

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
