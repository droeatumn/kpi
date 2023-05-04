# KPI
main.nf makes the predictions.

<h2>Dependencies</h2>
Install Java, Groovy, Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow
configuration (e.g., $HOME/.nextflow/config). Make sure Docker is running
and you are logged in to Docker Hub.

<h2>Running</h2>
<b>Input</b> <br>
There are two input options.<br>
1. An ID along with a folder of fasta or fastq files, optionally gzipped. (--raw and --id)<br>
2. A two-column text file, where the first column is an ID, and the second column is a path to a fasta or fastq file (--map). Each ID may have multiple rows. The paths to the files be absolute or relative, but the files must be in the same directory as the map file or under it. If using relative paths, the paths must start with the _parent_ folder of the map file.<br>
<br>
Option 1 is more efficient with respect to disk space. <br>
<br>
<b>Output</b> <br>
For each input ID, an output text file will be created named '<id>_prediction.txt'. Each ID's output file contains a header line and a second line with the haplotype pair predictions and gene predictions. <br>
Each haplotype within a pair is separated by a '+'. If the prediction is ambiguous, each pair of haplotypes is separated by '|'.
    e.g., <br> 'cA01&tilde;tA01+cA01&tilde;tB01|cA01&tilde;tA01+cB05&tilde;tB01|cA01&tilde;tB01+cB05&tilde;tA01' means haplotype <br>'cA01&tilde;tA01 and cA01&tilde;tB01' or 'cA01&tilde;tA01 and cB05&tilde;tB01' or 'cA01&tilde;tB01 and cB05&tilde;tA01'. <br>

The reference haplotypes are defined at https://github.com/droeatumn/kpi/blob/master/input/haps.txt <br>

<b>Running</b><br>
Use 'raw' to indicate the input directory, and 'output' to indicate the
directory to put the output. The defaults are 'raw' and 'output' under the
location where KPI was pulled.<br>
Use 'filetype' to indicated the input type; default is 'fq' (FASTQ).<br>
<code>    f<a/q/m/bam/kmc> - input in FASTA format (fa), FASTQ format (fq), multi FASTA (fm) or BAM (fbam) or KMC(fkmc); default: FASTQ</code><br>

Option 1: Provide and ID (--id) and a folder (--raw) with its raw data<br>
<code>    ./main.nf --id ID --raw inDir --output outDir --filetype fq</code><br>
e.g., <code>    ./main.nf --id id1 --raw ~/input --output ~/output</code><br>

Option 2: Provide a file with a map (--map) from IDs to their raw data<br>
<code>    ./main.nf --map mapFile.txt --output outDir --filetype fq</code><br>
e.g., <code>    ./main.nf --map ~/input/idstoRaw.txt --output ~/output</code><br>
In this example the path to files in idstoRaw.txt are somewhere under ~/input/.

<b>Example using data in the image, so no input is required.</b><br>
Example 1: cA01&tilde;tA01+cB01&tilde;tB01 with --raw.<br>
Run the following command for an example of interpreting synthetic reads created from sequences with Genbank IDs KP420439 and KP420440 (https://www.ncbi.nlm.nih.gov/nuccore/KP420439 and https://www.ncbi.nlm.nih.gov/nuccore/KP420440)). These two haplotypes contain all the genes, so the haplotype predictions are very ambiguous. <br>

<code>    ./main.nf --id ex1 --raw ~/git/kpi/input/example1 --output ~/output</code><br>
<br>
To run another example, replace 'example1' with 'example2'.<br>

Example 2: cA01&tilde;tA01+cA01&tilde;tB01 with --map and --id.<br>
Run the following command for an example of interpreting synthetic reads created from sequences with Genbank IDs KP420439 and KU645197 (https://www.ncbi.nlm.nih.gov/nuccore/KP420439 and https://www.ncbi.nlm.nih.gov/nuccore/KU645197)).<br>

<code>    ./main.nf --id ex2 --map ~/git/kpi/input/example2/example2.txt --output ~/output</code><br>
<br>
To run another example, replace 'example2' with 'example1'.<br>

Example 3: combine Example 1 and 2 with --map and --id.<br>
<code>    ./main.nf --id ex12 --map ~/git/kpi/input/example1-2.txt --output ~/output</code><br>
<br>

<b>Miscellaneous</b><br>
Hardware<br>
For targeted sequencing, kpi requires approximately 4 CPU, 8G RAM and 20G disk space. For WGS, it requires around 13 CPU, 16G RAM total and 100G temp disk space.<br>
<br>
Raw data<br>
The software assumes average coverage for both chromosomes is less than 255. If this is not the case for your data, please downsample before running. Support for high coverage data is a future enhancement.<br>
<br>
Containers<br>
To run without a container, use the --nocontainer parameter. To use a
container other than the default (droeatumn/kpi:latest), use the --container parameter.
<br><br>
To run in a self-contained environment with the --id parameter. Replace 'inDir' and 'outDir'.<br>
<code>docker run --rm -it -v inDir:/opt/kpi/raw/ -v outDir:/opt/kpi/output/ droeatumn/kpi:latest /opt/kpi/main.nf --id <output ID></code><br>
Or<br>
<code>docker run --rm -it -v $PWD/output:/opt/kpi/output droeatumn/kpi:latest /opt/kpi/main.nf --id ex1 --raw /opt/kpi/input/example1/</code><br>
Or<br>
<code>docker run --rm -it -v $PWD/output:/opt/kpi/output  droeatumn/kpi:latest /opt/kpi/main.nf --map /opt/kpi/input/example1/example1.txt</code><br>
Or, if your bam file (for one individual) is locally in ~/data<br>
<code>docker run --rm -it -v ~/data:/opt/kpi/raw -v $PWD/output:/opt/kpi/output  droeatumn/kpi:latest /opt/kpi/main.nf --filetype fbam --id testid</code><br>
Or, if a map to the bam files locally withing ~/data<br>
<code>docker run --rm -it -v ~/data:/opt/kpi/raw -v $PWD/output:/opt/kpi/output  droeatumn/kpi:latest /opt/kpi/main.nf --filetype fbam --map /opt/kpi/raw/map.txt</code><br>
<br>
<b>Citation</b><br>
Roe D, Kuang R. Accurate and Efficient KIR Gene and Haplotype Inference From Genome Sequencing Reads With Novel K-mer Signatures. Front Immunol (2020) 11:583013. (https://doi.org/10.3389/fimmu.2020.583013)
