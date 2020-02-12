# KPI
main.nf makes the predictions.

<h2>Dependancies</h2>
Install Nextflow, Docker, and Git.
Create accounts in GitHub and Docker Hub.
Add 'docker.enabled = true' and 'docker.fixOwnership = true' to your Nexflow
configuration (e.g., $HOME/.nextflow/config). Make sure Docker is running
and you are logged in to Docker Hub.

<h2>Running</h2>
There are two input options.<br>
1. An ID along with a folder of fasta or fastq files, optionally gzipped. (--raw)<br>
2. A two-column text file, where the first column is an ID, and the second column is a fully-qualified local path to a fasta or fastq file. Each ID may have multiple rows. (--map)<br>
<br>
Option 1 is more efficient with respect to disk space. <br>

<b>Output</b> <br>
For each input ID, an output text file will be created named '<id>_prediction.txt'. Each ID's output file contains a header line and a second line with the haplotype pair predictions and gene predictions. <br>
Each haplotype within a pair is separated by a '+'. If the prediction is ambiguous, each pair of haplotypes is separated by '|'.
    e.g., <br> 'cA01&tilde;tA01+cA01&tilde;tB01|cA01&tilde;tA01+cB05&tilde;tB01|cA01&tilde;tB01+cB05&tilde;tA01' means haplotype <br>'cA01&tilde;tA01 and cA01&tilde;tB01' or 'cA01&tilde;tA01 and cB05&tilde;tB01' or 'cA01&tilde;tB01 and cB05&tilde;tA01'. <br>

The reference haplotypes are defined at https://github.com/droeatumn/kpi/blob/master/input/haps.txt <br>

<b>Running</b><br>
Use the 'raw' arguments to indicate the input directory, and 'output' to indicate the directory to put the output. The defaults are 'raw' and 'output' under the directory KPI was installed.<br>

Option 1: Provide and ID (--id) and a folder (--raw) with its raw data<br>
<code>    ./main.nf --id ID --raw inDir --output outDir</code><br>
e.g., <code>    ./main.nf --id id1 --raw ~/input --output ~/output</code><br>

Option 2: Provide a file with a map (--map) from IDs to their raw data<br>
<code>    ./main.nf --map mapFile.txt --output outDir</code><br>
e.g., <code>    ./main.nf --map ~/input/idstoRaw.txt --output ~/output</code><br>

<b>Examples using data in the image, so no input is required.</b><br>
Example 1: cA01&tilde;tA01+cB01&tilde;tB01 with --raw.<br>
Run the following command for an example of interpreting synthetic reads created from sequences with Genbank IDs KP420439 and KP420440 (https://www.ncbi.nlm.nih.gov/nuccore/KP420439 and https://www.ncbi.nlm.nih.gov/nuccore/KP420440)). These two haplotypes contain all the genes, so the haplotype predictions are very ambiguous. <br>

<code>    ./main.nf --id ex1 --raw ~/git/kpi/input/example1 --output ~/output</code><br>

Example 2: cA01&tilde;tA01+cA01&tilde;tB01 with --map. <br>
Run the following command for an example of interpreting synthetic reads created from sequences with Genbank IDs KP420439 and KU645197 (https://www.ncbi.nlm.nih.gov/nuccore/KP420439 and https://www.ncbi.nlm.nih.gov/nuccore/KU645197).<br>

<code>    ./main.nf --map ~/git/kpi/input/example2/example2.txt --output ~/output</code>

Example 3: example 1 and example 2 with --map in the same file. <br>
Run example 1 and example 2 from one input map file. It will create example1_prediction.txt and example2_prediction.txt.<br>

<code>    ./main.nf --map ~/git/kpi/input/example1-2.txt --output ~/output</code>

<b>Miscellaneous</b><br>
Hardware
For targeted sequencing, kpi requires at least 5G RAM total and 1G temp disk space/ID. For WGS, it requires 25G RAM total and 15G temp disk space/ID. It will scale to the number of CPUs available, with 6-12 being most efficient in general for WGS.

Reference
A manuscript is under preparation.
