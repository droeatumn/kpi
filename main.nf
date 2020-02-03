#!/usr/bin/env nextflow

/*
 * kpi
 *
 * KIR structural interpretation for raw reads.
 *
 * Predict a (possibly-amiguous) pair of haplotypes given the
 * presence/absence (PA) genotype from one individual and a collection
 * of PA reference haplotypes.
 * 
 * The input is a directory of 2-column tab-separated unix-formatted text 
 * files. Each file maps IDs to one or more files. The first column of each
 * file is the ID of a set of reads (e.g., representing
 * an individual). The second column is non-space-containing path to a 
 * optionally-gzipped fastq or fasta file.
 *
 * The output contains genotype and haplotype-pair predictions.
 *
 * One individual with 60G of gzippped fastq takes about todo minutes 
 * on an 8 CPU computer with 20G memory.
 * 
 * @author Dave Roe
 */

inputSuffix = "txt"
nfKMCForks = 1 // run this many input text files in parallel
params.p = '/opt/kpi/raw/'
params.output = '/opt/kpi/output'
params.id = 'defaultID'
geneProbes  = '/opt/kpi/input/markers'
nfForks = 4 // run this many input text files in parallel
// input: kmc probe txt files
kmcNameSuffix = '_hits.txt'          // extension on the file name
bin1Suffix = 'bin1'
probeFile = '/opt/kpi/input/markers.fasta'
params.haps = '/opt/kpi/input/haps.txt'
params.m = null
//workflow.onComplete { file('work').deleteDir() }

// things that probably won't change per run
resultDir = params.output
haps = params.haps
if(!resultDir.trim().endsWith("/")) {
	resultDir += "/"
}
probeCmd = ""
mapDir = ""
if(params.m != null) {
    probeCmd = "probeFastqsKMC.groovy -m ${params.m} -o . -w ."
    mapDir - params.m
} else if(params.p != null) {
    probeCmd = "probeFastqsKMC.groovy -d ${params.id} -p ${params.p} -o . -w ."
    mapDir - params.p
}
if(!mapDir.trim().endsWith("/")) {
	mapDir += "/"
}
fqsIn = Channel.fromPath(mapDir).ifEmpty { error "cannot find anything in $mapDir" }

/* 
 * @todo handle both input options
 * @todo m option only publishes everything at the end
 */ 
process probeFastqs {
	//container = "droeatnmdp/kpi:latest"
	//publishDir resultDir, mode: 'copy', overwrite: true
//    maxForks 1
	
	input: file(f) from fqsIn
	output:
    	file('*.kmc_*') into kmcdb
	script:
		"""
        ${probeCmd}
		"""
} // probeFastqs

process probeDB {
    //publishDir resultDir, mode: 'copy', overwrite: true

	input: file(kmc) from kmcdb
	output:
		file{ "*_hits.txt"} into filterdb
	
	script:
		"""
        filterMarkersKMC2.groovy -d . -p ${geneProbes} -o . -w .
		"""
		
} // probeFastqs

/*
 * db2Locus
 *
 * Given a kmc output file, bin the hit reads into separate files based on locus.
 * 
 * e.g., ./kmc2Locus2.groovy -j 100a.txt -p kmers.txt -e bin1 -o output
 * 
 * Input files: e.g., 100a.fasta
 * Output files have an extension of 'bin1'.
 * @todo improve the memory usage here
 */
process db2Locus {
  //publishDir resultDir, mode: 'copy', overwrite: true
  maxForks nfForks

  input:
    file(hits) from filterdb.flatMap()
  output:
	file{"*.bin1"} into bin1Fastqs
	val(id) into idc

script: 
    // e.g., gonl-100a.fasta
    // todo: document this
	String dataset
	String idn
	id = hits.name.replaceFirst(kmcNameSuffix, "")
    """
    kmc2LocusAvg2.groovy -j ${hits} -p ${probeFile} -e ${bin1Suffix} -i ${id} -o .

if ls *.bin1 1> /dev/null 2>&1; then
    : # noop
else
    echo "combined: uninterpretable" > "${id}_prediction.txt"

    touch "${id}_uninterpretable.bin1"
fi
    """
} // db2Locus

/*
 * hapInterp
 * 
 * 1) Makes haplotype predictions from PA probes.
 *
 *
 * @todo document
 */
process hapInterp {
  publishDir resultDir, mode: 'copy', overwrite: true

  input:
	file(b1List) from bin1Fastqs
	val(idIn) from idc
  output:
	file{"*_prediction.txt"} into predictionChannel

  script:
    """
    FILES="*.bin1"
    fileList=""
    id=""
    ext="*bin1*"
    for bFile in \$FILES; do
        if [ -s \$bFile ]; then
            if [[ \$bFile == \$ext ]]; then
                id=\$(basename "\$bFile")
                # '%' Means start to remove after the next character;
				# todo: change this to kmcNameSuffix
                id="\${id%%_*}"
                if [ "\$id" == "" ]; then
                    id=\$(basename "\$bFile")
                fi
                #echo \$bFile
                if [ "\$fileList" == "" ]; then
                    :
                else
                    fileList+=","
                fi
                fileList+=\$bFile
            fi
        fi
    done
    outFile=${idIn}
    outFile+="_prediction.txt"
    pa2Haps.groovy -h ${haps} -q "\$fileList" -o "\$outFile"
    """
} // hapInterp

// get the per-sample name
def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  int end = name.indexOf(inputSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + inputSuffix + "'" );
  }
  end = end -1 // Remove the trailing '.'
  return name.substring(start, end)
} // sample
