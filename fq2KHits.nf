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

nfForks = 1 // run this many input text files in parallel
mapSuffix = "txt"
params.input = '/opt/kpi/raw/'
params.output = '/opt/kpi/output'
geneProbes  = '/opt/kpi/input/geneHapSigMarkers_v1-wRc'

mapDir = params.input
resultDir = params.output
mapPath = mapDir + '*' + mapSuffix
fqsIn = Channel.fromPath(mapPath).ifEmpty { error "cannot find any ${mapSuffix} files in ${mapDir}" }.map { path -> tuple(sample(path), path) }

process probeFastqs {
	//publishDir resultDir, mode: 'copy', overwrite: true
    maxForks nfForks

	input: set s, file(f) from fqsIn
	output:
		set s, file('*.kmc_*') into kmcdb
	script:
		"""
        probeFastqsKMC.groovy -m ${f} -o . -w .
		"""
		
} // probeFastqs

process probeDB {
	publishDir resultDir, mode: 'copy', overwrite: true

	input: set s, file(fList) from kmcdb
	output:
		set s, file{ "*_hits.txt"} into filterdb
	
	script:
		"""
        echo ${fList}
        filterMarkersKMC2.groovy -d . -p ${geneProbes} -o . -w .
		"""
		
} // probeFastqs

// get the per-sample name
def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  int end = name.indexOf(mapSuffix)
  if ( end <= 0 ) {
    throw new Exception( "Expected file " + name + " to end in '" + mapSuffix + "'" );
  }
  end = end -1 // Remove the trailing '.'
  return name.substring(start, end)
} // sample
