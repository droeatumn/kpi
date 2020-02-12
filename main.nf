#!/usr/bin/env nextflow

/*
 * KPI
 *
 * KIR structural interpretation for raw reads.
 *
 * Predict a (possibly-amiguous) pair of haplotypes given the
 * presence/absence (PA) genotype from one individual and a collection
 * of PA reference haplotypes.
 * 
 * @author Dave Roe
 * @todo test that required inputs are submitted
 * @todo combine probeDB and db2Locus so the kmc files can be deleted
 * @todo coverage of 255 or higher is considered off kir
 */

params.base = baseDir
base = params.base + '/'
params.raw = base + '/raw/'
raw = params.raw + "/"
params.output = base + '/output/'
params.id = ""
markerDBPrefix  = 'markers'
markerDBSuf = file("${baseDir}/input/markers.kmc_suf")
markerDBPre = file("${baseDir}/input/markers.kmc_pre")
// input: kmc probe txt files
kmcNameSuffix = '_hits.txt'          // extension on the file name
bin1Suffix = 'bin1'
markerFile = file("${baseDir}/input/markers.fasta")
params.haps = file("${baseDir}/input/haps.txt")
params.map = null
params.l = null // logging level (1=most to 5=least)

// things that probably won't change per run
fileSeparator = "/"
resultDir = params.output
haps = params.haps
if(!resultDir.endsWith("/")) {
	resultDir += "/"
}
probeCmd = ""
mapDir = ""
logIn = ""
if(params.l != null) {
    logIn = "-l ${params.l}"
}

inputSuffix = "*.txt"   // for --m
inOption = ""  // input type for probeFastqsKMC.groovy
dOption = "" // d option probeFastqsKMC.groovy
kpiIn = null
if(params.map != null) {
//    probeCmd = "probeFastqsKMC.groovy -m ${params.map} ${logIn} -o . -w . ${logOut}"
    inOption = "--m"
    mapFile = params.map
    kpiIn = Channel.fromPath(mapFile).ifEmpty { error "cannot find file $mapFile" }
} else if(raw != null) {
//    probeCmd = "${baseDir}/src/probeFastqsKMC.groovy -d ${params.id} -p ${raw} ${logIn} -o . -w . ${logOut}"
    inOption = "--p"
    dOption = "--d"
    if(params.id == null) {
        params.id = "defaultID"
    }
    mapDir = raw
    kpiIn = Channel.fromPath(mapDir).ifEmpty { error "cannot find fastq/fastq in $mapDir" }
//        fqsIn = Channel.fromPath(mapDir).map{ file -> tuple(file.baseName, file) }.ifEmpty { error "cannot find fastq/fastq in $mapDir" }
//    fqsIn = Channel.fromPath(["${mapDir}*.fq", "${mapDir}*.fastq","${mapDir}*.fq.gz", "${mapDir}*.fastq.gz", "${mapDir}*.fa", "${mapDir}*.fasta","${mapDir}*.fa.gz", "${mapDir}*.fasta.gz"] ).ifEmpty { error "cannot find fastq/fastq in $mapDir" }
}

/* 
 * @todo handle both input options
 * @todo m option only publishes everything at the end
 */ 
process makeKmerDB {
	container = "droeatumn/kpi:latest"
//    publishDir resultDir, pattern: '*.kmc_*', mode: 'copy', overwrite: true
// 	  publishDir resultDir, pattern: '*.log', mode: 'copy', overwrite: true
//    errorStrategy 'ignore'
//    validExitStatus 0,1
    
	input:
      path(f) from kpiIn
	output:
      file('*.kmc_*') optional true into kmcdb
//      file('*.log') optional true into kmcdbLog
	script:
		"""
        ${baseDir}/src/probeFastqsKMC.groovy ${dOption} ${params.id} ${inOption} ${f} ${logIn} -o . -w . 2> probeFastqsKMC.log
    #    ${probeCmd}
		"""
} // makeKmerDB

process queryDB {
    container = "droeatumn/kpi:latest"
    //publishDir resultDir, mode: 'copy', overwrite: true

	input:
      file(kmc) from kmcdb
      file(markerDBSuf)
      file(markerDBPre)
      val(markerDBPrefix)
	output:
		file{ "*_hits.txt"} into filterdb
	
	script:
		"""
        ${baseDir}/src/filterMarkersKMC2.groovy -d ${kmc[0]} -p ${markerDBPrefix} -o . -w . 2> filterMarkersKMC2.log
		"""
		
} // queryDB

/*
 * db2Locus
 *
 * Given a kmc output file, bin the hit reads into separate files based on locus.
 *
 * Input files: e.g., 100a.fasta
 * Output files have an extension of 'bin1'.
 * @todo improve the memory usage here
 */
process db2Locus {
    container = "droeatumn/kpi:latest"
    //publishDir resultDir, mode: 'copy', overwrite: true

    input:
      file(hits) from filterdb.flatMap()
      file(markerFile)
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
    ${baseDir}/src/kmc2LocusAvg2.groovy -j ${hits} -p ${markerFile} -e ${bin1Suffix} -i ${id} -o . 2> kmc2LocusAvg2.log

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
    container = "droeatumn/kpi:latest"
    publishDir resultDir, mode: 'copy', overwrite: true

    input:
  	  file(b1List) from bin1Fastqs
  	  val(idIn) from idc
      file(haps)
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
    ${baseDir}/src/pa2Haps.groovy -h ${haps} -q "\$fileList" -o "\$outFile" 2> pa2Haps.log
    """
} // hapInterp
