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

params.home = baseDir
home = params.home + '/'
params.raw = home + '/raw/'
raw = params.raw + "/"
params.output = home + '/output/'
params.id = "defaultID"
markerDBPrefix  = 'markers'
markerDBSuf = file("${home}/input/markers.kmc_suf")
markerDBPre = file("${home}/input/markers.kmc_pre")
// input: kmc probe txt files
kmcNameSuffix = '_hits.txt'          // extension on the file name
bin1Suffix = 'bin1'
markerFile = file("${home}/input/markers.fasta")
params.haps = file("${home}/input/haps.txt")
params.map = null
params.l = null // logging level (1=most to 5=least)
params.container = "droeatumn/kpi:latest"
makeKmerDBFile = file("${home}/src/probeFastqsKMC.groovy")
queryDBFile = file("${home}/src/filterMarkersKMC2.groovy")
db2LocusFile = file("${home}/src/kmc2LocusAvg2.groovy")
pa2HapsFile = file("${home}/src/pa2Haps.groovy")
srcDir = pa2HapsFile.parent
params.nocontainer = "null"

// things that probably won't change per run
fileSeparator = "/"
resultDir = params.output
haps = params.haps
if(!resultDir.endsWith("/")) {
	resultDir += "/"
}
logIn = ""
if(params.l != null) {
    logIn = "-l ${params.l}"
}

inputSuffix = "*.txt"   // for --m
inOption = ""  // input type for probeFastqsKMC.groovy
dOption = "" // d option probeFastqsKMC.groovy
kpiIn = null
mapDir = home
if(params.map != null) {
    inOption = "-m"
    mapFile = params.map
    mapDir = file(params.map).parent
    kpiIn = Channel.fromPath(mapFile).ifEmpty { error "cannot find file $mapFile" }
} else if(raw != null) {
    inOption = "-p"
    dOption = "-d " + params.id
    kpiIn = Channel.fromPath(raw).ifEmpty { error "cannot find fastq/fastq in $mapDir" }
//        fqsIn = Channel.fromPath(mapDir).map{ file -> tuple(file.baseName, file) }.ifEmpty { error "cannot find fastq/fastq in $mapDir" }
//    fqsIn = Channel.fromPath(["${mapDir}*.fq", "${mapDir}*.fastq","${mapDir}*.fq.gz", "${mapDir}*.fastq.gz", "${mapDir}*.fa", "${mapDir}*.fasta","${mapDir}*.fa.gz", "${mapDir}*.fasta.gz"] ).ifEmpty { error "cannot find fastq/fastq in $mapDir" }
}

/* 
 * @todo handle both input options
 * @todo m option only publishes everything at the end
 */ 
process makeKmerDB {
    if(params.nocontainer == "null") { 
	    container = params.container
    }
//    publishDir resultDir, pattern: '*.kmc_*', mode: 'copy', overwrite: true
// 	  publishDir resultDir, pattern: '*.log', mode: 'copy', overwrite: true
//    errorStrategy 'ignore'
//    validExitStatus 0,1
    
	input:
      path(f) from kpiIn
      path(makeKmerDBFile)
      path(markerDBSuf)
      path(markerDBPre)
      path(mapDir)
      path(queryDBFile)
      val(markerDBPrefix)
      env(JAVA_OPTS) from ('-Xms2G -Xmx50G')
	output:
  	  file{ "*_hits.txt"} into filterdb
//      file('*.log') optional true into kmcdbLog
	script:
	"""
    ./${makeKmerDBFile} ${dOption} ${inOption} ${f} ${logIn} -o . -w . 2> probeFastqsKMC.log
    ./${queryDBFile} -d . -p ${markerDBPrefix} -o . -w . 2> filterMarkersKMC2.log		
    find . -type f -size 0 | xargs rm # remove 0 length files, especially hits.txt files
    """
} // makeKmerDB

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
    if(params.nocontainer == "null") { 
	    container = params.container
    }
    //publishDir resultDir, mode: 'copy', overwrite: true

    input:
      file(hits) from filterdb.flatMap()
      file(markerFile)
      file(db2LocusFile)
      env(JAVA_OPTS) from ('-Xms2G -Xmx50G')
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
    ./${db2LocusFile} -j ${hits} -p ${markerFile} -e ${bin1Suffix} -i ${id} -o . 2> kmc2LocusAvg2.log

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
    if(params.nocontainer == "null") { 
	    container = params.container
    }
    publishDir resultDir, mode: 'copy', overwrite: true

    input:
  	  file(b1List) from bin1Fastqs
  	  val(idIn) from idc
      file(haps)
      file(pa2HapsFile)
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
    ./${pa2HapsFile} -h ${haps} -q "\$fileList" -o "\$outFile" 2> pa2Haps.log
    """
} // hapInterp
