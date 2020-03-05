#!/usr/bin/env nextflow
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
}

/* 
 * @todo handle both input options
 * @todo m option only publishes everything at the end
 */ 
process testKMC {
    if(params.nocontainer == "null") { 
	    container = params.container
    }
    errorStrategy 'ignore'
    publishDir resultDir, mode: 'copy', overwrite: true
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
  	  path{ "*.log"} into filterdb


	script:
	"""
for file in ${f}/*.gz
do
    base=\$(basename "\$file")
    echo "zcat \$file > \$base.fastq ... " > \$base.log
    if zcat \$file > \$base.fastq  2>> \$base.log; then
        echo "zcat \$file worked" >> \$base.log
        echo "kmc \$file \$base ..." >> \$base.log
        if kmc -k25 -ci3 \$file \$base . 2>> \$base.log; then
            echo "kmc \$file worked" >> \$base.log
        else
            echo "*ERROR*: kmc \$file failed" >> \$base.log
        fi
    else
        echo "*ERROR*: zcat \$file failed" >> \$base.log
    fi

    if du -hs \$file \$base.fastq \$base.kmc_pre >> \$base.log; then
        echo "du \$file worked"
    else
        echo "du \$file failed"
    fi
    rm -f \$base.fastq \$base.kmc_suf \$base.kmc_pre

    o=`env | grep -i java`
    echo \$o &>>\$base.log
    o=`cat /proc/meminfo`
    echo \$o &>>\$base.log

done
    """
} // testKMC
