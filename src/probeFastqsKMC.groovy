#!/usr/bin/env groovy

/*
 * probeFastqsKMC
 *
 * Given a set of FASTQ files, create a KMC 3 database for each ID.
 *
 * 
 * Requires KMC 3.
 *   http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc
 *   'kmc' to build database
 *   It requires the input format to be specified.
 *     -f<a/q/m/bam/kmc> - input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ
 *   e.g., kmc -k25 -ci2 -fq @./gonl-52b-cmd.txt ./gonl-52b work3
 *
 * @author Dave Roe
 * @todo handle blank lines in the input text files
 */

import groovy.io.*
import groovy.io.FileType
import org.apache.commons.cli.CliBuilder
import org.apache.commons.cli.OptionAccessor
# newer?
#import groovy.cli.commons.CliBuilder
#import groovy.cli.commons.OptionAccessor

// things that may change per run
debugging = 1 // TRACE=1, DEBUG=2, INFO=3
kmerSize = "25"
minKmers = "3" // reads hit less than this will be ignored

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')

OptionAccessor options = handleArgs(args)
if((options.l != null) && (options.l != false)) {
    debugging = options.l.toInteger()
}

// make list of fastq files for every individual
String id = options.d
String path = options.p
String mpath = options.m
String inFileType = options.t
if((path != null) && (path != "") && (path != "false")) {
    mpath = ""
} else {
    path = ""
}
if(debugging <= 4) {
    err.println "kmerSize=${kmerSize} minKmers=${minKmers}"
    err.println "id=${id} output=${options.o}"
    err.println "inFileType=${inFileType}"
}    
HashMap<String,ArrayList<String>> fqMap = loadFqMap(id, path, mpath)
if(debugging <= 2) {
    err.println "${fqMap.keySet().size()} IDs in the fastq map"
    firstKey = fqMap.keySet().iterator().next()
    err.println "${fqMap[firstKey].size()} fastq files for ${firstKey}"
}

probeHits = probeReads(options.o, options.w, inFileType, fqMap, kmerSize, minKmers)
err.println "done"
// end main

/*
 * loadFqMap
 *
 * @param fqMapFileName file containing a tab-delimited mapping between ids and fastq files
 *
 * if fpath is not null, append the file name to the path and use it;
 * if it is null, use all the files in the directory (fpath);
 * the id is take from the file or the directory
 */
HashMap<String,ArrayList<String>> loadFqMap(String inid,
											String fpath, String mapFile) { 
	// return value
    HashMap<String,ArrayList<String>> fqMap = new HashMap()

	if((fpath != null) && (fpath != "")) {
		new File(fpath).eachFileRecurse(FileType.FILES)  { inFile ->
            if(inFile.name.endsWith(".fq") || inFile.name.endsWith(".fq.gz") ||
               inFile.name.endsWith(".fastq") || inFile.name.endsWith(".fastq.gz") ||
               inFile.name.endsWith(".fa") || inFile.name.endsWith(".fa.gz") ||
               inFile.name.endsWith(".fasta") || inFile.name.endsWith(".fasta.gz") ||
               inFile.name.endsWith(".bam")) {
			    if(inid == null) {
				    // get the id from the directory name
				    id = fpath.split(System.getProperty('file.separator'))[-1]
			    } else {
                    id = inid
                }
			    //fileName = fpath + fileSeparator + inFile
			    fileName = inFile
			    idList = fqMap[id]
			    if(idList != null) { 
				    idList.add(fileName)
			    } else { 
				    ArrayList<String> l = new ArrayList()
				    l.add(fileName)
				    fqMap[id] = l
			    }
	        } // if correct file name
        } // each file
	} else {
		f = new File(mapFile)
		FileReader probeReader = new FileReader(f)
		probeReader.eachLine { line ->
			if(debugging <= 1) {
				//err.println line
			}
			(id, shortFileName) = line.split('\t')
			idList = fqMap[id]
			if(idList != null) { 
				idList.add(shortFileName)
			} else { 
				ArrayList<String> l = new ArrayList()
				l.add(shortFileName)
				fqMap[id] = l
			}    
		} // each line of file
		
		probeReader.close()
	}

    return fqMap
} // loadFqMap

/*
 * probeReads
 *
 * e.g., kmc -k25 -ci2 -fq @output/gonl-157a-cmd.txt output/gonl-157a work
 * 
 * For each individual, build a database from their fastq files.
 * 
 */
def void probeReads(String outputDir, String workDir, String inFileType,
                    HashMap<String,ArrayList<String>> fqMap,
                    String kmerSize, String minKmers) {
    if(debugging <= 1) {
        err.println "probeReads(workDir=${workDir}, outputDir=${outputDir})"
    }

    fqMap.keySet().sort().each { id ->
        outFile = outputDir + fileSeparator + id
        // make the command file
        defFileName = outputDir + fileSeparator + id + "-cmd.txt"
        outWriter = new PrintWriter(new File(defFileName).newOutputStream(),
                                    true)
        fqMap[id].each { fq ->
            outWriter.println fq
        }
        outWriter.close()

		//https://github.com/refresh-bio/KMC/issues/40 // todo (remove, na?)
        cmd = ["kmc", "-k${kmerSize}", "-ci${minKmers}", "-${inFileType}",
               "@${defFileName}", outFile, workDir]
        if(debugging <= 3) {
            err.println cmd
        }
        ret = cmd.execute()
        ret.waitFor()
        retVal = ret.exitValue()
        if(retVal) {
            err.println "makeKMCdb: *error* running ${cmd}"
        }

        if(debugging <= 2) {
            err.println "probeReads: returned ${retVal}"
        }
    } // each id
    if(debugging <= 1) {
        err.println "probeReads: return"
    }
} // probeReads

        
/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'probeFastqsKMC.groovy [options] ',
      header:'Options:')
    cli.help('print this message')

    cli.m(longOpt:'fastq name map or directory name (one ID only)', args:1,
		  argName:'map', 'fastq map', required: false)
    cli.d(longOpt:'ID for output', args:1,
		  argName:'id', 'id', required: false)
    cli.p(longOpt:'path to the sequences files', args:1,
		  argName:'path', 'fastq path', required: false)
    cli.l(longOpt:'log level, 1(most)-5(least)', args:1,
		  argName:'logging', 'logging level', required: false)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
		  'output directory', required: true)
    cli.t(longOpt:'input file type (for KMC)', args:1, argName:'type', 
		  'input in FASTA format (-fa), FASTQ format (-fq), multi FASTA (-fm) or BAM (-fbam) or KMC(-fkmc); default: FASTQ', required: true)
    cli.w(longOpt:'work directory', args:1, argName:'work', 
          'work directory', required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
