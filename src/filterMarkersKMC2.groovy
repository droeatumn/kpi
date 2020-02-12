#!/usr/bin/env groovy

/*
 * filterMarkersKMC2
 *
 * Given a single or directory of kmc databases (.kmc_pre), 
 * merge and dump each with a kmer db.
 * The output will include the counts in a tab-delimited text file.
 * It is in the format '_hit.txt'.
 * 
 * e.g., filterMarkersKMC2.groovy -d /project/compbioRAID1/davidr/gonl/kmc/db/todo -p $HOME/doc/kir/snp/markers_v1 -o . -w .
 * 
 * For now, place the work and output directories in a separate location from 
 * input.(?)
 * 
 * Requires KMC 3.
 *   http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc
 *   kmc_tools simple gonl-53a ~/doc/kir/snp/25mers_phv2 intersect gonl-53a_probes -ocleft
 *   kmc_tools dump -s gonl-53a_probes gonl-53a_hits.txt
 *
 * @author Dave Roe
 * @todo skip work and output folders?
 * @todo exclude probes that max out count (255)?
 * @todo remove the databases
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3

// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')

OptionAccessor options = handleArgs(args)
if(debugging <= 4) {
    err.println "filterMarkersKMC(options=${options})"
}
probeFile = options.p
TreeSet<String> inFileSet = makeInFileSet(options.d, fileSeparator)
outDir = options.o
if(!outDir.endsWith(fileSeparator)) {
    outDir += fileSeparator
}

// loop through each file within the input directory
TreeSet<String> dbSet = new TreeSet()
// loop through the fasta files
inFileSet.each { f ->
    //err.println f.toString() //todo
    // e.g., 2DS5-2DS1_unique.kmc_pre
    //(shortName, ext) = f.toString().split('\\.?') //todo:document and/or make better
	shortName = f.toString()[0..f.lastIndexOf('\\.')]
    dbSet.add(shortName)
}

//remove(todo) HashMap<String,String> fileMap
//remove(todo) fileMap = inFilesToSet(topdir)
dbSet.each { db ->
    if(debugging <= 4) {
        err.println "processing ${db} ..."
    }

    ret = filterKMC(db, dbSet, db, probeFile, outDir)
} // each subdirectory

err.println "done"

// end main

/*
 * filterKMC
 *
 * Run an KMC filter operation of the two inputs (db and fasta).
 *
 * e.g., todo
 *
 */
void filterKMC(String db, TreeSet dbSet, String id,
			   String probeFile, String outDir) {
    if(debugging <= 1) {
        err.println "filterKMC(db=${db}, dbSet=${dbSet}, id=${id}, probeFile=${probeFile}, outDir=${outDir})"
    }
    // make the command file
    String defFileName = null
    String dbIndex = null // input set index for db
    dbname = db.replace(".kmc_pre", "")
    extension = ".txt"
    dbNoUnique = db.replace(".kmc_pre", "_hits")
    dbNoUnique = dbNoUnique.replaceFirst(".*${fileSeparator}", "")
    resultDb = outDir + dbNoUnique
	resultFile = outDir + dbNoUnique + extension
	// e.g., kmc_tools simple gonl-53a ~/doc/kir/snp/25mers_phv2 intersect gonl-53a_probes -ocleft
    cmd = ["kmc_tools", "-hp", "simple", dbname, probeFile, "intersect", 
           resultDb, "-ocleft"]
    if(debugging <= 3) {
        err.println cmd
    }
    ret = cmd.execute()
    ret.waitFor()
    retVal = ret.exitValue()
    if(debugging <= 2) {
        err.println "filterKMC: returned ${retVal}"
    }

	// e.g., kmc_dump -s gonl-53a_probes gonl-53a_hits.txt
	cmd = ["kmc_dump", resultDb, resultFile]
    if(debugging <= 3) {
        err.println cmd
    }
    ret = cmd.execute()
    ret.waitFor()
    retVal = ret.exitValue()
    if(debugging <= 2) {
        err.println "filterKMC: returned ${retVal}"
    }

    // todo: not working
    pwd = System.getProperty("user.dir");
	rmFile1 = pwd + "/${dbname}.kmc_pre"
	rmFile2 = pwd + "/${dbname}.kmc_suf"
    try { 
        rmFile1S = Files.readSymbolicLink(rmFile1)
        rmFile2S = Files.readSymbolicLink(rmFile2)
	    cmd = ["rm", "-f", rmFile1S.toString(), rmFile2S.toString()]
        if(debugging <= 3) {
            err.println cmd
        }
        ret = cmd.execute()
        ret.waitFor()
        retVal = ret.exitValue()
        if(debugging <= 2) {
            err.println "filterKMC: returned ${retVal}"
        }
    } catch(NotLinkException) {
        err.println "not sym link"
	    cmd = ["rm", "-f", rmFile1.toString(), rmFile2.toString()]
        if(debugging <= 3) {
            err.println cmd
        }
        ret = cmd.execute()
        ret.waitFor()
        retVal = ret.exitValue()
        if(debugging <= 2) {
            err.println "filterKMC: returned ${retVal}"
        }
    }
    err.println "droe2"//todo
    if(debugging <= 1) {
        err.println "filterKMC: return"
    }
} // filterKMC

/*
 * inFilesToSet
 * 
 * Convert the input directory and file list to a set of Strings
 * containing the full file names.
 * 
 * Returns a map of short name to fully qualified name
 */ 
HashMap<String,String> inFilesToSet(File dir) {
    HashMap<String,String> fileMap = new HashSet()
    fp = ~/(.*).fasta(.*)/
    // loop through the fasta files
    dir.eachFileMatch(fp) { inFile ->
        if(debugging <= 1) {
            err.println "inFile = ${inFile}"
        }
        inFileName = inFile.getName()
        inFullName = inFile
        fileMap[(inFileName)] = inFullName
    } // each file
    return fileMap
}   //inFilesToSet

/*
 * Make a set of file names given a directory or file name.
 *
 */
TreeSet<String> makeInFileSet(String inDir, String fileSeparator) {
	if(debugging <= 2) { 
		err.println "inDir=${inDir}"
	}
	TreeSet<String> ret = new TreeSet()
	f = new File(inDir)
	if(f.isDirectory()) { 
		fp = ~/(.*).kmc_pre/
		f.eachFileMatch(fp) { f2 ->
            name = f2.toString()
			//name = inDir + fileSeparator + name
			name = name.replaceFirst("././", "")
			ret.add(name)
		}
	} else {
        name = inDir
		ret.add(name)
	}
	return ret
} // makeInFileSet

/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'filterMarkersKMC.groovy [options] ',
                                    header:'Options:')
    cli.help('print this message')

    cli.d(longOpt:'base directory for input file pattern', args:1, argName:'dir',
          'directory location',
      required:true)
    cli.o(longOpt:'single file or directory to .kmc_pre file(s)', args:1,
		  argName:'out', 'output directory', required: true)
	cli.p(longOpt:'KMC db for probes', args:1, argName:'probes', 
          'probe db', required: true)
    cli.w(longOpt:'work directory', args:1, argName:'work', 
          'work directory', required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
