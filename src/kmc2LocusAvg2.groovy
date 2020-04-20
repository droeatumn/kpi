#!/usr/bin/env groovy

/*
 * kmc2LocusAvg2
 *
 * Convert a kmc output file to probe hits per gene -- via average.
 * Only the refence markers are considered from the genotype markers.
 *
 * input: tsv of probe sequence to count
 * output: fasta file per gene containing sequences that hit (e.g. '2DL4.bin1')
 *   The loci with '/' in their names is converted to 'z'
 *     e.g., 7/24 -> 7z24 (e.g., 53a_7z24.bin1)
 *
 * usage: kmc2LocusAvg2.groovy [options]
 * Options:
 *  -e,--extension <extension>               extension for output files
 *  -help                                    print this message
 *  -i,--ID <id>                             ID for the individual
 *  -j,--kmc probe results <kmc>             input kmc counts
 *  -o,--directory to put the output <out>   output directory
 *  -p,--probes <probes>                     input probes
 *
 * Requires
 *   Apache Commons Math
 *
 * @author Dave Roe
 * @todo this is very slow
 */

import groovy.io.*
import groovy.util.CliBuilder.*
import groovy.util.OptionAccessor
// http://apache.mirrors.pair.com//commons/math/binaries/commons-math3-3.6.1-bin.tar.gz
import org.apache.commons.math3.stat.StatUtils

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
ignoreNonGenes = true
// ignore kmers with a count < this (currently not using)
minKmers = 3
kmcFasta = false   // kmc output format is fasta or text (false)
probeFasta = true   // probe file format is fasta (true) or text (false)
outputAll = false
outputProbeName = true
// kmc stops counting after this; assume off-kir hits if this or higher
maxHitCount = 255
// things that probably won't change per run
err = System.err
fileSeparator = System.getProperty('file.separator')
OptionAccessor options = handleArgs(args)
String kmcFileName = options.j
String id = options.i
String extension = options.e
String outputDir = options.o
String probeFileName = options.p

// convert kmc output file to fasta format (.bin1 extension)
kmc2Fasta(kmcFileName, probeFileName, outputDir, id, extension,
		  kmcFasta, probeFasta, ignoreNonGenes)
err.println "done"

/*
 * kmc2Fasta
 *
 * Convert Kmc output to fasta format -- one file per gene containing
 * all the hit probes.
 *
 * @param probesHitsFile String containing the name of the Kmc output file(e.g., '*hits.txt'): tab-separated probe DNA then count (including 0)
 * @todo modularize
 */
int kmc2Fasta(String probeHitsFile, String probeFileName, String outputDir, 
              String id, String extension, Boolean kmcFasta, Boolean probeFasta,
              Boolean ignoreNonGenes) { 
    if(debugging <= 1) { 
        err.println "kmc2Fasta(probeHitsFile=${probeHitsFile}, id=${id}, outputDir=${outputDir}, extension=${extension})"
    }

    int retval = 0
	// map: probe -> Set of loci
    HashMap<String, TreeSet<String>> locusProbeMap = loadProbeMap(probeFileName,
																  probeFasta,
                                                                  ignoreNonGenes) 
	if(debugging <= 3) {
		err.println "kmc2Fasta: ${locusProbeMap.keySet().size()} locus/probes in ${probeFileName} (locusProbeMap)"
	}
    // locus -> probe seq -> count
	// will build this as we loop through the kmc file 
   HashMap<String, HashMap<String, Integer>> locusProbeHitMap = new HashMap()
	// locus -> list of hit counts (including zeros)
	HashMap<String,ArrayList<Integer>> locusHitListMap = new HashMap()
    if(debugging <= 3) { 
        err.println "kmc2Fasta: probeHitsFile=${probeHitsFile}, outputDir=${outputDir}, extension=${extension}"
    }
	
    FileReader kmcReader = new FileReader(new File(probeHitsFile))
    Integer count = 0
	//populate the two data structures with the marker in the line (all ids)
	kmc2FastaLine(kmcReader, locusProbeMap, locusProbeHitMap,
				  locusHitListMap, kmcFasta)

    if(debugging <= 3) { 
        err.println "kmc2Fasta: ${locusProbeHitMap.size()} loci in locusProbeHitMap"
        err.println "kmc2Fasta: " + locusProbeHitMap.keySet().join(",")
    }
	// loop through all the probes and
    // add the zeros (non-hits) back to locusHitListMap
	locusProbeMap.keySet().each { probe ->
		locusSet = locusProbeMap[probe] // set of loci defined for this probe
		zeroCount = 0
		locusSet.each { loc ->
			Map pMap = locusProbeHitMap[loc]
			Integer pcount = 0
			Integer prcCount = 0
			if(pMap != null) {
				pcount = pMap[probe]
				prcCount = pMap[reverseComplement(probe)]
			}
                
			ArrayList plist = locusHitListMap[loc]
			if(plist == null) {
				plist = new ArrayList()
			}
			if(pcount == null) {
				pcount = 0
			}
			if(prcCount == null) {
				prcCount = 0
			}
			if(((pcount == 0) && (prcCount == 0) &&
               locusSet.contains(probe)) || ((pcount < minKmers) &&
                                             (prcCount < minKmers))) {
                if(debugging <= 2) { 
                    err.println "kmc2Fasta: adding zero to $loc $probe"
                }
				plist.add(0)
				zeroCount++
            } else if((pcount >= maxHitCount) ||
                      (prcCount >= maxHitCount) ||
                      !locusSet.contains(probe)) { // off-kir
                if(debugging <= 2) { 
                    err.println "kmc2Fasta: off-kir for $loc probe $probe"
                }
                // remove from locusProbeHitMap
                locusProbeHitMap.remove(probe)
                // remove from locusHitListMap later
                plist.removeAll(maxHitCount) // remove all the off-kir hits
			} else if((pcount < minKmers) ||
                      (prcCount < minKmers)) {
                
			}
		} // each locus defined for this probe
		if((zeroCount > 0) && (debugging <= 1)){ 
			err.println "kmc2Fasta: added $zeroCount zeros for $probe in " +
				locusSet.join(",")
		}
	} 

    // output hits per locus with 'extension' extension
	if(debugging <= 5) {
	    err.println "outputting ${extension} files ..."
    }
	loci = locusProbeHitMap.keySet()
    loci.each { loc ->
		/*if((loc != "2DL3") && (loc != "2DL1")) {
			return //todo(remove)
		}*debugging*/
		//err.println loc //todo
		// replace slash in locus (e.g., 7/24) with 'z'
		escapedLoc = loc.replaceAll("/", "z")
        fullFileName = "${outputDir}${fileSeparator}${id}_${escapedLoc}.${extension}"
        ArrayList plist = locusHitListMap[loc]
        if(plist == null) {
			if(debugging <= 2) {
				err.println "no hit list for ${loc}"
			}
			return
		}

		if(debugging <= 3) {
			err.println "kmc2Fasta: $loc plist=" + plist
		}
		double[] pListD = plist.toArray()
//todo(testing)		double[] avgList = StatUtils.mode(pListD)
        double[] avgList = Math.floor(StatUtils.mean(pListD))
		if(debugging <= 4) {
			err.println "$loc avgList=" + avgList
		}
        if((avgList == null)  || (avgList.size() == 0)){
			return
		}
		Float avg = 0
		if(avgList.size() > 1) {
			// mean the modes just to get one number
			avg = StatUtils.mean(pListD)
			if(debugging <= 2) {            
			    err.println "reducing the mode to one number: mean=${avg}"
            }
		} else {
			avg = avgList[0]
			if(debugging <= 2) {
				err.println "kmc2Fasta: setting to first $loc avgList: " + avgList[0]
			}
		}

		if((avg == 0) && (outputAll != true)){
			return
		}

        outWriter = new PrintWriter(new File(fullFileName).newOutputStream(),
									true)
		// non-null plist from here down
        plist.each { hitCount ->
			if(hitCount == 0) {
				return
			}
            if(debugging <= 2) { 
                err.println "kmc2Fasta: writing ${loc} ${hitCount}"
            }
			// this is in the file name
			if(outputProbeName) { 
				outWriter.println ">${loc}"
			}
            outWriter.println "${hitCount}"
        }
        outWriter.close()
    } // each locus

    if(debugging <= 1) { 
        err.println "kmc2Fasta: return ${retval}"
    }
    return retval
} // kmc2Fasta

/*
 * kmc2FastaLine
 *
 * Process each line of the kmc input file. Populate the two data structures.
 * 
 * @param line a line from the kmc input file
 *   e.g. marker  label   pvalue 100a    100b    100c ...
 * @param locusProbeMap Map: probe -> Set of loci
 * @param locusProbeHitMap locus -> probe seq -> count
 * @param locusHitListMap Map: locus -> list of hit counts (including zeros)
 */
def kmc2FastaLine(FileReader kmcReader, HashMap<String, TreeSet<String>> locusProbeMap,
				  HashMap<String, HashMap<String, Integer>> locusProbeHitMap,
				  HashMap<String,ArrayList<Integer>> locusHitListMap,
				  Boolean kmcFasta) {
	if(debugging <= 1) {
		err.println "kmc2FastaLine()"
	}

	// todo: modularize this
	// fasta format
	// e.g., >3DP1-2DL4
	//        AACATAAGCCAGTAGAATAGCATCT
	String probe = null
	Integer count = 0
    while(line = kmcReader.readLine()) {
		if(debugging <= 1) {
			err.println "kmc2FastaLine: line=$line"
		}
		if(!line || (line == "")) {
			continue
		}
		if(kmcFasta == true) {
			locus = line[1..-1].trim()
			if(debugging <= 1) {
				err.println "kmc2FastaLine: setting locus to $locus"
			}

			line = kmcReader.readLine()  // read the sequence row
			probe = line.trim()
			if(debugging <= 2) { 
				err.println "kmc2FastaLine: $locus $probe"
			}
			count = 10 // default to 10x for fasta input genotypes
		} else {
			// kmc text format
			(probe, countStr) = line.split('\t')
			count = countStr.toInteger()
			if(debugging <= 2) { 
				err.println "kmc2FastaLine: non-fasta: $probe $count"
			}
		}
		// here is where the kmc output is standardized to the
		// complementarity in the probe input file
        // this also subtracts the input reference probes (locusProbeMap)
        // from the genotyped probes (locusProbeHitMap and locusHitListMap)
		probeRc = reverseComplement(probe)
		locusSet = locusProbeMap[probe]
		if(locusSet == null) {		
			// check reverse complement
			locusSet = locusProbeMap[probeRc]
			if(locusSet == null) {
				if(debugging <= 2) { 
					err.println "kmc2FastaLine: no locus set for $probe"
				}
				continue
			} else {
				probe = probeRc
			}
		}
		if(debugging <= 2) { 
			err.println "kmc2FastaLine: $probe: " + locusSet.join(", ")
		}
		locusSet.each { locus ->
			if(debugging <= 1) {
				err.println "kmc2FastaLine: hit ${probe}"
			}

			// list of probes and their hits per locus
			ArrayList hitList = locusHitListMap[locus]
			Map probeHitMap = locusProbeHitMap[locus]
			/* this can get huge
			if(debugging <= 2) {
				err.println "hitList=$hitList"
				(huge) err.println "probeHitMap=" + probeHitMap
			}
  		    */
			if(probeHitMap == null) { 
				probeHitMap = new HashMap()
				hitList = new ArrayList()
			}
			currentCount = probeHitMap[probe]
			if(currentCount != null) {
				count = count + currentCount // forward and reverse complement
			}
			hitList.add(count)
			locusHitListMap[locus] = hitList
			
			probeHitMap[probe] = count
			locusProbeHitMap[locus] = probeHitMap
		} // each locus for a probe
	} // each line
	if(debugging <= 1) {
		err.println "kmc2FastaLine: return"
	}
} // kmc2FastaLine

// http://groovyconsole.appspot.com/script/29005
def String reverseComplement(String seq) { 
    if(debugging) { 
        //err.println "reverseComplement(seq=${seq})"
    }
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]

    seq = seq.toUpperCase().replaceAll( /./ ) { complements."$it" ?: 'X' } 
    seq = seq.reverse()
    if(debugging) { 
        //err.println "reverseComplement: return ${seq}"
    }
    return seq
} // reverseComplement

/*
 * loadProbeMap
 *
 * @param probeListFileName a reduced and sorted version of the exact output matrix: marker  label   pvalue 100a    100b    100c ...
 * @return Map: probe -> Set of loci
 */
HashMap<String,TreeSet<String>> loadProbeMap(String probeListFileName,
											 Boolean probeFasta,
                                             Boolean ignoreNonGenes) {
	if(debugging <= 1) {
		err.println "loadProbeMap(probeListFileName=${probeListFileName}, probeFasta=$probeFasta"
	}
    // open file with probes
    FileReader probeReader = new FileReader(new File(probeListFileName))
	// map: probe -> Set of loci
    HashMap<String,TreeSet<String>> probeMap = new HashMap()
 	String locus = null
	String probe = null
    probeReader.eachLine { line ->
		if(probeFasta == true) {
			start = line.startsWith('>')
			if(start == true) {
				locus = line[1..-1].trim()
			} else {
				probe = line.trim()
				/*if(debugging <= 2) { 
					err.println "fasta: $locus $probe"
				}*/
				addToProbeMap(probeMap, locus, probe, ignoreNonGenes)
				locus = null
				probe = null
			}
		} else {
			// kmc text format
			(locus, probe) = line.split('\t')
			if(debugging <= 2) { 
				err.println "non-fasta: $locus $probe"
			}
			addToProbeMap(probeMap, locus, probe, ignoreNonGenes)
		}
    } // each line of file

    probeReader.close()
	if(debugging <= 3) {
		err.println "loadProbeMap: return ${probeMap.keySet().size()} probe/regions"
	}

	//err.println probeMap //todo
    return probeMap
} // loadProbeMap

void addToProbeMap(HashMap<String,TreeSet<String>> probeMap,
				   String locus, String probe, Boolean ignoreNonGenes) {
    if(ignoreNonGenes && (locus.contains('-') || locus.contains('~'))) {
        return
    }

	s = probeMap[probe]
	if(s == null) {
		s = new TreeSet()
		probeMap[probe] = s
	}
	s.add(locus)

} // addToProbeMap
/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'kmc2LocusAvg2.groovy [options] ',
      header:'Options:')
    cli.help('print this message')

    cli.j(longOpt:'kmc probe results', args:1, argName:'kmc', 'input kmc counts',
      required: true)
    cli.o(longOpt:'directory to put the output', args:1, argName:'out', 
      'output directory', required: true)
    cli.i(longOpt:'ID', args:1, argName:'id', 'ID for the individual',
      required: true)
    cli.e(longOpt:'extension', args:1, argName:'extension', 'extension for output files',
      required: true)
    cli.p(longOpt:'probes', args:1, argName:'probes', 'input probes',
      required: true)

    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs
