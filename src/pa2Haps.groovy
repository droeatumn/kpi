#!/usr/bin/env groovy

/*
 * pa2Haps
 *
 * Fit PA genotypes to (potentially ambiguous) haplotype pairs.
 * best to worst: gene+hap, inter-gene reduced, gene only, hap only
 * Explicitly doubles cA01~tA01 (1).
 * 
 * todo: update example
 * e.g., pa2Haps.groovy -h $HOME/doc/kir/snp/all_haps_v4.txt -q 2DL2.bin1,3DL2.bin1,2DL3.bin1,3DL3.bin1,2DL4.bin1,3DP1-2DL4.bin1,2DP1.bin1,3DP1.bin1,2DS2.bin1,cA01~tA01.bin1,2DS4.bin1,cB02~tA01.bin1,3DL1.bin1,cB02~tB01.bin1 -o prediction.txt
 *
 * 
 * Requires
 *  guava.jar: https://github.com/google/guava
 *    http://google.github.io/guava/releases/19.0/api/docs/com/google/common/collect/Table.html
 *    export CLASSPATH=$HOME/bin/guava/guava/target/guava-27.1-jre.jar:$CLASSPATH
 *
 * @author Dave Roe
 * @todo use numeric sorting for haplotypes in output
 * @todo add absence-only testing for haplotypes
 */

import groovy.io.*
import org.apache.commons.cli.CliBuilder
import org.apache.commons.cli.OptionAccessor
// newer?
//import groovy.cli.commons.CliBuilder
//import groovy.cli.commons.OptionAccessor
import com.google.common.collect.Table
import com.google.common.collect.HashBasedTable

// things that may change per run
debugging = 3 // TRACE=1, DEBUG=2, INFO=3
// for all_haps_v6.xlsx, loci start at column index 5
// haplotype	nomenclature	freq	structure	3DL3
Integer startLocusIndex = 4

// thing that probably won't change per run
err = System.err
ArrayList<String> loci = ["3DL3", "2DS2", "2DL2", "2DL3", "2DP1", "2DL1", "3DP1", "2DL4", "3DL1", "3DS1", "2DL5", "2DS3", "2DS5", "2DS4", "2DS1", "3DL2"]

// can call absent, but no present
HashSet<String> absentOnlyRefLociSet = new HashSet()
// these two haplotypes have present/absent markers
HashSet<String> refHapLociSet = new HashSet(Arrays.asList())
HashSet<String> absentOnlyRefHapLociSet = new HashSet()
//HashSet<String> absentOnlyRefHapLociSet = new HashSet(Arrays.asList('cB02~tB01', 'cA01~tB01', 'cB01~tA01', 'cB01~tB01', 'cB02~tA01', 'cA01~tA01', 'cA01~tA02'))
// not relevant for locus determination
//HashSet<String> alleleSet = new HashSet(Arrays.asList('2DS3', '2DS5'))
OptionAccessor options = handleArgs(args)

if(debugging <= 3) {
    err.println "pa2Haps -h ${options.h} -q ${options.q} -o -${options.o}"
}
// open file with haplotype definitions
FileReader hapReader = new FileReader(new File(options.h))
// string with probe query hits
// e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
String qString = new String(options.q)

// load input files into two tables, one for present/absent genotypes
// and one for absent-only genotypes
// Table<hap number, locus names, gene count>
HashBasedTable<String, String, Boolean> absHapTable
HashBasedTable<String, String, Boolean> paHapTable
HashMap<String, Float> freqMap // hap number -> frequency
HashMap<String, String> nomenclatureMap // hap name -> hap number
(absHapTable, paHapTable, freqMap, nomenclatureMap) =
    readReferenceHaplotypes(absentOnlyRefLociSet, hapReader,
                            startLocusIndex)
if(debugging <= 1) {
    err.println "nomenclatureMap=" + nomenclatureMap
}
// loci that can be used for present and absent testing
// the companion to absentOnlyRefLociSet
paRefLociSet = paHapTable.columnKeySet()

// Table<hap pair names, locus names, gene count>
HashBasedTable<String, String, Boolean> paHapPairTable =
    makeHapPairTable(paHapTable)
HashBasedTable<String, String, Boolean> absHapPairTable =
    makeHapPairTable(absHapTable)
log1(paHapTable, absHapPairTable, debugging)

HashSet<String> hitSet = parseGenotypes(qString)

// _all_ P/A gene genotype from genes -> Boolean
HashMap<String,Boolean> genPAMap = new HashMap()
// genHitSet: a set of genes _that hit_
// it is like hitSet, except there may be some nomenclature changes
HashSet<String> genHitSet = new HashSet()
makeGenotypeMaps(hitSet, paHapTable.columnKeySet(),
                 paHapPairTable, nomenclatureMap, genPAMap, genHitSet)
// first, match genotype with pairs of PA haplotypes
// i.e, without the absent only
// haplotype pairs -> Map[locus: boolean]
// make hap pair predictions from gene probe hits
HashMap<String, HashMap<String,Boolean>> interpPAMap = null
interpPAMap = interpret(paHapPairTable, genPAMap, new Boolean(false))
// now narrow the results by the absent-only loci
interpAbsMap = interpret(absHapPairTable, genPAMap, new Boolean(true))
HashSet<String> nonHapPredictionSet = interpPAMap.keySet()
if((interpAbsMap != null) && (interpAbsMap.size() > 0)) {
    nonHapPredictionSet = interpPAMap.keySet().intersect(interpAbsMap.keySet()) // ~
}
if(debugging <= 3) { 
    err.println "abs" + interpAbsMap.keySet()
    err.println "first ${interpPAMap.keySet().size()} PA marker interpretation(s) " +
        interpPAMap.keySet().join('|')
    err.println "narrowed with absent-only markers:  ${nonHapPredictionSet.size()} PA marker interpretation(s) " +
        nonHapPredictionSet.join('|')
}

// make hap pair predictions from haplotype probe hits
HashSet<String> interpHapSet = null
interpHapSet = interpretHapMarkers(genHitSet, nomenclatureMap.values().toSet())
if(debugging <= 3) { 
    err.println "${interpHapSet.size()} haplotype marker interpretation(s)"
}

writeOutput(options, genPAMap, genHitSet, nonHapPredictionSet,
            interpHapSet, refHapLociSet, absHapPairTable, absentOnlyRefLociSet,
            absentOnlyRefHapLociSet, paHapPairTable, loci)

// end main

/*
 * writeOutput
 *
 * Make homozygous if the haplotype is cA01~tA01 (1).
 * 
 * Send the haplotype pair predictions to the output file.
 * @param options OptionAccessor contains option 'o' with String with output file name
 * @param genPAMap HashMap<String,Boolean> of input of _all_ loci genotype hit map
 * @param genHitSet HashSet<String> of loci that _hit_ (present)
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @param paHapPairTable Table<hap number, locus names, gene count>
 * @param refHapLociSet HashSet<String> of the reference haplotypes that are tested
 * @param absHapPairTable Table of hap pair names and loci to gene counts
 *
 */
def void writeOutput(OptionAccessor options, Map genPAMap,
                     Set genHitSet, 
                     HashSet<String> interpPASet,
                     HashSet<String> interpHapSet, HashSet<String> refHapLociSet,
                     HashBasedTable<String, String, Boolean> absHapPairTable,
                     HashSet<String> absentOnlyRefLociSet,
                     HashSet<String> absentOnlyRefHapLociSet,
                     HashBasedTable<String, String, Boolean> paHapPairTable,
                     ArrayList<String> loci) {
    if(debugging <= 1) {
        err.println "writeOutput(interpPASet=${interpPASet}, interpHapSet=${interpHapSet})"
    }
    String outFileName = options.o
    // open output file
    outWriter = new PrintWriter(new File(outFileName).newOutputStream(), true)
    // output header
    header = "id\thaplotypes"
    loci.each { locus ->
        header += "\t${locus}"
    }
    outWriter.println header

    startIdx = outFileName.indexOf(System.getProperty('file.separator')) + 1
    if(startIdx < 0) {
        startIdx = 0
    }
    // trim things after the last "_"
    //endIdx = outFileName.lastIndexOf("_") - 1
    //if(endIdx < 0) {
        endIdx = outFileName.length() - 1
    //}
    id = outFileName[startIdx..endIdx]
    //id = outFileName.replaceFirst("_prediction.txt", "")
	
    // output genotypes
    TreeSet<String> intergeneHitSet = new TreeSet()

    // todo: move this before writeOutput?
    HashSet<String> pacombinedSet = paReduceByHap(genHitSet,
												  interpPASet,
												  interpHapSet, refHapLociSet,
                                                  absHapPairTable,
                                                  absentOnlyRefLociSet,
                                                  absentOnlyRefHapLociSet)
	if(debugging <= 3) {
		err.println "haplotype-reduced gene interp from ${interpPASet.size()} to ${pacombinedSet.size()}"
	}    
	if(debugging <= 2) {
        err.println "${pacombinedSet.size()} combined predictions"
		err.println "genotype: " + genHitSet.sort().join("+")
        err.println "gene: " + interpPASet.sort().join('|')
        err.println "haplotype: " + interpHapSet.sort().join('|')
        err.println "pa combined: " + pacombinedSet.sort().join('|')
        err.println options.a
    }
    if(pacombinedSet.size() == 0) {
        outStr = "${id}\tuninterpretable"
        loci.each { locus ->
            genePA = "N"
            if(genHitSet.contains(locus)) { 
                genePA = "Y"
            }
            outStr += "\t${genePA}"
        }
    } else {
        hapPair = pacombinedSet[0]
        locusMap = paHapPairTable.row(hapPair)
        //todo print out locusMap?
        // this is isn't right; absent-only could be present or absent
        absLocusMap = absHapPairTable.row(hapPair) 
        outStr = "${id}\t${pacombinedSet.sort().sort().join('|')}"
        if(debugging <= 2) {
            err.println "writeOutput: hapPair=$hapPair"
            err.println "writeOutput: absent-only $hapPair = " + absHapPairTable.row(hapPair)
            err.println "writeOutput: absLocusMap = " + absLocusMap
        }
        locusMap = paHapPairTable.row(hapPair)
        loci.each { locus ->
            genePA = locusMap[locus] ?  "Y" : "N"
            outStr += "\t${genePA}"
        }
    } // interp or uninterp
    
    outWriter.println outStr
	outWriter.close()
	if(debugging <= 1) {
        err.println "writeOutput: return"
    }
} // writeOutput

/* 
 * reduceGeneByInterGene
 *
 * @param geneHitSet full set of gene hits
 * @param genHitSet full set of intergene hits
 * @param hapTable Table<hap number, locus names, gene count>
 */
HashSet<String> reduceGeneByInterGene(HashSet<String> interpPASet,
									  ArrayList<String> intergeneHitSet,
									  HashBasedTable<String, String, Boolean> paHapPairTable) {
 	if(debugging <= 1) {
		err.println "reduceGeneByInterGene(${interpPASet.size()} haplotype pairs)"
	}
	HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0

    Iterator e = interpPASet.iterator()
	Map<String,String> hapPairRowMap = paHapPairTable.rowMap()

	// loop through the PA hap pairs (e.g., 1+4)
	// and check if all the expected hit are found in the genotype
	// note, the opposite is not checked: all the genotype hits
	// are found in the hap pairs
    while ( e.hasNext() ) {
        String paPair = (String)e.next().toString();
		Integer paPairCount = 0
 		if(debugging <= 2) {
			err.println "reduceGeneByInterGene: paPair=$paPair"
		}
		// get the intergene regions for this haplotype
		hapPairColMap = hapPairRowMap[paPair]

		// count how many are really in the genotype
		hapPairColMap.each { geneOrInter, present ->
			if(debugging <= 2) {
				err.println "reduceGeneByInterGene: geneOrInter=$geneOrInter"
			}
			// if intergene present in hap pair
			if(geneOrInter.contains("-") && (present)) {
				if(debugging <= 2) {
					err.println "reduceGeneByInterGene: $paPair should have $geneOrInter ..."
				}
				if(intergeneHitSet.contains(geneOrInter)) {
					paPairCount++
					if(debugging <= 2) {
						err.println "reduceGeneByInterGene: and it does: total count now: $paPairCount"
					}
				}
			}
		} // each gene/intergene for this region 
		
		if(paPairCount > highestHitcount) { 
            if(debugging <= 2) {
                err.println "reduceGeneByInterGene: new highest: $paPair $paPairCount"
            }
            outSet = new HashSet()
            outSet.add(paPair)
            highestHitcount = paPairCount
        } else if(paPairCount == highestHitcount) {
            if(debugging <= 2) {
                err.println "reduceGeneByInterGene: adding $paPair to $paPairCount"
            }
            outSet.add(paPair)
		}
	} // each paPair

 	if(debugging <= 1) {
		err.println "reduceGeneByInterGene: return ${outSet.size()} haplotype pairs"
	}
    return outSet
} // reduceGeneByInterGene

/* 
 * Combine the two (gene and haplotype) interpretations
 * by going with one if the other doesn't exist, or reducing
 * the PA haplotype list by haplotype markers if possible.
 * 
 * @param genHitSet Set of loci (genes and haps) that hit; *this is altered*
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation (~)
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @param refHapLociSet HashSet<String> of the reference haplotypes that are tested
 * @param absHapPairTable Table of hap pair names and loci to gene counts
 * @param absentOnlyRefLociSet Set of genes that have absent-only markers
 * @return Set of the new haplotype pair predictions; one pair per String
 */
HashSet<String> paReduceByHap(Set genHitSet,
							  HashSet<String> interpPASet,
							  HashSet<String> interpHapSet,
                              HashSet<String> refHapLociSet,
                              HashBasedTable<String, String, Boolean> absHapPairTable,
                              HashSet<String> absentOnlyRefLociSet,
                              HashSet<String> absentOnlyRefHapLociSet) {
	if(debugging <= 1) {
        err.println "paReduceByHap(interpPASet=${interpPASet.join(',')}, interpHapSet=${interpHapSet.join(',')})"
    }

    HashSet<String> outSet = new HashSet()
    Integer highestScore = 0
    Iterator e = interpPASet.iterator()
    // loop through each PA haplotype pair prediction
    while ( e.hasNext() ) {
        // all single haplotypes in e
        TreeSet interpPASingleHapsSet = new TreeSet()
        String paPair = (String)e.next().toString();
	    if(debugging <= 2) {
            err.println "paReduceByHap: evaluating pa pair ${paPair}"
        }
        // for this PA hap pair prediction, score how it fits the ref haps
        score = 0
        // todo: loop through absentOnlyRefHapLociSet and eliminate that way?
        paPair.split("\\+").each { paHap ->
            interpPASingleHapsSet.add(paHap)
            //paHapNoTilde = paHap.replaceFirst('~', '-')
            if(debugging <= 2) { 
                //err.println paHapNoTilde
                err.println genHitSet
                err.println genHitSet.contains(paHap)
                //err.println interpHapSet
                //err.println interpHapSet.contains(paHap)
            }

            // if the pa hap call is also in the haplotype call
            if(genHitSet.contains(paHap)) {
                if(refHapLociSet.contains(paHap)) {
                    // if the paHap is in the PA haplotype hit list, score +2
                    score += 2
                } // else neutral score
            } else { // pa call is not in the haplotype call
                if(refHapLociSet.contains(paHap)) {
                    // if the paHap is in the PA haplotype hit list, score -2
                    score -= 2
                } else { 
                    // if the paHap is in the absent-only list, score -2
                    if(absentOnlyRefHapLociSet.contains(paHap)) {
                        score -= 2
                    }
                } // else leave the score neutral
            }
        } // each haplotype in a pair
        if(debugging <= 2) {
            err.println "first count: ${score}"
        }
        // bump score for this predicted haplotype pair if the markable haplotypes
        // don't hit in both the predicted pa and haplotype calls
        if(debugging <= 2) {
            err.println "paReduceByHap: score ${score} for ${paPair}"
        }
        if(score > highestScore) {
            if(debugging <= 2) {
                err.println "paReduceByHap: setting highestScore"
            }
            outSet = new HashSet()
            outSet.add(paPair)
            highestScore = score
        } else if(score == highestScore) {
            if(debugging <= 2) {
                err.println "paReduceByHap: adding highestScore"
            }
            outSet.add(paPair)
        }
    }
	
    if(debugging <= 1) {
        err.println "paReduceByHap: genHitSet=${genHitSet}"
        err.println "paReduceByHap: returning ${outSet}"
    }
    return outSet
} // paReduceByHap

/* 
 * Combine the two (gene and haplotype) interpretations
 * by going with one if the other doesn't exist, or reducing
 * the haplotype marker list if possible using the PA interpreted pairs.
 * It only checks pairs right now.
 * 
 * @param genHitSet Set of loci (genes and haps) that hit
 * @param interpPASet Set of Strings with haplotype predictions from gene interpretation
 * @param interpHapSet Set of Strings with haplotype predictions from haplotype interpretation
 * @return Set of the new haplotype pair predictions; one pair per String
 * @todo split the pairs and check each haplotype
 */
HashSet<String> hapReduceByPA(Set genHitSet,
							  HashSet<String> interpPASet,
							  HashSet<String> interpHapSet) {
    HashSet<String> outSet = new HashSet()
    Integer highestHitcount = 0
    Iterator e = interpHapSet.iterator()
    // loop through each hap-marker haplotype pair prediction
    while ( e.hasNext() ) {
        String hapPair = (String)e.next().toString();
        // for this hap-marker hap pair prediction, count how many of
        // those haplotypes typed as present in pa data
        hitcount = 0 
        if(interpPASet.contains(hapPair)) {
            hitcount++
        }
		// automatically rule in any pair of haplotypes
		// where at least one is a haplotype without markers
        hapPair.split("\\+").each { hap ->
		    if(!hapWithMarkers(hap)) {
			    outSet.add(hapPair)
		    }
        }
        if(debugging <= 2) {
            err.println "hapReduceByPA: hitcount ${hitcount} for ${hapPair}"
        }
        if(hitcount > highestHitcount) {
            if(debugging <= 2) {
                err.println "hapReduceByPA: setting highestHitcount"
            }
            outSet = new HashSet()
            outSet.add(hapPair)
            highestHitcount = hitcount
        } else if(hitcount == highestHitcount) {
            if(debugging <= 2) {
                err.println "hapReduceByPA: adding highestHitcount"
            }
            outSet.add(hapPair)
        }
    } // while
	
    return outSet
} // hapReduceByPA

/*
 * Returns true if the input haplotype name/number is 
 * one for which we have markers.
 */
Boolean hapWithMarkers(String hap) {
	Boolean ret = false
	hapsWithMarkers = ['1', '3/11', '4', '6/25', '7/9', '8/17', '32', '98/99', 'cA01~tA01', 'cA01~tB01', 'cB02~tA01', 'cB01~tA01', 'cB01~tB01', 'cB02~tB01', 'cB04~tB03', 'cB01~tB05', 'cA01~tB04', 'cB05~tB01', 'cA01~tB05', 'cB05~tA01', 'cA01~tB06', 'cA01~tA02', 'cA03~tB02', 'cB03~tA01']
	if(hapsWithMarkers.contains(hap)) {
		ret = true
	}
	return ret
} // hapWithMarkers

/*
 * makeHapPairTable
 *
 * @param hapTable single haplotype table of <hap, gene, Boolean>
 * @return Table all haplotype-pairs of <hap pair, gene, Boolean>
 */
HashBasedTable<String, String, Boolean> makeHapPairTable(Table hapTable) {
    if(debugging <= 1) { 
        err.println "makeHapPairTable()"
    }
    HashBasedTable<String, String, Boolean> paHapPairTable = HashBasedTable.create()
    HashSet<String> hapList = hapTable.rowKeySet()
    if(debugging <= 2) { 
        err.println "hapList=${hapList}"
    }
    Set<String> locusList = hapTable.columnKeySet()

//    for(int i=0; hap1 = hapList[i++]; i < hapList.size()) {
    for(hap1 in hapList) {
//        for(int j=0; hap2 = hapList[j++]; j < hapList.size()) {
        for(hap2 in hapList) {
            if(hap1.compareTo(hap2) > 0) { // sort and eliminate duplicates
                continue
            }
            String hapPair = "${hap1}+${hap2}"
            locusList.each { locus ->
                Boolean cellVal = false
                if(hapTable.get(hap1, locus) || hapTable.get(hap2, locus)) {
                    cellVal = true
                }
                paHapPairTable.put(hapPair, locus, cellVal)
            } // each gene
        }
    }
    if(debugging <= 1) { 
        err.println "makeHapPairTable: return"
    }
    return paHapPairTable
} // makeHapPairTable

/*
 * readReferenceHaplotypes
 *
 * Read the haplotype definitons and return their details in
 * two tables: one for absent-only loci and one for both.
 *
 * haplotype	nomenclature	freq	structure	3DL3	2DS2	2DL2	2DL3	2DP1	2DL1	3DP1	2DL4	3DL1	3DS1	2DL5	2DS3	2DS5	2DS4	2DS1	3DL2
 * 1	cA01~tA01	55.40%		1	0	0	1	1	1	1	1	1	0	0	0	0	1	0	1
 *
 * @param absentOnlyRefLociSet Set of loci w/probes that can only determine absence.
 * @param hapFile tab-separated file containing the reference haplotypes, their nomenclature name and frequencies
 * @param startLocusIndex index of first allele (3DL3) in spreadsheet
 * @return HashBasedTable containing the haplotype definitions for the absent-only loci
 *    rows headers are the haplotypes
 *    column headers are the loci
 *    cell values are the locus count
 * @return HashBasedTable containing the haplotype definitions for the present/absent loci
 *    rows headers are the haplotypes
 *    column headers are the loci
 *    cell values are the locus count
 *  HashMap<haplotype number, Float> for the frequencies
 *  HashMap<haplotype number, String> for the nomenclature names
 */
def List readReferenceHaplotypes(Set absentOnlyRefLociSet,
                                 FileReader reader,
                                 Integer startLocusIndex) {
    ArrayList retList = new ArrayList()
    // hap number -> frequency
    HashMap<String, Float> freqMap = [:]
    // hap number -> hap name
    HashMap<String, String> nomenclatureMap = [:]

    // one table each for presence/absence and absent-only
    // row haplotypes, column loci, cell locus counts
    HashBasedTable<String, String, Boolean> paHapTable = HashBasedTable.create()
    HashBasedTable<String, String, Boolean> absHapTable = HashBasedTable.create()
    // first row are loci
    TreeSet<Integer> loci = new TreeSet()
    ArrayList header = reader.readLine().split('\t')
    header[startLocusIndex..-1].each { locus ->
        if(debugging <= 2) {
            err.println "readReferenceHaplotypes: adding ${locus} to loci list"
        }
        loci.add(locus)
    }

    // non-header rows
    reader.eachLine { line ->
        if(debugging <= 1) {
            err.println "readReferenceHaplotypes: line=$line"
        }
        ArrayList cols = line.split('\t')
        String haplotype = cols[1].trim() // haplotype name
        if((haplotype == null) || (haplotype == "")) { // skip blank rows
            return
        }
        freqMap[haplotype] = cols[2].toFloat()
        //old nomenclatureMap[haplotype] = cols[1]
        nomenclatureMap[cols[1]] = haplotype // new
        // the loci are from here to the end
        cols[startLocusIndex..-1].eachWithIndex { cell, i ->
            i+=+startLocusIndex
            Boolean genePresent = false
            if(cell.toInteger() > 0) {
                genePresent = true
            }
            locus = header[i]
            // skip columns wo headers and  haplotypes
            if(locus == null) {
                return
            }
            if(debugging < 2) {
                err.println "readReferenceHaplotypes: adding ${locus}=${genePresent} for ${haplotype}"
            }
            if(absentOnlyRefLociSet.contains(locus)) {
                absHapTable.put(haplotype, locus, genePresent)
                if(debugging < 2) {
                    err.println "readReferenceHaplotypes: added absent only ${haplotype},${locus}=${genePresent}"
                }
            } else { 
                paHapTable.put(haplotype, locus, genePresent)
                if(debugging < 2) {
                    err.println "readReferenceHaplotypes: added pa ${haplotype},${locus}=${genePresent}"
                }
            }
        }
    } // each line/haplotype

    if(debugging < 1) {
        err.println "readReferenceHaplotypes: return ${paHapTable.rowKeySet().size()} rows, ${paHapTable.columnKeySet().size()} columns"
    }
    retList.add(absHapTable)
    retList.add(paHapTable)
    retList.add(freqMap)
    retList.add(nomenclatureMap)
    return retList
} // readReferenceHaplotypes

/*
 * readProbes
 *
 * @param reader tab-separated file containing the loci and their probes; no column headers; locus in first column, probe in second
 * @return a Map of probe to gene/locus
 */
def HashMap<String,String> readProbes(FileReader reader) {
    HashMap<String,String> probeMap = new HashMap()
    reader.eachLine { line ->
        (gene, probe) = line.split('\t')
        if((gene == null) || (gene == "")) {
            return
        }
        probeMap[probe] = gene
    } // each line

    return probeMap
} // readProbes

/*
 * parseGenotypes
 *
 * Returns a Set of the locus names based of the file names
 * with the '.bin1' extension.
 *
 * e.g., 2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3.bin1,2DL3-2DL5B_2DL3-2DP1.bin1,2DL3-2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1.bin1,2DP1-2DL1_2DP1-2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1.bin1,2DS1-3DL2.bin1,3DL2.bin1,3DL2.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3.bin1,3DL3-2DL3.bin1,1G.bin1,2G.bin1
 * @param hitString String containing loci that hit
 * @return Set of only the probes that hit
 */
def HashSet<String> parseGenotypes(String reader) {
    if(debugging <= 1) {
        err.println "parseGenotypes(reader=${reader})"
    }
    
    HashSet<String> probeSet = new HashSet()
    reader.split(',').each { multiLocus ->
        multiLocus.split('_').each { locus ->
            lIn = locus.replaceFirst("\\.bin1", "") //todo pass this in
			/*if(lIn =~ /2DL5/) {
				// convert 2DL5A and 2DL5B to 2DL5
				lIn = "2DL5"
			}remove(todo)*/
            if((lIn != "") && !(lIn =~ /unmapped/)) {
                probeSet.add(lIn)
            }
        } // each part of multi locus
    } // each comma-separated list of loci

    if(debugging <= 1) {
        err.println "parseGenotypes: return ${probeSet.size()} loci: ${probeSet}"
    }
    return probeSet
} // parseGenotypes

/*
 * makeGenotypeMaps
 * Populate paProbeMap and allProbeMap
 *
 * @param lociSet Set containing Strings from the '.bin1' loci list
 * @param geneSet Set containing all reference loci 
 * @param hapPair Table<hap pair names, locus names, gene count>
 * @return List of two objects HashMap<String,Boolean> gene PA and HashSet<String> gene all
 */
List makeGenotypeMaps(HashSet<String> lociSet, Set<String> geneSet,
                      HashBasedTable<String, String, Boolean> paHapPairTable,
                      HashMap<String, String> nomenclatureMap,
                      HashMap<String,Boolean> paProbeMap,
                      HashSet<String> allProbeMap) {
    if(debugging <= 1) {
        err.println "makeGenotypeMaps()"
    }

    // initialze paProbeMap to false for all genes
    for(String gene : geneSet) {
        paProbeMap[gene] = new Boolean(false)
    } // initialze paProbeMap

    lociSet.each { locus ->
        paProbeMap[locus] = new Boolean(true)
        nomLocus = nomenclatureMap[locus]
        if(debugging < 1) { 
            err.println "nomLocus=${nomLocus}, locus=${locus}"
        }
        if(nomLocus != null) {
            locus = nomLocus
        }
        allProbeMap.add(locus)
    } // each locus probe

    if(debugging <= 1) {
        err.println "makeGenotypeMaps: return ${paProbeMap}, ${allProbeMap}"
    }
    return [paProbeMap, allProbeMap]
} // makeGenotypeMaps

/*
 * Given a gene, return it.
 * Given an intergene,
 *   if it doesn't contain a framework gene, return it
 *   if it does contain a framework gene, return it 
 *     and the framework gene
 *
 */
ArrayList<String> checkInterGeneForFramework(String locus) {
	ArrayList<String> ret = new ArrayList()
	ret.add(locus)
	if(locus.contains("-")) {
		(gene1, gene2) = locus.split('-')
		if(gene1.contains("3DL3") || gene1.contains("3DP1") ||
		   gene1.contains("2DL4") || gene1.contains("3DL2")) {
			ret.add(gene1)
		}
		if(gene2.contains("3DL3") || gene2.contains("3DP1") ||
		   gene2.contains("2DL4") || gene2.contains("3DL2")) {
			ret.add(gene2)
		}
	} // intergene
	return ret
} // checkInterGeneForFramework

/*
 * hapPairs2Genotypes
 *
 * Generate synthetic genotypes from pairs of reference haplotypes.
 * 49 choose 2 is 1176.
 * Due to the properties of the Groovy Map '+' operation, a locus from one
 * haplotype overwrites the other, so false values must be null in the Maps.
 *
 * @param hapMap Map of haplotype names to a Map of loci and their Boolean values
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values
 */
def HashMap<String, HashMap<String,Boolean>> hapPairs2Genotypes(HashMap<String,HashMap<String,Boolean>> hapMap, HashMap<String,String> nomMap) {
    estimate = hapMap.size() * 2
    HashMap<String, HashMap<String,Boolean>> genMap = new HashMap(estimate)
    hapMap.each { h1Name, h1Map ->
        hapMap.each { h2Name, h2Map ->
            if(h1Name > h2Name) { // avoid dups and sort by name
                return
            }
            (h1Name, h2Name) = [h1Name, h2Name].sort()
            genMap["${h1Name}+${h2Name}"] = h1Map + h2Map //todo quotes here?
        }
    }
    return genMap
} // hapPairs2Genotypes (todo: remove)

/*
 * interpret
 *
 * Interpret a genotype in the context of haplotype pairs.
 * It assumes that the genotypes just include the genes that are present.
 *
 * @param hapPairTable Table<hap pair names, locus names, gene count>
 * @param genPAMap HashMap _all_ P/A gene genotype from genes -> Boolean
 * @param absentOnly Boolean if genPAMap contains absent-only markers
 * @return a Map from haplotype pairs to a Map of loci and their Boolean values; this is a potentially ambiiguous P/A interpretPAation
 * todo make a more elegant (e.g., vector-based) approach
 */
HashMap<String, HashMap<String,Boolean>> interpret(HashBasedTable<String, String, Boolean> hapPairTable,
                  HashMap<String,Boolean> genPAMap, Boolean absentOnly) {
    if(debugging <= 1) {
        err.println "interpret()"
        err.println "genPAMap = " + genPAMap
    }
    // haplotype pairs -> Map[locus] -> Boolean
    HashMap<String, HashMap<String,Boolean>> retMap = new HashMap()
    // hap name -> Map[locus] -> Boolean
    Map<String, Map<String, Boolean>> rowHapMap = hapPairTable.rowMap()
    // 'haplotype' is really a haplotype pair
    rowHapMap.each { haplotype, hapColMap ->
        if(debugging <= 1) {
            err.println "interpret: comparing with haplotype ${haplotype}"
			//todo err.println "interpret: ${haplotype} ${hapColMap}"
        }

        Boolean allFound = true
        hapColMap.each { locus, value ->
            if(allFound == false) { // hap pair is false if any locus is false
                return
            }
			// don't intepret haplotypes
			if(locus.contains("~") ||
			   locus.startsWith("c") || locus.startsWith("t")) {
				return
			}

            genVal = genPAMap[locus]
            same = false
            // test the haplotype value with the genotype value
            if(((value == true) && (genVal == true)) ||
               (value == false) && ((genVal == null) || (genVal == false))) {
                    same = true
            } else {
                // absent-only means you can trust genotype if absent, but not present
                // so accept the case when the genotype is present but not in the haplotype
                if(absentOnly && (genVal == true) && (value == false)) {
                    same = true
                } else { 
				    if(debugging <= 1) {
					    err.println "interpret: $locus not same: $value (hap) and $genVal (genotype)"
				    }
                }
			}
            allFound = same
        }

        if(allFound) { 
            if(debugging <= 1) {
                err.println "interpret: true"
            }
            retMap[haplotype] = hapColMap
        } else {
            if(debugging <= 1) {
                err.println "interpret: false"
            }
        }
    } // each row/haplotype

    if(debugging <= 1) {
        err.println "interpret: return ${retMap}"
    }
    return retMap
} // interpret


/*
 * interpretHapMarkers
 *
 * Interpret pairs of haplotypes from the hap (not PA) markers.
 * For now, it just makes all combinations.
 *
 * @param g HashSet of all markers _that hit_
 * @param refHaps Set of all names of the reference haplotypes from the ref hap input file
 * @return Set of haplotype pair predictions (e.g., cA01~tA01+cB01~tA01, ...)
 */
HashSet<String> interpretHapMarkers(HashSet<String> g,
									Set<String> refHaps) {
    if(debugging <= 1) {
        err.println "interpretHapMarkers(g=${g.join(",")}, refHaps=${refHaps.join(",")}"
    }
    HashSet<String> retSet = new HashSet()

    // put all genotype hit into retSet
	g.each { gHit ->
        gHit = gHit.replaceFirst('-', '~')
		if(refHaps.contains(gHit)) {
			retSet.add(gHit)
		}
	}

    HashSet<String> retSetPairs = new HashSet()
	if(retSet.size() == 1) {
		hap = retSet.toList()[0]
		// don't homozygousify (for now)
        //retSetPairs.add([hap, hap].join('+'))
		retSetPairs.add(hap)
    } else { 
		retSet.each { hap1 ->
			retSet.each { hap2 ->
				if(hap2 <= hap1) {
					return
				}
				hl = [hap1, hap2].sort().join("+")
				retSetPairs.add(hl)
			}
		}
	}

    if(debugging <= 1) {
        err.println "interpretHapMarkers return: ${retSetPairs}"
    }
    return retSetPairs
} // interpretHapMarkers

def log1(Table paHapPairTable, Table absHapPairTable,
         Integer debugging) {
    if(debugging <= 3) {
        err.print "loaded ${paHapPairTable.columnKeySet().size()} pa loci, "
        err.println "${paHapPairTable.rowKeySet().size()} reference pa haplotypes, "
        err.println "${absHapPairTable.rowKeySet().size()} reference abs haplotypes, "
        err.println "${paHapPairTable.rowKeySet().size()} reference pa haplotype pairs"
        err.println "${absHapPairTable.rowKeySet().size()} reference abs haplotype pairs"
        err.println "pa haplotypes loci " + paHapPairTable.columnKeySet()
    }
} // log1


/*
 * handleArgs
 * 
 * Parse and return the input.
 *
 * @todo e.g., here
 *
 * @param args List of Strings containing command line arguments.
 * @return Option
 */
OptionAccessor handleArgs(String[] args) { 
    CliBuilder cli = new CliBuilder(usage:'pa2Haps.groovy [options] ', header:'Options:')
    cli.help('print this message')
    cli.h(longOpt:'haps', args:1, argName:'file', 'file with haplotype definitions',
      required: true)
    cli.a(longOpt:'all', args:1, argName:'boolean', 'all output if set (debugging)',
      required: false)
    cli.q(longOpt:'qout', args:1, argName:'file', 'file with results of probe queries',
      required: true)
    cli.o(longOpt:'output', args:1, argName:'file', 'output file containing haplotype pair predictions',
      required: true)
    OptionAccessor options = cli.parse(args)
    return options
} // handleArgs

