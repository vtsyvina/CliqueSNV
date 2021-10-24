package edu.gsu.algorithm;

import edu.gsu.algorithm.em.IlluminaEM;
import edu.gsu.algorithm.util.AlgorithmUtils;
import edu.gsu.algorithm.util.CommonReadsIlluminaParallelTask;
import edu.gsu.algorithm.util.DeletionEdgesParallelTask;
import edu.gsu.algorithm.util.FindIlluminaEdgesParallelTask;
import edu.gsu.algorithm.util.TrueFrequencyEstimator;
import edu.gsu.model.Clique;
import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.Interval;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.SNVStructure;
import edu.gsu.start.Start;
import edu.gsu.util.Utils;
import edu.gsu.util.builders.SNVStructureBuilder;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SNVIlluminaMethod extends AbstractSNV {
    public static final int minorCount = al.length() - 2;
    public static final int MINIMUM_COVERAGE_FOR_HAPLOTYPE_SNP = 5;
    public static final int MINIMUM_COVERAGE_FOR_CONSENSUS_SNP = 3 * MINIMUM_COVERAGE_FOR_HAPLOTYPE_SNP;
    public static final int MAX_HAPLOTYPES_BEFORE_EM = 100;
    public int MAX_EDGES_NUMBER;
    private int[][] commonReads;
    private boolean maxEdgesLimitReached = false;
    // we assume them to be sorted
    private List<FindIlluminaEdgesParallelTask.EdgeSummary> noLimitEdges;

    private final Map<String, Integer> knownHits = new HashMap<>();

    // if the consensus was removed once in one window we don't want to try to add it again and again
    private boolean consensusWasRemoved = false;

    @Override
    protected Technology technology() {
        return Technology.ILLUMINA;
    }

    public final IlluminaSNVSample sample;
    private SNVStructure struct;
    private double[][] profile;
    private String consensus;

    //private Map<String, CorrelationContainer> correlationMap;

    /**
     * @param sample Given input reads
     * @param log    if lof is enables
     */
    public SNVIlluminaMethod(IlluminaSNVSample sample, boolean log) {
        this.log = log;
        this.sample = sample;
        this.MAX_EDGES_NUMBER = Start.settings.get("-el") == null ? 17000 : Integer.parseInt(Start.settings.get("-el"));
        this.MAX_READ_ERROR = Start.settings.get("-re") == null ? 0.2 : Double.parseDouble(Start.settings.get("-re"));
        initParameters(100, Double.parseDouble(Start.settings.getOrDefault("-tf", "0.05")));
    }

    /**
     * @param sample       Given input reads
     * @param minThreshold threshold for O22 value to consider alleles as SNPs
     * @param minFreq      minimum frequency (relative to reads' coverage) for O22 value to consider alleles as SNP
     * @param log          if log is enabled
     */
    public SNVIlluminaMethod(IlluminaSNVSample sample, int minThreshold, double minFreq, boolean log) {
        this.log = log;
        this.sample = sample;
        this.MAX_EDGES_NUMBER = Start.settings.get("-el") == null ? 17000 : Integer.parseInt(Start.settings.get("-el"));
        this.MAX_READ_ERROR = Start.settings.get("-re") == null ? 0.2 : Double.parseDouble(Start.settings.get("-re"));
        initParameters(minThreshold, minFreq);
    }

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public List<SNVResultContainer> getHaplotypes() {
        long start = System.currentTimeMillis();
        processOverlaps(sample.reads);
        sample.reads.sort((r1, r2) -> r1.lOffset == r2.lOffset ? Integer.compare(r1.rOffset, r2.rOffset) : Integer.compare(r1.lOffset, r2.lOffset));
        log("Start 2SNV method");
        logSame("Compute profile");
        profile();
        log(" - DONE " + (System.currentTimeMillis() - start));
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        log(" - DONE " + (System.currentTimeMillis() - start));
        readAnswerHaplotypes();
        log("Compute cliques");
        //getMergedCluques first time to get cliques
        List<Set<Integer>> adjacencyList = getUnprocessedSnpsAdjecencyList();
        adjacencyList = rotateFrequentMinors(adjacencyList);
        List<Set<Integer>> cliques = getResolvedCliques(adjacencyList);

        log("Found cliques: " + cliques.size());
        if (cliques.size() < 50 && !cliques.isEmpty())
            log(cliques.stream().map(s -> new Clique(s, consensus()).toString()).reduce((t, s) -> t + "\n" + s).get());
        log(" - DONE " + (System.currentTimeMillis() - start));
        log("Start getting haplotypes");
        // divide by clusters and find haplotypes
        List<SNVResultContainer> totalResults = getHaplotypesFromCliques(start, cliques);
        log("Haplotypes before return " + totalResults.size());
        if (totalResults.size() == 0) {
            totalResults = getDefaultHaplotype();
            Start.errorCode = 5;
            Start.errorMessage = "All haplotypes got too low freqeuncy";
        }
        outputAnswerChecking(totalResults);
        //just for debug

        return totalResults;
    }

    /**
     * Unlike in previous version we will reconstruct haplotypes within windows of fragment length size
     * 1) We start finding haplotypes in a window SP, SE,
     * 2) add all edges in found haplotypes to form a clique
     * 3) extend window to include more edges until we cover the whole sample
     *
     * @return reconstructed haplotypes
     */
    public List<SNVResultContainer> getHaplotypesV2() {
        long start = System.currentTimeMillis();
        processOverlaps(sample.reads);
        sample.reads.sort((r1, r2) -> r1.lOffset == r2.lOffset ? Integer.compare(r1.rOffset, r2.rOffset) : Integer.compare(r1.lOffset, r2.lOffset));

        log("Start 2SNV method");
        logSame("Compute profile");
        profile();
        log(" - DONE " + (System.currentTimeMillis() - start));
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        log(" - DONE " + (System.currentTimeMillis() - start));
        int mostCoveredPosition = mostCoveredPosition();
        log("Most covered position is " + mostCoveredPosition + " with coverage " + struct.readsAtPosition[mostCoveredPosition].length);
        readAnswerHaplotypes();
        log("Compute cliques");
        sampleFragmentLength = calculateFragmentLength();
        if (sampleFragmentLength > 0.8 * (END_POSITION-START_POSITION)) {
            sampleFragmentLength = sample.referenceLength;
        }
        if (Start.settings.containsKey("-fl")) {
            sampleFragmentLength = Start.tryParseInt(Start.settings.get("-fl"), sampleFragmentLength);
        }
        log("Fragment length is " + sampleFragmentLength);
        int workingWindowSize = sampleFragmentLength;
        MAX_EDGES_NUMBER = 1_000_000; // we want to capture all the edges
        List<Set<Integer>> fullAdjacencyList = getUnprocessedSnpsAdjecencyList();
        fullAdjacencyList = rotateFrequentMinors(fullAdjacencyList);
        readAnswerHaplotypes();

        if (answer != null) {
            log("FN edges");
            Set<String> FNedges = new HashSet<>();
            for (Clique clique : answer) {
                for (Integer i : clique.splittedSnps) {
                    if (i / minorCount < START_POSITION || i / minorCount > END_POSITION) {
                        continue;
                    }
                    for (Integer j : clique.splittedSnps) {
                        if (j / minorCount < START_POSITION || j / minorCount > END_POSITION) {
                            continue;
                        }
                        if (i >= j || (j - i) / minorCount > sampleFragmentLength) {
                            continue;
                        }
                        if (al.charAt(getAllele(i)) == '-' && al.charAt(getAllele(j)) == '-') {
                            continue;
                        }
                        if (!FNedges.contains(i + "_" + j) && !fullAdjacencyList.get(i).contains(j)) {
                            FNedges.add(i + "_" + j);
                            int reads = commonReads[i / minorCount][j / minorCount];
                            log(i / minorCount + " " + j / minorCount + " " + al.charAt(getAllele(i)) + " " + al.charAt(getAllele(j)) + Arrays.toString(getOs(i, j)) + " " + getP(i, j)+ ((double)getOs(i, j)[3]/reads));
                        }
                        if (Start.settings.containsKey("-addFnedges")) {
                            fullAdjacencyList.get(i).add(j);
                            fullAdjacencyList.get(j).add(i);
                        }
                    }
                }
            }
            log("Total FN edges " + FNedges.size());
        }

        List<SNVResultContainer> totalResults = new ArrayList<>();
        //mostCoveredPosition = 0;
        int processedWindowStart = mostCoveredPosition; // we start from the same position and will extend either to the left or to the right
        int processedWindowEnd = mostCoveredPosition;
        // extend the processed window gradually
        while (processedWindowStart > START_POSITION || processedWindowEnd < END_POSITION) {
            int leftCoverage = processedWindowStart - workingWindowSize < START_POSITION ? 0 : struct.readsAtPosition[processedWindowStart - workingWindowSize].length;
            int rightCoverage = processedWindowEnd + workingWindowSize > END_POSITION ? 0 : struct.readsAtPosition[processedWindowEnd + workingWindowSize].length;
            boolean extendLeft = extendLeft(leftCoverage, rightCoverage, processedWindowStart, processedWindowEnd);
            if (extendLeft) { // extend to where we see higher coverage
                processedWindowStart -= workingWindowSize;
                processedWindowStart = Math.max(START_POSITION, processedWindowStart);
                int additionalLength = 0;
                if (processedWindowStart - START_POSITION < 20) { // no need to do tiny windows for last positions
                    additionalLength = processedWindowStart - START_POSITION;
                    processedWindowStart = START_POSITION;
                }
                workingWindowStart = processedWindowStart; // extend to the left
                workingWindowEnd = processedWindowStart + workingWindowSize + additionalLength;
                workingWindowEnd = Math.min(END_POSITION, workingWindowEnd);
                if (workingWindowEnd > processedWindowEnd) { // for the wirst iteration PWE = PWS so we need to extend both start and end
                    processedWindowEnd = workingWindowEnd;
                }
            } else {
                processedWindowEnd += workingWindowSize;
                processedWindowEnd = Math.min(END_POSITION, processedWindowEnd);
                int additionalLength = 0;
                if (END_POSITION - processedWindowEnd < 20) {
                    additionalLength = END_POSITION - processedWindowEnd;
                    processedWindowEnd = END_POSITION;
                }
                workingWindowEnd = processedWindowEnd;
                workingWindowStart = processedWindowEnd - workingWindowSize - additionalLength;
                workingWindowStart = Math.max(START_POSITION, workingWindowStart);
                if (workingWindowStart < processedWindowStart) {
                    processedWindowStart = workingWindowStart;
                }
            }
            log("Processing window [" + processedWindowStart + ", " + processedWindowEnd + "] window [" + workingWindowStart + ", " + workingWindowEnd + "] " + (extendLeft ? " extend left" : " extend right"));
            int splitWWStart = workingWindowStart * minorCount;
            int splitWWEnd = workingWindowEnd * minorCount;
            int splitPWStart = processedWindowStart * minorCount;
            int splitPWEnd = processedWindowEnd * minorCount;
            // create a new adjacency list for each window
            List<Set<Integer>> adjacencyList = new ArrayList<>();
            for (int i = 0; i < fullAdjacencyList.size(); i++) {
                adjacencyList.add(new HashSet<>());
                List<SNVResultContainer> fResults = totalResults;
                int fI = i;
                // add edges only in the range
                adjacencyList.get(i).addAll(fullAdjacencyList.get(i).stream().filter(
                        j -> (insideInterval(j, splitWWStart, splitWWEnd) && insideInterval(fI, splitWWStart, splitWWEnd)) // both inside interval
                                || edgeInHaplotype(fI, j, fResults)  //&& insideInterval(j, splitPWStart, splitPWEnd)// both in haplotype
                                || insideHaplotype(fI, fResults) && insideInterval(j, splitWWStart, splitWWEnd) // i in haplotype, j in window
                                || insideHaplotype(j, fResults) && insideInterval(fI, splitWWStart, splitWWEnd)) // j in haploty, i in window
                        .collect(Collectors.toList()));
            }
            List<Set<Integer>> cliques = getResolvedCliques(adjacencyList);
            log("Found resolved cliques: " + cliques.size());
            if (cliques.size() < 100 && !cliques.isEmpty())
                log(cliques.stream().map(s -> new Clique(s, consensus()).toString()).reduce((t, s) -> t + "\n" + s).get());
            if (Start.answer != null) log("Merged cliques missmatches");
            checkMissmatchesWithAnswer(cliques);
            log(" - DONE " + (System.currentTimeMillis() - start));
            log("Start getting haplotypes");

            totalResults = getHaplotypesFromCliques(start, cliques, splitWWStart, splitWWEnd, splitPWStart, splitPWEnd);
            // add all edges from found haplotypes inside the graph to form a clique in next iteration
            for (SNVResultContainer container : totalResults) {
                for (Integer i : container.sourceClique.splittedSnps) {

                    for (Integer j : container.sourceClique.splittedSnps) {
                        if (i.equals(j)) {
                            continue;
                        }
                        fullAdjacencyList.get(i).add(j);
                        fullAdjacencyList.get(j).add(i);
                    }
                }
            }
            outputAnswerChecking(totalResults, processedWindowStart, processedWindowEnd);
            log("Processed window [" + processedWindowStart + ", " + processedWindowEnd + "] window [" + workingWindowStart + ", " + workingWindowEnd + "] " + (extendLeft ? " extend left" : " extend right"));
            //log("Haplo cliques ("+totalResults.stream().+")");
            totalResults.forEach(t -> log(t.haploClique.toString()));
            log("");
        }
        totalResults = filterHaplotypeFrequencies(totalResults, HAPLOTYPE_CUT_THRESHOLD);
        outputAnswerChecking(totalResults, processedWindowStart, processedWindowEnd);
        return totalResults;
    }

    /**
     * Checks if the window should be extended left or right
     *
     * @param leftCoverage         coverage on the left position candidate to extend
     * @param rightCoverage        coverage on the right position candidate to extend
     * @param processedWindowStart already processed start position
     * @param processedWindowEnd   already processed end position
     * @return true if we should extend left, false otherwise
     */
    private boolean extendLeft(int leftCoverage, int rightCoverage, int processedWindowStart, int processedWindowEnd) {
        return leftCoverage > rightCoverage || processedWindowEnd >= END_POSITION; // if we reached the end we go to the left
    }

    private boolean insideInterval(int position, int left, int right) {
        return position >= left && position <= right;
    }

    /**
     * Checks if split position belongs to any haplotye
     *
     * @param position   split position
     * @param haplotypes a list of haplotypes
     * @return true if the position is present in at least on of the haplotypes
     */
    private boolean insideHaplotype(int position, List<SNVResultContainer> haplotypes) {
        for (SNVResultContainer haplotype : haplotypes) {
            if (haplotype.sourceClique.splittedSnps.contains(position)) {
                return true;
            }
        }
        return false;
    }


    /**
     * Tells if the edge exists in the given list of haplotypes. Since we work with windows we want to keep remote edges from existing haplotypes
     *
     * @param i          split position
     * @param j          split position
     * @param haplotypes list of previously found haplotypes
     * @return true if (i,j) belongs to any haplotype
     */
    private boolean edgeInHaplotype(int i, int j, List<SNVResultContainer> haplotypes) {
        for (SNVResultContainer haplotype : haplotypes) {
            if (haplotype.sourceClique.splittedSnps.contains(i) && haplotype.sourceClique.splittedSnps.contains(j)) {
                return true;
            }
        }
        return false;
    }

    private List<SNVResultContainer> getHaplotypesFromCliques(long start, List<Set<Integer>> cliques) {
        return getHaplotypesFromCliques(start, cliques, 0, sample.referenceLength * minorCount, 0, sample.referenceLength * minorCount);
    }


    private List<SNVResultContainer> getHaplotypesFromCliques(long start, List<Set<Integer>> cliques, int splitWorkingWindowStart, int splitWorkingWindowEnd, int splitProcessingWindowStart, int splitProcessingWindowEnd) {
        List<Set<Integer>> adjacencyList;// divide by clusters and find haplotypes
        List<SNVResultContainer> snvResultContainers = processCliques(cliques);
        log(" - DONE getting haplotypes " + (System.currentTimeMillis() - start));


        // calculate haplotypes for different edge frequency limits
        // so that we will make sure that higher frequency haplotypes won't be contaminated with the noise
        double[] edgesLimit = new double[0];//{0.1, 0.05, 0.01}; for now we don't want to have thise correction in any case
//        if (Start.settings.get("-fc").equals("false")) {
//            edgesLimit = new double[0];
//        }
        List<SNVResultContainer> totalResults = new ArrayList<>();
        log("Haplotypes before " + snvResultContainers.size());
        //outputAnswerChecking(snvResultContainers);
        // start with adding high frequency haplotypes
        for (double v : edgesLimit) {
            adjacencyList = getTopNEdgesAdjacencyList(v, splitWorkingWindowStart, splitWorkingWindowEnd);
            cliques = getResolvedCliques(adjacencyList);
            List<SNVResultContainer> haplotypes = processCliques(cliques);
            log("Haplo cliques:");
            for (SNVResultContainer haplotype : haplotypes) {
                log(haplotype.haploClique.toString());
                boolean added = false;
                for (SNVResultContainer resultHaplotype : totalResults) {
                    // we assume that if the distance between haplotypes is 1 then they are the same, maybe one SNP fell of at
                    // cluster building step. We trust the haplotype from higher frequency limit more so we keep it
                    if (Utils.hammingDistance(haplotype.haplotype, resultHaplotype.haplotype) <= 1) {
                        resultHaplotype.frequency += haplotype.frequency;
                        added = true;
                        break;
                    }
                }
                if (!added) {
                    totalResults.add(haplotype);
                }
            }
            log("Found haplotypes " + haplotypes.size());
            log("Haplotypes after " + v + " " + totalResults.size());
        }
        // add haplotypes obtained with all allowed edges
        for (SNVResultContainer haplotype : snvResultContainers) {
            boolean added = false;
            for (SNVResultContainer resultHaplotype : totalResults) {
                if (Utils.hammingDistance(haplotype.haplotype, resultHaplotype.haplotype) <= 1) {
                    resultHaplotype.frequency += haplotype.frequency;
                    added = true;
                    break;
                }
            }
            if (!added) {
                totalResults.add(haplotype);
            }
        }
        log("Haplotypes after no freq limit " + totalResults.size());
        int startPosition = splitProcessingWindowStart / minorCount;
        int endPosition = splitProcessingWindowEnd / minorCount;
        outputAnswerChecking(totalResults, startPosition, endPosition);
        log("");
        double[] freq = new double[totalResults.size()];
        for (int i = 0; i < freq.length; i++) {
            freq[i] = totalResults.get(i).frequency;
        }
        log("Freq before filtering" + Arrays.toString(freq));

        log("Freq before second filtering" + Arrays.toString(freq));
        //if (Start.settings.get("-ch").equals("true")) {
        totalResults = filterHaplotypeFrequencies(totalResults, HAPLOTYPE_CUT_THRESHOLD / 2);
        //}
        List<SNVResultContainer> filteredRecombinations = new ArrayList<>();
        Set<Set<Integer>> existingHaplotypes = new HashSet<>();
        //proccessing window start and end
        int PWS;
        int PWE;
        if (splitWorkingWindowStart == splitProcessingWindowStart) { // we extended window left
            PWS = splitWorkingWindowEnd; // from working end
            PWE = splitProcessingWindowEnd; // to processing end
        } else {
            PWS = splitProcessingWindowStart; // from processing start
            PWE = splitWorkingWindowStart; // to working start
        }
        if (PWS != PWE && Start.settings.containsKey("-filter")) { // not the first ineration
            for (SNVResultContainer container : totalResults) {
                Set<Integer> processedClique = container.sourceClique.splittedSnps.stream().filter(i -> i >= PWS && i <= PWE).collect(Collectors.toSet());
                if (!existingHaplotypes.contains(processedClique)) {
                    existingHaplotypes.add(processedClique);
                    filteredRecombinations.add(container);
                    log(" Haplotype " + " with freq " + container.frequency + " passed " + container.sourceClique);
                } else {
                    log(" Haplotype " + " with freq " + container.frequency + " was filtered out" + container.sourceClique);
                }

            }
            totalResults = filteredRecombinations;
        }
        return totalResults;
    }

    /**
     * Method for Variant Calling. Finds all SNPs and writes them into Clique object
     *
     * @return CLique with all found SNPs
     */
    public Clique getVC() {
        long start = System.currentTimeMillis();
        processOverlaps(sample.reads);
        sample.reads.sort((r1, r2) -> r1.lOffset == r2.lOffset ? Integer.compare(r1.rOffset, r2.rOffset) : Integer.compare(r1.lOffset, r2.lOffset));
        log("Start 2SNV method");
        logSame("Compute profile");
        log(" - DONE " + (System.currentTimeMillis() - start));
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        log(" - DONE " + (System.currentTimeMillis() - start));
        List<Set<Integer>> adjacencyList = getUnprocessedSnpsAdjecencyList();
        adjacencyList = rotateFrequentMinors(adjacencyList);
        log(" - DONE");
        //get haplotypes based on pure cliques, i.e. not merged.
        List<Set<Integer>> cliques = AlgorithmUtils.findCliques(adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());
        List<SNVResultContainer> snvResultContainers = processCliques(cliques);
        //IntStream.range(0, adjacencyList.size()).filter(i -> adjacencyList.get(i).size() > 0).boxed().collect(Collectors.toSet())
        return getVCCliqie(snvResultContainers.stream().map(c -> c.haploClique).flatMap(cl -> cl.splittedSnps.stream()).collect(Collectors.toSet()),
                TrueFrequencyEstimator.ILLUMINA_EPS);
    }

    @Override
    public double getP(int i, int j) {
        long[] os = getOs(i, j);
        if (os[0] == 0) {
            long tmp = os[0];
            os[0] = os[3];
            os[3] = tmp;
        }
        long reads = commonReads[i / minorCount][j / minorCount];
        double p = (os[1] * os[2]) / ((double) os[0] * reads);
        if (p > 1) {
            return 1;
        }

        if (os[3] / (MAX_READ_ERROR * os[0] + (1 - MAX_READ_ERROR) * (os[1] + os[2])) > MAX_READ_ERROR) {
            return 0;
        }
        if (p < 1E-12) {
            return 0;
        } else {
            return Utils.binomialOneMinusPvalue((int) os[3], p, (int) reads);
        }
    }

    @Override
    public double getNonEndgeP(int i, int j) {
        long[] os = getOs(i, j);
        if (os[0] == 0) {
            long tmp = os[0];
            os[0] = os[3];
            os[3] = tmp;
        }
        long reads = commonReads[i / minorCount][j / minorCount];
        double p = NON_EDGE_T_FREQ; //TODO change
        return Utils.binomialPvalue((int) os[3], p, (int) reads);

    }

    @Override
    protected long[] getOsImp(int i, int j) {
        int first = i / minorCount;
        int second = j / minorCount;
        long reads = commonReads[first][second];
        /*
         * false 1 means that in actual sample it has another minor or N in given position
         */
        long o22 = AlgorithmUtils.getSortedArraysIntersectionCount(struct.rowMinors[i], struct.rowMinors[j]);
        long o21 = struct.rowMinors[i].length; //all 2*
        long o12 = struct.rowMinors[j].length; //all *2
        // subtract 2N and false 21 from o21
        o21 = calculateO21(o21, second, i);
        //subtract N2 and false 12 from o12
        o12 = calculateO12(o12, first, j);
        //subtract 1N from reads
        long o11 = getO11(i, j, first, second, o22, o21, o12, reads);
        return new long[]{o11, o12, o21, o22};
    }

    @Override
    public long getCommonReadsCount(int i, int j) {
        return commonReads[i / minorCount][j / minorCount];
    }

    /**
     * In case when we just what to know snps, onlySnps argument should be true - useful in case we want just to clean reads and need snp positions
     */
    private List<Set<Integer>> getResolvedCliques(List<Set<Integer>> adjacencyList) {
        log("Edges in total: " + adjacencyList.stream().mapToInt(Set::size).sum() / 2);
        log("Not isolated vertices " + adjacencyList.stream().filter(l -> !l.isEmpty()).count());
        List<Set<Integer>> mergedCliques = getMergedCliques(adjacencyList);
        int conflictsResolved = 0;
        //resolve conflicts in cliques, leave only snps on the same position with higher frequency
        for (Set<Integer> clique : mergedCliques) {
            Iterator<Integer> it = clique.iterator();
            while (it.hasNext()) {
                Integer v = it.next();
                int row = v - v % minorCount;
                int allele = getAllele(v);
                for (int i = 0; i < minorCount; i++) {
                    if (row + i == v) {
                        continue;
                    }
                    if (clique.contains(row + i)) {
                        if (profile()[allele][v / minorCount] < profile()[getAllele(row + i)][v / minorCount]) {
                            it.remove();
                            conflictsResolved++;
                            break;
                        }
                    }
                }
            }
        }
        log("Conflicts in cliques resolved " + conflictsResolved);
        //List<Clique> collect = mergedCliques.stream().map(c -> new Clique(c, struct)).sorted(Comparator.comparingInt(c -> c.minors.length())).collect(Collectors.toList());
        return mergedCliques;
    }

    /**
     * In some cases we may need to change minor with major since they have close frequencies
     * Methods finds all minors with frequency > 45% swap with majors and check how many edges we will have
     * for both variants. Variant with most edges is kept. Consensus is updated accordingly
     *
     * @param adjacencyList adjacency list to rotate
     * @return new adjacency list with rotated minors
     */
    private List<Set<Integer>> rotateFrequentMinors(List<Set<Integer>> adjacencyList) {
        if (maxEdgesLimitReached) {
            log("No rotation because edges limit exceeded");
            return adjacencyList;
        }
        List<Integer> positions = new ArrayList<>();
        List<Character> oldChars = new ArrayList<>();
        List<Character> newChars = new ArrayList<>();
        List<Integer> oldSplitted = new ArrayList<>();
        List<Integer> newSplitted = new ArrayList<>();
        for (int i = 0; i < sample.referenceLength; i++) {
            //won't find correlation anyway
            if (struct.readsAtPosition[i].length < 40) {
                continue;
            }
            double[] tmp = new double[profile().length];
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = profile()[j][i];
            }
            Arrays.sort(tmp);
            if (tmp[tmp.length - 2] < 0.45) {
                continue;
            }
            char newChar = '_';
            for (int j = 0; j < profile().length; j++) {
                if (Math.abs(profile()[j][i] - tmp[tmp.length - 2]) < 0.000_000_001) {
                    newChar = al.charAt(j);
                    newChars.add(al.charAt(j));
                    break;
                }
            }
            if (newChar == 'N' || newChar == '-') {
                continue;
            }
            positions.add(i);
            oldSplitted.add(splitPosition(i, newChar)); // newChar is a minor before swap and oldChar will be minor after swap
            char oldChar = consensus.charAt(i);
            oldChars.add(oldChar);
            log("Exchange " + oldChar + " to " + newChar + " at " + i + " position");
            newSplitted.add(splitPosition(i, oldChar));
        }
        if (positions.size() == 0) {
            log("No rotation needed");
            return adjacencyList;
        }
        StringBuilder st = new StringBuilder(consensus);
        for (int i = 0; i < positions.size(); i++) {
            st.setCharAt(positions.get(i), newChars.get(i));
        }
        consensus = st.toString();
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        log(" - DONE ");
        List<Set<Integer>> newAdjacencyList = getUnprocessedSnpsAdjecencyList();
        for (int i = 0; i < positions.size(); i++) {
            int s1 = adjacencyList.get(oldSplitted.get(i)).size();
            int s2 = newAdjacencyList.get(newSplitted.get(i)).size();
            if (s1 < s2) {
                log(String.format("Keep exchange of %s to %s at %d position %d vs %d", oldChars.get(i), newChars.get(i), positions.get(i), s2, s1));
                adjacencyList = newAdjacencyList;
            } else {
                log(String.format("Don't keep exchange of %s to %s at %d position %d vs %d", oldChars.get(i), newChars.get(i), positions.get(i), s2, s1));
                st.setCharAt(positions.get(i), oldChars.get(i));
                consensus = st.toString();
            }
        }
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        adjacencyList = getUnprocessedSnpsAdjecencyList();
        log(" - DONE ");
        return adjacencyList;
    }

    private List<Set<Integer>> getUnprocessedSnpsAdjecencyList() {
        computeCommonReads();
        List<FindIlluminaEdgesParallelTask> tasks = new ArrayList<>();
        for (int i = 0; i < sample.referenceLength * minorCount; i++) {
            tasks.add(new FindIlluminaEdgesParallelTask(struct, i, this, sample, commonReads, Start.settings.containsKey("-ignoreDeletion")));
        }
        List<FindIlluminaEdgesParallelTask.EdgeSummary> allEdges = new ArrayList<>();
        boolean oldLog = log;
        if (Start.settings.containsKey("-noSNPlog")) {
            log = false;
        }
        try {
            List<Future<List<FindIlluminaEdgesParallelTask.EdgeSummary>>> futures = Start.executor.invokeAll(tasks);
            for (Future<List<FindIlluminaEdgesParallelTask.EdgeSummary>> future : futures) {
                allEdges.addAll(future.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        if (Start.settings.containsKey("-noSNPlog")) {
            log = oldLog;
        }
        log("Total edges before filtering " + allEdges.size());
        List<Set<Integer>> adjacencyList = convertEdgesIntoAdjacencyList(allEdges);
        if (!Start.settings.containsKey("-ignoreDeletion")){
            addDeletionEdges(adjacencyList);
        }
        noLimitEdges = allEdges;
        removeEdgesForSecondMinors(adjacencyList, struct);
        return adjacencyList;
    }

    /**
     * Some haplotypes have long deletions (a factor of 3). Since we have a restriction on the minimum distance between SNPs it might be a problem
     * So we preemptively add all edges to create a clique between deletions that are supported by enough reads
     *
     * @param adjacencyList adjacency list
     */
    private void addDeletionEdges(List<Set<Integer>> adjacencyList) {
        log("Start adding deletion edges");
        Map<Interval, Integer> deletionMap = new HashMap<>();
        List<Callable<Map<Interval, Integer>>> tasks = new ArrayList<>();
        int threads = Start.threadsNumber();
        for (int i = 0; i < threads; i++) {
            tasks.add(new DeletionEdgesParallelTask(i, threads, sample, MIN_O22_THRESHOLD));
        }
        try {
            List<Future<Map<Interval, Integer>>> futures = Start.executor.invokeAll(tasks);
            for (Future<Map<Interval, Integer>> future : futures) {
                Map<Interval, Integer> pairIntegerMap = future.get();
                pairIntegerMap.forEach((k, v) -> {
                    Integer count = deletionMap.getOrDefault(k, 0);
                    deletionMap.put(k, count + v);
                });
            }
        } catch (InterruptedException | ExecutionException e) {
            System.out.println("Finding deletion edges cased an exception");
            e.printStackTrace();
        }
        AtomicInteger totalDeletionRegions = new AtomicInteger();
        //deletionMap.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach();
        deletionMap.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach((entry) -> {
            Integer nReads = entry.getValue();
            Interval region = entry.getKey();
            if (nReads >= MIN_O22_THRESHOLD) {
                totalDeletionRegions.getAndIncrement();
                log("Deletion region " + region + " supported by " + nReads + " reads, length " + (region.e - region.s + 1));
                for (int i = region.s; i < region.e; i++) {
                    for (int j = i + 1; j < region.e; j++) {
                        int l = splitPosition(i, '-');
                        int r = splitPosition(j, '-');
                        adjacencyList.get(l).add(r);
                        adjacencyList.get(r).add(l);
                    }
                }
            }
        });
    }

    /**
     * Returns the adjacency list composed of n the most frequent edges
     */
    private List<Set<Integer>> getTopNEdgesAdjacencyList(double freq) {
        return convertEdgesIntoAdjacencyList(noLimitEdges.stream().filter(e -> e.relFreq > freq).collect(Collectors.toList()));
    }

    /**
     * Returns the adjacency list composed of n the most frequent edges
     */
    private List<Set<Integer>> getTopNEdgesAdjacencyList(double freq, int splitWorkingWindowStart, int splitWorkingWindowEnd) {
        return convertEdgesIntoAdjacencyList(noLimitEdges.stream().filter(e -> e.relFreq > freq && e.i >= splitWorkingWindowStart && e.i <= splitWorkingWindowEnd
                && e.j >= splitWorkingWindowStart && e.j <= splitWorkingWindowEnd).collect(Collectors.toList()));
    }

    /**
     * We don't want to process more than 700 edges, so we pick only 700 with the best support( highest o22/coverage value)
     *
     * @param allEdges
     */
    private List<Set<Integer>> convertEdgesIntoAdjacencyList(List<FindIlluminaEdgesParallelTask.EdgeSummary> allEdges) {
        List<Set<Integer>> adjacencyList = new ArrayList<>();
        for (int i = 0; i < sample.referenceLength * minorCount; i++) {
            adjacencyList.add(new HashSet<>());
        }
        allEdges.sort((e1, e2) -> -Double.compare(e1.relFreq, e2.relFreq));

        int maxEdges = Math.min(allEdges.size(), MAX_EDGES_NUMBER);
        maxEdgesLimitReached = maxEdges == MAX_EDGES_NUMBER;
        for (int i = 0; i < maxEdges; i++) {
            FindIlluminaEdgesParallelTask.EdgeSummary edge = allEdges.get(i);
            adjacencyList.get(edge.i).add(edge.j);
            adjacencyList.get(edge.j).add(edge.i);
        }
        return adjacencyList;
    }

    private void computeCommonReads() {
        if (commonReads == null) {
            commonReads = new int[sample.referenceLength][sample.referenceLength];
            log("Common reads matrix calculation");
            List<Callable<Boolean>> tasks = new ArrayList<>();

            for (int i = 0; i < sample.referenceLength; i++) {
                tasks.add(new CommonReadsIlluminaParallelTask(i, commonReads, struct, sample, START_POSITION, END_POSITION, log));
            }
            try {
                List<Future<Boolean>> futures = Start.executor.invokeAll(tasks);
                futures.forEach(future -> {
                    try {
                        future.get();
                    } catch (InterruptedException | ExecutionException e) {
                        System.err.println("Error! Parallel tasks were not successful on get");
                        e.printStackTrace();
                    }
                });
            } catch (InterruptedException e) {
                System.err.println("Error! Parallel tasks were not successful on invoke");
                e.printStackTrace();
            }
            log("");
        }
    }

    @Override
    public String consensus() {
        if (consensus == null) {
            consensus = Utils.consensus(profile(), al, false);// we don't want to have deletions in consensus since it causes problems when some haplotypes have deletions
            log(" Reference length = " + consensus.length());// we want deletion to always be a minor
        }
        return consensus;
    }

    @Override
    public double[][] profile() {
        if (profile == null) {
            profile = Utils.profile(sample, al);
        }
        return profile;
    }

    public long getO11(int i, int j, int first, int second, long o22, long o21, long o12, long reads) {
        long o11 = reads - o12 - o21 - o22; //all 11(and some false positive 11),12
        //1 1'
        if (struct.readsAtPosition[second].length == 0 || struct.readsAtPosition[first].length == 0) {
            return 0;
        }
        for (int k = 0; k < minorCount; k++) {
            if (k == j % minorCount) {
                continue;
            }
            int currentMinor = j - j % minorCount + k;
            for (int m = AlgorithmUtils.binarySearch(struct.rowMinors[currentMinor], struct.readsAtPosition[second][0]); m < struct.rowMinors[currentMinor].length; m++) {
                PairEndRead read = sample.reads.get(struct.rowMinors[currentMinor][m]);
                //find in what part of read is first
                if (read.lOffset <= first && read.lOffset + read.l.length() > first) {
                    if (read.l.charAt(first - read.lOffset) == consensus().charAt(first)) {
                        o11--;
                    }
                } else if (read.rOffset <= first && read.rOffset + read.r.length() > first) {
                    if (read.r.charAt(first - read.rOffset) == consensus().charAt(first)) {
                        o11--;
                    }
                }
            }

        }
        //1' 1
        for (int k = 0; k < minorCount; k++) {
            if (k == i % minorCount) {
                continue;
            }
            int currentMinor = i - i % minorCount + k;
            for (int m = AlgorithmUtils.binarySearch(struct.rowMinors[currentMinor], struct.readsAtPosition[first][0]); m < struct.rowMinors[currentMinor].length; m++) {
                PairEndRead read = sample.reads.get(struct.rowMinors[currentMinor][m]);
                //find in what part of read is second
                if (read.lOffset <= second && read.lOffset + read.l.length() > second) {
                    if (read.l.charAt(second - read.lOffset) == consensus().charAt(second)) {
                        o11--;
                    }
                } else if (read.rOffset <= second && read.rOffset + read.r.length() > second) {
                    if (read.r.charAt(second - read.rOffset) == consensus().charAt(second)) {
                        o11--;
                    }
                }
            }
        }
        //1' 1', 1' 2, 2 1'
        for (int k = 0; k < minorCount; k++) {
            if (struct.rowMinors[i - i % minorCount + k].length == 0) {
                continue;
            }
            for (int l = 0; l < minorCount; l++) {
                if ((l == j % minorCount && k == i % minorCount) || struct.rowMinors[j - j % minorCount + l].length == 0) {
                    continue;
                }
                o11 -= AlgorithmUtils.getSortedArraysIntersectionCount(struct.rowMinors[i - i % minorCount + k], struct.rowMinors[j - j % minorCount + l]);
            }
        }
        return o11;
    }

    /**
     * If left and right reads overlap - merge them into just left read
     */
    private void processOverlaps(List<PairEndRead> reads) {
        for (PairEndRead read : reads) {
            if (read.lOffset + read.l.length() > read.rOffset && read.rOffset != -1) {
                StringBuilder newL = new StringBuilder(read.l);
                //maybe r lies entirely in l, then no copying
                if (read.lOffset + read.l.length() - read.rOffset < read.r.length()) {
                    newL.append(read.r.substring(read.lOffset + read.l.length() - read.rOffset));
                }
                read.l = newL.toString();
                read.r = "";
                read.rOffset = -1;
            }
        }
    }

    private boolean readHasPosition(PairEndRead read, int position) {
        return (read.lOffset <= position && read.lOffset + read.l.length() > position)
                || (read.rOffset <= position && read.rOffset + read.r.length() > position);
    }

    private boolean readHasCharAtPosition(PairEndRead read, int position, char m) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position) {
            return read.l.charAt(position - read.lOffset) == m;
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position) {
            return read.r.charAt(position - read.rOffset) == m;
        }
        return false;
    }

    private char readCharAtPosition(PairEndRead read, int position) {
        if (read.lOffset <= position && read.lOffset + read.l.length() > position) {
            return read.l.charAt(position - read.lOffset);
        } else if (read.rOffset <= position && read.rOffset + read.r.length() > position) {
            return read.r.charAt(position - read.rOffset);
        }
        return 'N';
    }

    /**
     * Methods process computed cliques. At first it separate reads into clusters and then compute consensus for each cluster and some additional information
     *
     * @param cliques Set of cliques in form of splitted SNPs(where position and minor are encoded in a single number)
     * @return Set of containers with results. Each container contains haplotype itself and some additional helpful information
     */
    private List<SNVResultContainer> processCliques(List<Set<Integer>> cliques) {
        //add consensus clique
        if (!consensusWasRemoved) {
            cliques.add(new HashSet<>());
        }
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, consensus())));
        log("Start build clusters");
        Map<String, List<Pair<PairEndRead, Integer>>> clustersRelative = buildClustersRelative(allPositionsInCliques, allCliquesCharacters);
        log(" - DONE");
        /*
            simple estimation of haplotypes frequencies through the number of assigned reads. We don't want too much hapltoypes since it doesn't make much sense
            take only top MAX_HAPLOTYPES_BEFORE_EM
         */
        if (clustersRelative.size() > MAX_HAPLOTYPES_BEFORE_EM) {
            List<Pair<Double, Map.Entry<String, List<Pair<PairEndRead, Integer>>>>> freqList = new ArrayList<>();
            for (Map.Entry<String, List<Pair<PairEndRead, Integer>>> entry : clustersRelative.entrySet()) {
                double freq = 0;
                for (Pair<PairEndRead, Integer> p : entry.getValue()) {
                    freq += 1. / p.getValue();
                }
                freqList.add(new Pair<>(freq, entry));
            }
            // sort descending
            freqList.sort((l, r) -> r.getKey().compareTo(l.getKey()));
            clustersRelative = new HashMap<>();
            for (int i = 0; i < MAX_HAPLOTYPES_BEFORE_EM; i++) {
                clustersRelative.put(freqList.get(i).getValue().getKey(), freqList.get(i).getValue().getValue());
            }
        }
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        List<SNVResultContainer> haplotypes = clustersRelative.entrySet().parallelStream().filter(s -> s.getValue().size() > 10).map(s -> {
            List<Pair<PairEndRead, Integer>> cluster = s.getValue();
            List<PairEndRead> reads = cluster.stream().map(Pair::getKey).collect(Collectors.toCollection(ArrayList::new));
            IlluminaSNVSample snvSample = new IlluminaSNVSample("tmp", reads, sample.referenceLength, sample.singleRead);
            //if there is no reads, put consensus there
//            double[][] haploProfile = Utils.profile(snvSample,  al);
            double[][] haploProfile = Utils.profile(snvSample, cluster.stream().mapToInt(Pair::getValue).toArray(), al);
            for (int i = 0; i < haploProfile[0].length; i++) {
                double max = 0;
                for (double[] aProfile : haploProfile) {
                    if (aProfile[i] > max) {
                        max = aProfile[i];
                    }
                }
                if (!(max > 0.01)) {
                    int i1 = al.indexOf(consensus().charAt(i)) == -1 ? 'N' : al.indexOf(consensus().charAt(i));
                    haploProfile[i1][i] = 1;
                }
            }
            String haplotype = Utils.consensus(haploProfile, al);
            Clique sourceClique = getSourceClique(allPositionsInCliques, cliquesSet, s.getKey());
            //rare case where all reads have N on a certain position
            StringBuilder str = new StringBuilder(haplotype);
            int[][] count = Utils.countCoverage(snvSample, al);
            int[] coverage = new int[str.length()];
            for (int i = 0; i < coverage.length; i++) {
                for (int j = 0; j < count.length; j++) {
                    coverage[i] += count[j][i];
                }
            }
            for (int i = 0; i < str.length(); i++) {
                if (str.charAt(i) == 'N') {
                    str.setCharAt(i, consensus().charAt(i));
                }
                // some reads may be aligned wrongly or just have minor allele with tiny coverage(1-2 reads) that will give an SNP there
                // for consensus haplotype we want an additional check since reads  with consensus base can go to another clique and for consensus we will assign only "trash" reads
                if (str.charAt(i) != consensus.charAt(i) && coverage[i] < MINIMUM_COVERAGE_FOR_HAPLOTYPE_SNP
                        || (str.charAt(i) != consensus.charAt(i) && sourceClique.snps.size() == 0 && count[al.indexOf(str.charAt(i))][i] < MINIMUM_COVERAGE_FOR_CONSENSUS_SNP)
                ) {
                    //log("prevented at pos " + i + " to get " + str.charAt(i) + " with coverage " + coverage[i]);
                    str.setCharAt(i, consensus().charAt(i));
                }
            }
            haplotype = str.toString();
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < haplotype.length(); i++) {
                if (haplotype.charAt(i) != consensus().charAt(i)) {
                    snps.add(splitPosition(i, haplotype.charAt(i)));
                }
            }
            Clique haplotypeClique = new Clique(snps, consensus());
            SNVResultContainer container = new SNVResultContainer(s.getKey(), reads, haplotypeClique, haplotype);
            container.sourceClique = sourceClique;
            return container;
        }).collect(Collectors.toList());
        Map<String, SNVResultContainer> merge = new HashMap<>();
        haplotypes.forEach(h -> {
            if (!merge.containsKey(h.haplotype)) {
                merge.put(h.haplotype, h);
            } else {
                merge.get(h.haplotype).illuminaCluster.addAll(h.illuminaCluster);
            }
        });
        List<SNVResultContainer> result = new ArrayList<>(merge.values());
        List<String> h = result.stream().map(s -> s.haplotype).collect(Collectors.toList());
        log("Start EM ");
        IlluminaEM illuminaEM = new IlluminaEM();
        List<Double> frequencies = illuminaEM.frequencies(h, sample);
        // there is no evidence that consensus haplotype exists if it's frequency is less than 25%
        // remove it and recalculate frequencies
        for (int i = 0; i < result.size(); i++) {
            if ((result.get(i).sourceClique.splittedSnps.isEmpty() || result.get(i).haploClique.splittedSnps.isEmpty()) && frequencies.get(i) < 0.25) {
                result.remove(i);
                h.remove(i);
                log("Consensus haplotype was removed");
                frequencies = new IlluminaEM().frequencies(h, sample);
                consensusWasRemoved = true;
                break;
            }
        }
        log(" - DONE ");
        for (int i = 0; i < frequencies.size(); i++) {
            result.get(i).frequency = frequencies.get(i);
        }
        return result.stream().filter(ha -> ha.frequency > 8e-4).sorted((s1, s2) -> -Double.compare(s1.frequency, s2.frequency)).collect(Collectors.toList());
    }

    /**
     * Build read clusters based on cliques. Each read will go to nearest clique in terms of Hamming distance divided by the number of nearest cliques
     *
     * @param allPositionsInCliques sorted array or all positions with at least one clique
     * @param allCliquesCharacters  characters in cliques according to allPositionsInCliques.
     *                              Has consensus allele if clique doesn't include particular position from allPositionsInCliques
     * @return Map with clusters, where key is string of clique characters, value is a set of reads
     */
    private Map<String, List<Pair<PairEndRead, Integer>>> buildClustersRelative(List<Integer> allPositionsInCliques, List<String> allCliquesCharacters) {
        Map<String, List<Pair<PairEndRead, Integer>>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new ArrayList<>()));
        if (allCliquesCharacters.size() == 1) {
            for (PairEndRead read : sample.reads) {
                clusters.get(allCliquesCharacters.get(0)).add(new Pair<>(read, 1));
            }

            return clusters;
        }
        int consensusCliqueIndex = 0;
        for (int i = 0; i < allCliquesCharacters.size(); i++) {
            boolean fl = true;
            String s = allCliquesCharacters.get(i);
            for (int j = 0; j < s.length(); j++) {
                if (s.charAt(j) != consensus().charAt(allPositionsInCliques.get(j))) {
                    fl = false;
                    break;
                }
            }
            if (fl) {
                consensusCliqueIndex = i;
            }
        }
        List<Callable<List<Integer>>> tasks = new ArrayList<>();
        AtomicInteger iter = new AtomicInteger();
        for (PairEndRead read : sample.reads) {
            int finalConsensusCliqueIndex = consensusCliqueIndex;
            tasks.add(() -> {
                int min = 1_000_000;
                List<Integer> minI = new ArrayList<>();
                // for each clique find the distance and overlap with the read
                for (int i = 0; i < allCliquesCharacters.size(); i++) {
                    int d = 0;
                    int coincidences = 0; //overlap in positions
                    boolean isConsensusClique = i == finalConsensusCliqueIndex;
                    String c = allCliquesCharacters.get(i);
                    for (int j = 0; j < allPositionsInCliques.size(); j++) {
                        if (allPositionsInCliques.get(j) < read.lOffset) {
                            continue;
                        }
                        if (allPositionsInCliques.get(j) > read.lOffset + read.l.length() && allPositionsInCliques.get(j) > read.rOffset + read.r.length()) {
                            break;
                        }
                        char charAtPosition = readCharAtPosition(read, allPositionsInCliques.get(j));
                        //increase distance
                        if (charAtPosition != 'N') {
                            if (charAtPosition != c.charAt(j)) {
                                d++;
                            }
                            //coincidence with current clique
                            if (c.charAt(j) != consensus().charAt(allPositionsInCliques.get(j))) {
                                coincidences++;
                            }
                            //coincidence with any cluque snp. Only for consensus clique
                            if (isConsensusClique) {
                                coincidences++;
                            }
                        }
                    }
                    if (coincidences != d && d < min && coincidences > 0) {
                        min = d;
                        minI = new ArrayList<>();
                    }
                    if (d == min && coincidences > 0) {
                        minI.add(i);
                    }
                }
                if (iter.incrementAndGet() % 100000 == 0) {
                    logSame("\r" + iter.get());
                }
                return minI;
            });

        }
        try {
            List<Future<List<Integer>>> futures = Start.executor.invokeAll(tasks);
            for (int i = 0; i < sample.reads.size(); i++) {
                List<Integer> minI = futures.get(i).get();
                PairEndRead read = sample.reads.get(i);
                if (!minI.isEmpty()) {
                    minI.forEach(e -> clusters.get(allCliquesCharacters.get(e)).add(new Pair<>(read, minI.size())));
                }
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        log("");
        return clusters;
    }

    public long calculateO12(long o12, int first, int j) {
        //subtract N2 and false 12 from o12
        for (int k = 0; k < struct.rowMinors[j].length; k++) {
            PairEndRead read = sample.reads.get(struct.rowMinors[j][k]);
            //find in what part of read is first
            if (read.lOffset <= first && read.lOffset + read.l.length() > first) {
                if (read.l.charAt(first - read.lOffset) != consensus().charAt(first)) {
                    o12--;
                }
            } else if (read.rOffset <= first && read.rOffset + read.r.length() > first) {
                if (read.r.charAt(first - read.rOffset) != consensus().charAt(first)) {
                    o12--;
                }
            } else {
                o12--;//first doesn't have this position at all
            }
        }
        return o12;
    }

    public long calculateO21(long o21, int second, int i) {
        for (int k = 0; k < struct.rowMinors[i].length; k++) {
            PairEndRead read = sample.reads.get(struct.rowMinors[i][k]);
            //find in what part of read is second
            if (read.lOffset <= second && read.lOffset + read.l.length() > second) {
                if (read.l.charAt(second - read.lOffset) != consensus().charAt(second)) {
                    o21--;
                }
            } else if (read.rOffset <= second && read.rOffset + read.r.length() > second) {
                if (read.r.charAt(second - read.rOffset) != consensus().charAt(second)) {
                    o21--;
                }
            } else {
                o21--;//second doesn't have this position at all
            }

        }
        return o21;
    }


    private int hitsInKnownHaplotypes(int first, char m1, int second, char m2) {
        String key = first + "_" + second + m1 + m2;
        if (knownHits.containsKey(key)) {
            return knownHits.get(key);
        }
        int result = 0;
        for (int i = 0; i < Start.coronaHaplotypes.reads.length; i++) {
            String read = Start.coronaHaplotypes.reads[i];
            if (second < read.length() - 1 && read.charAt(first) == m1 && read.charAt(second) == m2) {
                result++;
            }
        }
        knownHits.put(key, result);
        return result;
    }

    /**
     * Calculates the fragment length of the sample - the length of read fragment in 90-th percentile
     *
     * @return fragment length
     */
    protected int calculateFragmentLength() {
        int[] fragmentCount = new int[sample.referenceLength + 1];
        for (PairEndRead read : sample.reads) {
            if (read.rOffset != -1) {
                int length = Math.min(sample.referenceLength - 1, read.rOffset + read.r.length() - read.lOffset);
                fragmentCount[length]++;
            } else {
                fragmentCount[read.l.length()]++;
            }
        }
        int count = 0;
        int readsCount = sample.reads.size();
        for (int i = 0; i < fragmentCount.length; i++) {
            count += fragmentCount[i];
            if (count > readsCount * 0.9) {
                return i;
            }
        }
        return sample.referenceLength - 1;
    }

    /**
     * Calculates median most covered position
     *
     * @return median most covered position
     */
    private int mostCoveredPosition() {
        int maxCoverege = Arrays.stream(struct.readsAtPosition).skip(START_POSITION).limit(END_POSITION - START_POSITION + 1).mapToInt(r -> r.length).max().getAsInt();
        int[] maxCoveragePositions = IntStream.range(0, struct.readsAtPosition.length).filter(i -> struct.readsAtPosition[i].length == maxCoverege).toArray();
        return maxCoveragePositions[maxCoveragePositions.length / 2];
    }

    @Override
    protected List<PairEndRead> getIlluminaReads() {
        return sample.reads;
    }

    @Override
    protected Map<Integer, String> getPacBioCluster() {
        return null;
    }
}
