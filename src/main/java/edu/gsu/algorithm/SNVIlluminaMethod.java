package edu.gsu.algorithm;

import edu.gsu.algorithm.em.IlluminaEM;
import edu.gsu.algorithm.util.AlgorithmUtils;
import edu.gsu.algorithm.util.CommonReadsIlluminaParallelTask;
import edu.gsu.algorithm.util.FindIlluminaEdgesParallelTask;
import edu.gsu.algorithm.util.TrueFrequencyEstimator;
import edu.gsu.model.Clique;
import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.SNVStructure;
import edu.gsu.start.Start;
import edu.gsu.util.Utils;
import edu.gsu.util.builders.SNVStructureBuilder;

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

public class SNVIlluminaMethod extends AbstractSNV {
    public static final int minorCount = al.length() - 2;
    public int MAX_EDGES_NUMBER;
    private long[][] commonReads;
    private boolean maxEdgesLimitReached = false;
    // we assume them to be sorted
    private List<FindIlluminaEdgesParallelTask.EdgeSummary> noLimitEdges;

    @Override
    protected Technology technology() {
        return Technology.ILLUMINA;
    }

    private final IlluminaSNVSample sample;
    private SNVStructure struct;
    private double profile[][];
    private String consensus;

    /*private class CorrelationContainer {
        long o11;
        long o12;
        long o21;
        long o22;
        long reads;

        CorrelationContainer(long o11, long o12, long o21, long o22, long reads) {
            this.o11 = o11;
            this.o12 = o12;
            this.o21 = o21;
            this.o22 = o22;
            this.reads = reads;
        }
    }*/

    //private Map<String, CorrelationContainer> correlationMap;

    /**
     * @param sample Given input reads
     * @param log
     */
    public SNVIlluminaMethod(IlluminaSNVSample sample, boolean log) {
        super(100, Double.parseDouble(Start.settings.getOrDefault("-tf", "0.05")));
        this.log = log;
        this.sample = sample;
        this.MAX_EDGES_NUMBER = Start.settings.get("-el") == null ? 700 : Integer.parseInt(Start.settings.get("-el"));
        this.MAX_READ_ERROR = Start.settings.get("-re") == null ? 0.2 : Double.parseDouble(Start.settings.get("-re"));
    }

    /**
     * @param sample       Given input reads
     * @param minThreshold threshold for O22 value to consider alleles as SNPs
     * @param minFreq      minimum frequency (relative to reads' coverage) for O22 value to consider alleles as SNP
     * @param log
     */
    public SNVIlluminaMethod(IlluminaSNVSample sample, int minThreshold, double minFreq, boolean log) {
        super(minThreshold, minFreq);
        this.log = log;
        this.sample = sample;
        this.MAX_EDGES_NUMBER = Start.settings.get("-el") == null ? 700 : Integer.parseInt(Start.settings.get("-el"));
        this.MAX_READ_ERROR = Start.settings.get("-re") == null ? 0.2 : Double.parseDouble(Start.settings.get("-re"));
    }

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public List<SNVResultContainer> getHaplotypes() {
        long start = System.currentTimeMillis();
        sample.reads = processOverlaps(sample.reads);
        sample.reads.sort((r1, r2) -> r1.lOffset == r2.lOffset ? Integer.compare(r1.rOffset, r2.rOffset) : Integer.compare(r1.lOffset, r2.lOffset));
        log("Start 2SNV method");
        logSame("Compute profile");
        profile();
        log(" - DONE " + (System.currentTimeMillis() - start));
        logSame("Compute SNV data structure");
        struct = SNVStructureBuilder.buildIllumina(sample, consensus(), profile());
        log(" - DONE " + (System.currentTimeMillis() - start));
        log("Compute cliques");
        //getMergedCluques first time to get cliques
        List<Set<Integer>> adjacencyList = getUnprocessedSnpsAdjecencyList();
        adjacencyList = rotateFrequentMinors(adjacencyList);
        Set<Set<Integer>> cliques = getResolvedCliques(adjacencyList);

        log("Found cliques: " + cliques.size());
        if (cliques.size() < 50 && !cliques.isEmpty())
            log(cliques.stream().map(s -> new Clique(s, consensus()).toString()).reduce((t, s) -> t + "\n" + s).get());
        log(" - DONE " + (System.currentTimeMillis() - start));
        log("Start getting haplotypes");
        // divide by clusters and find haplotypes
        List<SNVResultContainer> snvResultContainers = processCliques(cliques);
        log(" - DONE getting haplotypes " + (System.currentTimeMillis() - start));


        // calculate haplotypes for different edge frequency limits
        // so that we will make sure that higher frequency haplotypes won't be contaminated with the noise
        double[] edgesLimit = {0.1, 0.05, 0.01};
        List<SNVResultContainer> totalResults = new ArrayList<>();
        log("Haplotypes before " + snvResultContainers.size());
        readAnswerHaplotypes();
        outputAnswerChecking(snvResultContainers);
        // start with adding high frequency haplotypes
        for (int i = 0; i < edgesLimit.length; i++) {
            log = false;
            adjacencyList = getTopNEdgesAdjacencyList(edgesLimit[i]);
            cliques = getResolvedCliques(adjacencyList);
            log = true;
            List<SNVResultContainer> haplotypes = processCliques(cliques);

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
            log("Haplotypes after " + edgesLimit[i] + " " + totalResults.size());
            outputAnswerChecking(totalResults);
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
        outputAnswerChecking(totalResults);
        log("");
        double[] freq = new double[totalResults.size()];
        for (int i = 0; i < freq.length; i++) {
            freq[i] = totalResults.get(i).frequency;
        }
        freq = Utils.normalize(freq);
        for (int i = 0; i < freq.length; i++) {
            totalResults.get(i).frequency = freq[i];
        }
        log("Freq before filtering" + Arrays.toString(freq));
        totalResults = totalResults.stream().filter(ha -> ha.frequency > HAPLOTYPE_CUT_THRESHOLD).sorted((s1, s2) -> -Double.compare(s1.frequency, s2.frequency)).collect(Collectors.toList());
        log("Haplotypes before return " + totalResults.size());
        outputAnswerChecking(totalResults);
        //just for debug

        return totalResults;
    }

    /**
     * Method for Variant Calling. Finds all SNPs and writes them into Clique object
     *
     * @return CLique with all found SNPs
     */
    public Clique getVC() {
        long start = System.currentTimeMillis();
        sample.reads = processOverlaps(sample.reads);
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
        Set<Set<Integer>> cliques = AlgorithmUtils.findCliques(adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toSet());
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
            return USE_LOG_PVALUE ? Utils.binomialLogPvalue((int) os[3], p, (int) reads) : Utils.binomialPvalue((int) os[3], p, (int) reads);
        }
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
    private Set<Set<Integer>> getResolvedCliques(List<Set<Integer>> adjacencyList) {
        log("Edges in total: " + adjacencyList.stream().mapToInt(Set::size).sum());
        log("Not isolated vertices " + adjacencyList.stream().filter(l -> !l.isEmpty()).count());
        Set<Set<Integer>> mergedCliques = getMergedCliques(adjacencyList);
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
     * for both variants. Variant with most edges is kept. COnsensus is updated accordingly
     *
     * @param adjacencyList
     * @return
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
            if (newChar == 'N') {
                continue;
            }
            positions.add(i);
            oldSplitted.add(splittedPosition(i, newChar)); // newChar is a minor before swap and oldChar will be minor after swap
            char oldChar = consensus.charAt(i);
            oldChars.add(oldChar);
            log("Exchange " + oldChar + " to " + newChar + " at " + i + " position");
            newSplitted.add(splittedPosition(i, oldChar));
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
        try {
            List<Future<List<FindIlluminaEdgesParallelTask.EdgeSummary>>> futures = Start.executor.invokeAll(tasks);
            for (Future<List<FindIlluminaEdgesParallelTask.EdgeSummary>> future : futures) {
                allEdges.addAll(future.get());
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        log("Total edges before filtering " + allEdges.size());
        List<Set<Integer>> adjacencyList = convertEdgesIntoAdjacencyList(allEdges);
        noLimitEdges = allEdges;
        removeEdgesForSecondMinors(adjacencyList, struct);
        return adjacencyList;
    }

    /**
     * Returns the adjacency list composed of n the most frequent edges
     */
    private List<Set<Integer>> getTopNEdgesAdjacencyList(double freq) {
        return convertEdgesIntoAdjacencyList(noLimitEdges.stream().filter(e -> e.relFreq > freq).collect(Collectors.toList()));
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
        log("Top edges");
        allEdges.sort((e1, e2) -> -Double.compare(e1.relFreq, e2.relFreq));
        if (allEdges.size() > 10) {
            int edges = Math.min(allEdges.size(), MAX_EDGES_NUMBER);
            for (int i = 0; i < edges; i++) {
                log(allEdges.get(i).toString());
            }
            log("Bottom edges");
            for (int i = 0; i < 10; i++) {
                log(allEdges.get(allEdges.size() - 10 + i).toString());
            }
        }

        int maxEdges = Math.min(allEdges.size(), MAX_EDGES_NUMBER);
        maxEdgesLimitReached = maxEdges == MAX_EDGES_NUMBER;
        for (int i = 0; i < maxEdges; i++) {
            FindIlluminaEdgesParallelTask.EdgeSummary edge = allEdges.get(i);
            adjacencyList.get(edge.i).add(edge.j);
            adjacencyList.get(edge.j).add(edge.i);
        }
        return adjacencyList;
    }

    public void processZeroO(Set<Integer> adjacencyNode, int i, int l, int first, int[] hits, int j, int second, long reads, char m1, char m2, long o22, long o21, long o12, long o11) {
        if (o22 / (MAX_READ_ERROR * o11 + (1 - MAX_READ_ERROR) * (o12 + o21)) < MAX_READ_ERROR) {
            return;
        }
        log(String.format("%d %d %c %c m1=%d m2=%d hits=%d p012=%.3e reads=%d %d %d %d %d",
                first, second, m1, m2, l, struct.rowMinors[j].length, hits[j], 0.0, reads, o11, o12, o21, o22));
        //correlationMap.put(getCorrelationKey(first, second, m1, m2), new CorrelationContainer(o11, o12, o21, o22, reads));
        adjacencyNode.add(j);
    }

    private void computeCommonReads() {
        if (commonReads == null) {
            commonReads = new long[sample.referenceLength][sample.referenceLength];
            log("Common reads matrix calculation");

            int cores = Start.threadsNumber();
            List<Callable<Boolean>> tasks = new ArrayList<>();

            for (int i = 0; i < sample.referenceLength; i++) {
                tasks.add(new CommonReadsIlluminaParallelTask(i, commonReads, struct, sample, log));
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
            consensus = Utils.consensus(profile(), al);
            log(" Reference length = "+consensus.length());
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
    private List<PairEndRead> processOverlaps(List<PairEndRead> reads) {
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
        return reads;
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
    private List<SNVResultContainer> processCliques(Set<Set<Integer>> cliques) {
        //add consensus clique
        cliques.add(new HashSet<>());
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(struct, cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, consensus())));
        log("Start build clusters");
        Map<String, List<PairEndRead>> clusters = buildClusters(allPositionsInCliques, allCliquesCharacters);
        log(" - DONE");
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        List<SNVResultContainer> haplotypes = clusters.entrySet().parallelStream().filter(s -> s.getValue().size() > 10).map(s -> {
            List<PairEndRead> cluster = s.getValue();
            IlluminaSNVSample snvSample = new IlluminaSNVSample("tmp", new ArrayList<>(cluster), sample.referenceLength);
            //if there is no reads, put consensus there
            double[][] haploProfile = Utils.profile(snvSample, al);
            for (int i = 0; i < haploProfile[0].length; i++) {
                double max = 0;
                for (double[] aProfile : haploProfile) {
                    if (aProfile[i] > max) {
                        max = aProfile[i];
                    }
                }
                if (!(max > 0.01)) {
                    int i1 = al.indexOf(consensus().charAt(i)) == -1 ? '-' : al.indexOf(consensus().charAt(i));
                    haploProfile[i1][i] = 1;
                }
            }
            String haplotype = Utils.consensus(haploProfile, al);
            //rare case where all reads have N on a certain position
            StringBuilder str = new StringBuilder(haplotype);
            for (int i = 0; i < str.length(); i++) {
                if (str.charAt(i) == 'N') {
                    str.setCharAt(i, consensus().charAt(i));
                }
            }
            haplotype = str.toString();
            Set<Integer> snps = new HashSet<>();
            for (int i = 0; i < haplotype.length(); i++) {
                if (haplotype.charAt(i) != consensus().charAt(i)) {
                    snps.add(splittedPosition(i, haplotype.charAt(i)));
                }
            }
            Clique haplotypeClique = new Clique(snps, consensus());
            SNVResultContainer container = new SNVResultContainer(s.getKey(), haplotypeClique, haplotype, cluster);
            container.sourceClique = getSourceClique(allPositionsInCliques, cliquesSet, s.getKey(), container);
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
        List<Double> frequencies = new IlluminaEM().frequencies(h, sample);
        // there is no evidence that consensus haplotype exists if it's frequency is less than 25%
        // remove it and recalculate frequencies
        for (int i = 0; i < result.size(); i++) {
            if (result.get(i).haploClique.splittedSnps.isEmpty() && frequencies.get(i) < 0.25) {
                result.remove(i);
                h.remove(i);
                log("Consensus haplotype was removed");
                frequencies = new IlluminaEM().frequencies(h, sample);
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
     * Build read clusters based on cliques. Each read will go to nearest clique in terms of Hamming distance
     *
     * @param allPositionsInCliques sorted array or all positions with at least one clique
     * @param allCliquesCharacters  characters in cliques according to allPositionsInCliques.
     *                              Has consensus allele if clique doesn't include particular position from allPositionsInCliques
     * @return Map with clusters, where key is string of clique characters, value is a set of reads
     */
    private Map<String, List<PairEndRead>> buildClusters
    (List<Integer> allPositionsInCliques, List<String> allCliquesCharacters) {
        Map<String, List<PairEndRead>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new ArrayList<>()));
        if (allCliquesCharacters.size() == 1) {
            clusters.get(allCliquesCharacters.get(0)).addAll(sample.reads);
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
                for (int i = 0; i < allCliquesCharacters.size(); i++) {
                    int d = 0;
                    int coincidences = 0;
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
                    if (coincidences != d && d < min) {
                        min = d;
                        minI = new ArrayList<>();
                    }
                    if (d == min) {
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
                    minI.forEach(e -> clusters.get(allCliquesCharacters.get(e)).add(read));
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

    private String getCorrelationKey(int first, int second, char m1, char m2) {
        return first + m1 + "_" + second + m2;
    }
}
