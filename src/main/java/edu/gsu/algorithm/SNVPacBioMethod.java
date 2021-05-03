package edu.gsu.algorithm;

import edu.gsu.algorithm.em.PacBioEM;
import edu.gsu.algorithm.util.AlgorithmUtils;
import edu.gsu.algorithm.util.CommonReadsPacbioParallelTask;
import edu.gsu.algorithm.util.TrueFrequencyEstimator;
import edu.gsu.model.Clique;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.SNVStructure;
import edu.gsu.model.Sample;
import edu.gsu.start.Start;
import edu.gsu.util.HammingDistance;
import edu.gsu.util.Utils;
import edu.gsu.util.builders.SNVStructureBuilder;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.stream.Collectors;

public class SNVPacBioMethod extends AbstractSNV {
    private final Sample sample;
    public int MINIMUM_READS_NUMBER_IN_CLUSTER = 10;
    private SNVStructure srcStruct;
    private double[][] profile;
    private String consensus;
    private final int[][] commonReads;
    private boolean isForImputation;
    // How many reads among all we assume to have bad quality
    public double READS_PORTION_TO_FILTER = 0.1;
    public int MISTAKES_TO_FILTER_LIMIT = 3000;

    @Override
    protected Technology technology() {
        return Technology.PACBIO;
    }

    /**
     * @param sample Given input reads with N in the beginning or end of read to make them all equal size
     * @param log    if log is enabled
     */
    public SNVPacBioMethod(Sample sample, boolean log) {

        this.sample = sample;
        this.log = log;
        commonReads = new int[sample.reads[0].length()][sample.reads[0].length()];
        for (int i = 0; i < commonReads.length; i++) {
            for (int j = 0; j < commonReads.length; j++) {
                commonReads[i][j] = -1;
            }
        }
        sampleFragmentLength = sample.reads[0].length();
        CLOSE_POSITIONS_THRESHOLD = 5;
        initParameters(10, 0.05);
    }

    /**
     * @param sample       Given input reads with N in the beginning or end of read to make them all equal size
     * @param minThreshold threshold for O22 value to consider alleles as SNPs
     * @param minFreq      minimum frequency (relative to reads' coverage) for O22 value to consider alleles as SNP
     * @param log          if log is enabled
     */
    public SNVPacBioMethod(Sample sample, int minThreshold, double minFreq, boolean log) {

        this.sample = sample;
        this.log = log;
        commonReads = new int[sample.reads[0].length()][sample.reads[0].length()];
        for (int i = 0; i < commonReads.length; i++) {
            for (int j = 0; j < commonReads.length; j++) {
                commonReads[i][j] = -1;
            }
        }
        sampleFragmentLength = sample.reads[0].length();
        CLOSE_POSITIONS_THRESHOLD = 5;
        initParameters(minThreshold, minFreq);
    }

    /**
     * Main method for SNV method on PacBio reads input
     *
     * @return Return set of containers with result info, such as haplotypes them self,
     * human-friendly representation of haplotypes, clique and cluster of reads from which it was obtained
     */
    public List<SNVResultContainer> getHaplotypes() {
        log("Start 2SNV method");
        logSame("Compute profile");
        profile();
        //just for debug
        readAnswerHaplotypes();
        log(" - DONE");
        logSame("Compute split columns");
        log(" - DONE");
        logSame("Compute SNV data structure");
        srcStruct = SNVStructureBuilder.buildPacBio(sample, profile, consensus());
        log(" - DONE");
        log("Compute cliques");
        //run first time to get cliques
        List<Set<Integer>> cliques = run(srcStruct, sample);
        log(" - DONE");

        //remove bad reads( >23 mistakes outside of cliques positions)
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / 4)).distinct().sorted().collect(Collectors.toList());
        if (!isForImputation && "true".equals(Start.settings.getOrDefault("-lowQFilter", "true"))) {
            logSame("Remove reads with low quality ");
            Sample newSample = getFilteredSample(allPositionsInCliques);
            log(" - DONE");
            logSame("Compute profile");
            //after removing all bad reads, rerun whole process again
            double[][] newProfile = Utils.profile(newSample, al);
            log(" - DONE");
            logSame("Compute split columns");
            log(" - DONE");
            logSame("Compute SNV data structure");
            SNVStructure newStructure = SNVStructureBuilder.buildPacBio(newSample, newProfile, consensus());
            log(" - DONE");
            logSame("Compute cliques");
            cliques = run(newStructure, newSample);
            log(" - DONE");
            log("Found cliques: " + cliques.size() + "\n"
                    + (cliques.isEmpty() ?
                    "" :
                    cliques.stream().map(s -> new Clique(s, consensus()).toString()).reduce((t, s) -> t + "\n" + s).get()));
        }
        logSame("Start getting haplotypes");
        // divide by clusters and find haplotypes
        List<SNVResultContainer> snvResultContainers = processCliques(cliques, sample);
        log(" - DONE");
        //if (Start.settings.get("-ch").equals("true")) {
            snvResultContainers = filterHaplotypeFrequencies(snvResultContainers, HAPLOTYPE_CUT_THRESHOLD);
        //}
        if (snvResultContainers.size() == 0) {
            snvResultContainers = getDefaultHaplotype();
            Start.errorCode = 5;
            Start.errorMessage = "All haplotypes got too low freqeuncy";
        }
        outputAnswerChecking(snvResultContainers);
        return snvResultContainers;
    }

    public Clique getVC() {
        log("Start 2SNV method");
        logSame("Compute profile");
        profile();
        log(" - DONE");
        logSame("Compute split columns");
        log(" - DONE");
        logSame("Compute SNV data structure");
        srcStruct = SNVStructureBuilder.buildPacBio(sample, profile, consensus());
        log(" - DONE");
        logSame("Compute cliques");
        //run first time to get cliques
        List<Set<Integer>> cliques = run(srcStruct, sample);
        log(" - DONE");
        logSame("Remove reads with low quality ");
        //remove bad reads( >23 mistakes outside of cliques positions)
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / 4)).distinct().sorted().collect(Collectors.toList());
        Sample newSample = getFilteredSample(allPositionsInCliques);
        log(" - DONE");
        logSame("Compute profile");
        //after removing all bad reads, rerun whole process again
        double[][] newProfile = Utils.profile(newSample, al);
        log(" - DONE");
        logSame("Compute split columns");
        log(" - DONE");
        logSame("Compute SNV data structure");
        SNVStructure newStructure = SNVStructureBuilder.buildPacBio(newSample, newProfile, consensus());
        log(" - DONE");
        logSame("Compute SNPs");
        List<Set<Integer>> adjacencyList = getUnprocessedSnpsAdjecencyList(sample.reads[0].length() * minorCount, newStructure, newSample);
        removeEdgesForSecondMinors(adjacencyList, srcStruct);
        log(" - DONE");
        //get haplotypes based on pure cliques, i.e. not merged.
        cliques = AlgorithmUtils.findCliquesIgnoreIsolated(adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());
        List<SNVResultContainer> snvResultContainers = processCliques(cliques, sample);
        //IntStream.range(0, adjacencyList.size()).filter(i -> adjacencyList.get(i).size() > 0).boxed().collect(Collectors.toSet())
        return getVCCliqie(snvResultContainers.stream().map(c -> c.haploClique).flatMap(cl -> cl.splittedSnps.stream()).collect(Collectors.toSet()),
                TrueFrequencyEstimator.PACBIO_EPS);
    }

    /**
     * Filters out 10% of reads with biggest number of mistakes in positions out of SNPs
     *
     * @param allPositionsInCliques all positions with SNPs
     * @return filtered sample
     */
    private Sample getFilteredSample(List<Integer> allPositionsInCliques) {
        int[] mistakes = new int[MISTAKES_TO_FILTER_LIMIT];
        List<Integer>[] m = new ArrayList[MISTAKES_TO_FILTER_LIMIT];
        for (int i = 0; i < MISTAKES_TO_FILTER_LIMIT; i++) {
            m[i] = new ArrayList<>();
        }
        List<String> newSequences = new ArrayList<>();
        List<String> newNames = new ArrayList<>();
        String[] reads = sample.reads;
        for (int j = 0; j < reads.length; j++) {
            String sequence = reads[j];
            int apply = new HammingDistance().apply(sequence, consensus());
            for (int i = 0; i < sequence.length(); i++) {
                if (sequence.charAt(i) == 'N' && consensus().charAt(i) != 'N') {
                    apply--;
                }
            }
            for (Integer i : allPositionsInCliques) {
                if (sequence.charAt(i) != 'N' && sequence.charAt(i) != consensus().charAt(i)) {
                    apply--;
                }
            }
            if (apply >= MISTAKES_TO_FILTER_LIMIT) {
                apply = MISTAKES_TO_FILTER_LIMIT - 1;
            }
            mistakes[apply]++;
            m[apply].add(j);
        }
        outer:
        for (List<Integer> idxs : m) {
            for (Integer r : idxs) {
                if (newSequences.size() < (1 - READS_PORTION_TO_FILTER) * sample.reads.length) {
                    newSequences.add(sample.reads[r]);
                    newNames.add(sample.readNames[r]);
                } else {
                    break outer;
                }
            }
        }
        return new Sample(sample.name, newSequences.toArray(new String[0]), newNames.toArray(new String[0]));
    }

    @Override
    public double getP(int i, int j) {
        long[] os = getOs(i, j);
        if (os[0] == 0) {
            long tmp = os[0];
            os[0] = os[3];
            os[3] = tmp;
        }
        long reads = getCommonReadsCount(sample.reads.length, i, j, srcStruct.rowN[i / minorCount], srcStruct.rowN[j / minorCount]);
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
        long reads = getCommonReadsCount(sample.reads.length, i, j, srcStruct.rowN[i / minorCount], srcStruct.rowN[j / minorCount]);
        double p = NON_EDGE_T_FREQ;
        return Utils.binomialPvalue((int) os[3], p, (int) reads);
    }

    @Override
    protected long[] getOsImp(int i, int j) {
        int first = i / minorCount;
        int second = j / minorCount;
        long o22 = AlgorithmUtils.getSortedArraysIntersectionCount(srcStruct.rowMinors[i], srcStruct.rowMinors[j]);
        long o21 = srcStruct.rowMinors[i].length; //all 2*
        long o12 = srcStruct.rowMinors[j].length; //all *2
        // subtract all that a not 21 from o21
        for (int k = 0; k < srcStruct.rowMinors[i].length; k++) {
            if (sample.reads[srcStruct.rowMinors[i][k]].charAt(second) != consensus().charAt(second)) {
                o21--;
            }
        }
        // subtract all that a not 12 from o12
        for (int k = 0; k < srcStruct.rowMinors[j].length; k++) {
            if (sample.reads[srcStruct.rowMinors[j][k]].charAt(first) != consensus.charAt(first)) {
                o12--;
            }
        }
        long o11 = getO11(srcStruct, sample, first, second);
        return new long[]{o11, o12, o21, o22};
    }

    /**
     * Auxiliary method that calculates all hits(and p-value) between all alleles and create cliques according to them
     * It will merge cliques that share more than 50% of edges
     *
     * @param struct SNV structure that helps to count hits between minors {@link SNVStructureBuilder} {@link SNVStructure}
     * @param src    Source sample with given PacBio reads
     * @return Set of SNPs position for all found cliques
     */
    private List<Set<Integer>> run(SNVStructure struct, Sample src) {
        computeCommonReads();
        List<Set<Integer>> adjacencyList = getUnprocessedSnpsAdjecencyList(src.reads[0].length() * minorCount, struct, src);

        /*
         *   if for any 2 positions we have edges like  X <-> Y and X <-> Z,
         *   then we delete edge with less frequency of second allele(to avoid false positive cliques)
         */
        removeEdgesForSecondMinors(adjacencyList, struct);
        log("Edges found " + adjacencyList.stream().mapToInt(Set::size).sum());
        return getMergedCliques(adjacencyList);

    }

    private List<Set<Integer>> getUnprocessedSnpsAdjecencyList(int splittedLength, SNVStructure struct, Sample src) {
        List<Set<Integer>> adjacencyList = new ArrayList<>();
        for (int i = 0; i < splittedLength; i++) {
            adjacencyList.add(new HashSet<>());
            int l = struct.rowMinors[i].length;
            if (l < 10) {
                continue;
            }
            int first = i / minorCount;
            if (first < START_POSITION || first > END_POSITION) {
                continue;
            }
            int[] hits = getHits(struct, struct.rowMinors[i], splittedLength);
            for (int j = 0; j < hits.length; j++) {
                int second = j / minorCount;
                if (second < START_POSITION || second > END_POSITION) {
                    continue;
                }
                //skip small amount of hits
                if (hitsFitThreshold(hits[j], getCommonReadsCount(i, j)) && Math.abs(first - second) > CLOSE_POSITIONS_THRESHOLD) {
                    //get unsplitted columns, minors, o_kl
                    int allele1 = i % minorCount >= Utils.getMajorAllele(consensus(), al, first) ?
                            i % minorCount + 1 : i % minorCount;
                    int allele2 = j % minorCount >= Utils.getMajorAllele(consensus(), al, second) ?
                            j % minorCount + 1 : j % minorCount;

                    char m1 = al.charAt(allele1);
                    char m2 = al.charAt(allele2);
                    if (m1 == '-' || m2 == '-') {
                        continue;
                    }
                    /*
                     * false 1 means that in actual sample it has another minor or N in given position
                     */
                    long o22 = hits[j];
                    long o21 = struct.rowMinors[i].length; //all 2*
                    long o12 = struct.rowMinors[j].length; //all *2
                    // subtract all that a not 21 from o21
                    for (int k = 0; k < struct.rowMinors[i].length; k++) {
                        if (src.reads[struct.rowMinors[i][k]].charAt(second) != consensus().charAt(second)) {
                            o21--;
                        }
                    }
                    // subtract all that a not 12 from o12
                    for (int k = 0; k < struct.rowMinors[j].length; k++) {
                        if (src.reads[struct.rowMinors[j][k]].charAt(first) != consensus.charAt(first)) {
                            o12--;
                        }
                    }
                    long o11 = getO11(struct, src, first, second);
                    //amount of common reads for i and j column (first two summands to get amount of 1 in j(major and all minors that turned to 1)), then add 21 and 22. We need inly to subtract 1N
                    long reads = getCommonReadsCount(src.reads.length, i, j, struct.rowN[first], struct.rowN[second]);
                    //start calculate p-value, starting with p
                    double p = (o12 * o21) / ((double) o11 * reads);
                    if (hasO22Edge(o11, o12, o21, o22, reads, 0.0000001, src.reads[0].length())) {
                        log(String.format("%d %d %c %c m1=%d m2=%d hits=%d p=%.3e %d %d %d %d",
                                first, second, m1, m2, l, struct.rowMinors[j].length, hits[j],
                                Utils.binomialOneMinusPvalue((int) o22, p, (int) reads),
                                o11, o12, o21, o22));
                        adjacencyList.get(i).add(j);
                    }
                }
            }
        }
        return adjacencyList;
    }

    private long getO11(SNVStructure struct, Sample src, int first, int second) {
        long o11 = struct.majorsInRow[first]; //all 11(and some false positive 11),12,1N

        //subtract 1N from o11
        for (int k = 0; k < struct.rowN[second].length; k++) {
            if (src.reads[struct.rowN[second][k]].charAt(first) == consensus().charAt(first)) {//make sure that 1 in minColumn is true 1
                o11--;
            }
        }
        //subtract 12 and false positive 11 from o11
        for (int k = 0; k < minorCount; k++) {
            int[] rowMinor = struct.rowMinors[second * minorCount + k];
            for (int m : rowMinor) {
                if (src.reads[m].charAt(first) == consensus().charAt(first)) { //make sure that 1 in minColumn is true 1
                    o11--;
                }
            }
        }
        return o11;
    }

    private void computeCommonReads() {
        if (commonReads[0][0] == -1) {
            //commonReads = new long[sample.referenceLength][sample.referenceLength];
            log("Common reads matrix calculation");
            List<Callable<Boolean>> tasks = new ArrayList<>();

            for (int i = 0; i < commonReads.length; i++) {
                tasks.add(new CommonReadsPacbioParallelTask(i, commonReads, srcStruct, sample, START_POSITION, END_POSITION, log));
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

    private int getCommonReadsCount(int totalReads, int i, int j, int[] rowNFirst, int[] rowNSecond) {
        if (commonReads[i / minorCount][j / minorCount] == -1) {
            commonReads[i / minorCount][j / minorCount] = totalReads - rowNFirst.length - rowNSecond.length + AlgorithmUtils.getSortedArraysIntersectionCount(rowNFirst, rowNSecond);
            commonReads[j / minorCount][i / minorCount] = commonReads[i / minorCount][j / minorCount];
        }
        return commonReads[i / minorCount][j / minorCount];
    }


    @Override
    public String consensus() {
        if (consensus == null) {
            consensus = Utils.consensus(sample.reads, al);
            log(" Reference length = " + consensus.length());
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

    @Override
    public long getCommonReadsCount(int i, int j) {
        return getCommonReadsCount(sample.reads.length, i, j, srcStruct.rowN[i / minorCount], srcStruct.rowN[j / minorCount]);
    }

    /**
     * Methods process computed cliques. At first it separate reads into clusters and then compute consensus for each cluster and some additional information
     *
     * @param cliques Set of cliques in form of splitted SNPs(where position and minor are encoded in a single number)
     * @param src     Source sample
     * @return Set of containers with results. Each container contains haplotype itself and some additional helpful information
     */
    private List<SNVResultContainer> processCliques(List<Set<Integer>> cliques, Sample src) {
        cliques.add(new HashSet<>());
        List<Integer> allPositionsInCliques = cliques.stream().flatMap(s -> s.stream().map(c -> c / minorCount)).distinct().sorted().collect(Collectors.toList());
        List<String> allCliquesCharacters = getAllCliquesCharacters(cliques, allPositionsInCliques);

        Set<Clique> cliquesSet = new HashSet<>();
        cliques.forEach(c -> cliquesSet.add(new Clique(c, consensus())));
//        Map<String, Map<Integer, String>> clusters = buildClusters(src, allPositionsInCliques, allCliquesCharacters);
        Map<String, List<Pair<Integer, Integer>>> clustersRelative = buildClustersRelative(src, allPositionsInCliques, allCliquesCharacters);
        //skip clusters with less than 10 reads. Do some stuff for transforming output into human-friendly format
        List<SNVResultContainer> haplotypes = clustersRelative.entrySet().stream().filter(s -> s.getValue().size() > MINIMUM_READS_NUMBER_IN_CLUSTER).map(c -> {
            String[] reads = c.getValue().stream().map(e -> sample.reads[e.getKey()]).toArray(String[]::new);
            // the number of cliques in which each read goes
            int[] splitPortions = c.getValue().stream().mapToInt(Pair::getValue).toArray();
            Sample tmpSample = new Sample("tmp", reads);
            double[][] profile = Utils.profile(tmpSample, splitPortions, al);
            String haplotype = Utils.consensus(profile, al);
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
                    snps.add(splitPosition(i, haplotype.charAt(i)));
                }
            }
            Clique haplotypeClique = new Clique(snps, consensus());
            Map<Integer, String> pacBioCluster = new HashMap<>();
            c.getValue().forEach(e -> pacBioCluster.put(e.getKey(), sample.reads[e.getKey()]));
            SNVResultContainer container = new SNVResultContainer(c.getKey(), pacBioCluster, haplotypeClique, haplotype);
            container.sourceClique = getSourceClique(allPositionsInCliques, cliquesSet, c.getKey());
            return container;
        }).collect(Collectors.toList());
        Map<String, SNVResultContainer> merge = new HashMap<>();
        haplotypes.forEach(h -> {
            if (!merge.containsKey(h.haplotype)) {
                merge.put(h.haplotype, h);
            } else {
                merge.get(h.haplotype).pacBioCluster.putAll(h.pacBioCluster);
            }
        });
        List<SNVResultContainer> result = new ArrayList<>(merge.values());
        List<String> h = result.stream().map(s -> s.haplotype).collect(Collectors.toList());
        List<Double> frequencies = new PacBioEM(consensus().length() + 2000).frequencies(h, src);
        // there is no evidence that consensus haplotype exists if it's frequency is less than 25%
        // remove it and recalculate frequencies
        for (int i = 0; i < result.size(); i++) {
            if (result.get(i).haploClique.splittedSnps.isEmpty() && frequencies.get(i) < 0.25) {
                result.remove(i);
                h.remove(i);
                log("Consensus haplotype was removed");
                frequencies = new PacBioEM(consensus().length() + 2000).frequencies(h, sample);
                break;
            }
        }
        for (int i = 0; i < frequencies.size(); i++) {
            result.get(i).frequency = frequencies.get(i);
        }
        return result.stream().sorted((s1, s2) -> -Double.compare(s1.frequency, s2.frequency)).collect(Collectors.toList());
    }

    /**
     * Build read clusters based on cliques. Each read will go to nearest clique in terms of Hamming distance divided by the number of nearest cliques
     *
     * @param src                   Source sample
     * @param allPositionsInCliques sorted array or all positions with at least one clique
     * @param allCliquesCharacters  characters in cliques according to allPositionsInCliques.
     *                              Has consensus allele if clique doesn't include particular position from allPositionsInCliques
     * @return Map with clusters, where key is string of clique characters, value is a set of reads
     */
    private Map<String, List<Pair<Integer, Integer>>> buildClustersRelative(Sample src, List<Integer> allPositionsInCliques, List<String> allCliquesCharacters) {
        Map<String, List<Pair<Integer, Integer>>> clusters = new HashMap<>();
        allCliquesCharacters.forEach(s -> clusters.put(s, new ArrayList<>()));

        //dirty hack to find consensus clique since it should be processed separately
        String consensusClique = "";
        for (String characters : allCliquesCharacters) {
            boolean fl = true;
            for (int i = 0; i < characters.length(); i++) {
                if (characters.charAt(i) != consensus().charAt(allPositionsInCliques.get(i))) {
                    fl = false;
                    break;
                }
            }
            if (fl) {
                consensusClique = characters;
            }
        }
        // put all reads to appropriate container (Voronoi regions)
        String[] reads = src.reads;
        for (int readIndex = 0; readIndex < reads.length; readIndex++) {
            String s = reads[readIndex];
            List<Integer> minI = new ArrayList<>();
            int min = 1_000_000;
            for (int i = 0; i < allCliquesCharacters.size(); i++) {
                String c = allCliquesCharacters.get(i);
                int d = 0;
                int coincidences = 0;
                boolean isConsensusClique = c.equals(consensusClique);
                for (int j = 0; j < allPositionsInCliques.size(); j++) {
                    char charAtPosition = s.charAt(allPositionsInCliques.get(j));
                    if (charAtPosition == 'N') {
                        continue;
                    }
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
                if (coincidences != d && d < min && coincidences > 0) {
                    min = d;
                    minI = new ArrayList<>();
                }
                if (d == min && coincidences > 0) {
                    minI.add(i);
                }
            }

            if (!minI.isEmpty()) {
                int finalReadIndex = readIndex;
                List<Integer> finalMinI = minI;
                minI.forEach(i -> clusters.get(allCliquesCharacters.get(i)).add(new Pair<>(finalReadIndex, finalMinI.size())));
            }
        }
        return clusters;
    }

    public void setForImputation(boolean forImputation) {
        isForImputation = forImputation;
    }

    @Override
    protected List<PairEndRead> getIlluminaReads() {
        return null;
    }

    @Override
    protected Map<Integer, String> getPacBioCluster() {
        Map<Integer, String> map = new HashMap<>();
        for (int i = 0; i < sample.reads.length; i++) {
            map.put(i, sample.reads[i]);
        }
        return map;
    }
}
