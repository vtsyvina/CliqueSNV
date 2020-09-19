package edu.gsu.algorithm;

import edu.gsu.algorithm.util.AlgorithmUtils;
import edu.gsu.algorithm.util.NotEdgeParallelTask;
import edu.gsu.algorithm.util.TrueFrequencyEstimator;
import edu.gsu.model.Clique;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.SNVStructure;
import edu.gsu.start.Start;
import edu.gsu.util.HammingDistance;
import edu.gsu.util.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class AbstractSNV {
    public int CLOSE_POSITIONS_THRESHOLD = 5;
    public static final double NO_EDGE_THRESHOLD = 0.1;
    public static final int MINIMUM_READS_FOR_NO_EDGE = 50;
    protected final double HAPLOTYPE_CUT_THRESHOLD;
    //    public static final double GREY_ZONE_NO_EDGE_THROSHOLD = 0.1;
    public static String al = "ACGT-N";
    public static int minorCount = al.length() - 2;
    public boolean log = false;
    public int MIN_O22_THRESHOLD;
    private double MIN_O22_FREQ;
    Map<String, long[]> osCache = new ConcurrentHashMap<>();
    List<Clique> answer = new ArrayList<>();
    public double MAX_READ_ERROR = 0.1;
    public double LOG_THRESHOLD;
    public boolean USE_LOG_PVALUE;

    AbstractSNV(int minThreshold, double minFreq) {
        this.MIN_O22_THRESHOLD = minThreshold;
        this.MIN_O22_FREQ = minFreq;
        this.LOG_THRESHOLD = Double.parseDouble(Start.settings.getOrDefault("-lt", "200"));
        this.USE_LOG_PVALUE = Start.settings.get("-lp") != null;
        this.HAPLOTYPE_CUT_THRESHOLD = Double.parseDouble(Start.settings.getOrDefault("-tf", "0.05"));
    }

    /**
     * Method to build all cliques based on adjacencyList for SNPs.
     * Also it merges cliques if they have more than 50% of edges in common into 'pseudoclique'
     *
     * @param adjacencyList Matrix with edges (in splitted sample)
     * @return set of cliques (clique - set of splitted positions that encode position + minor)
     */
    Set<Set<Integer>> getMergedCliques(List<Set<Integer>> adjacencyList) {
        if (Start.settings.getOrDefault("-cm", "accurate").equals("transitive")) {
            return newCliquesFinding(adjacencyList);
        }
        if (Start.settings.getOrDefault("-cm", "accurate").equals("accurate")) {
            return slowCliquesMerging(adjacencyList);
        }
        return fastCliquesMerging(adjacencyList);
    }

    public void log(String s) {
        if (log) System.out.println(s);
    }

    void logSame(String s) {
        if (log) System.out.print(s);
    }

    /**
     * Compute amount of common edges between the two given cliques. Returns 0 if cliques have common position but different minors
     */
    private int mergeScore(Set<Integer> c1, Set<Integer> c2, List<Set<Integer>> adjacencyList) {
        int matches = 0;
        for (Integer i : c1) {
            for (Integer i2 : c2) {
                if (i / minorCount == i2 / minorCount && i % minorCount != i2 % minorCount) {
                    return 0;
                }
                if (adjacencyList.get(i).contains(i2)) {
                    matches++;
                }
            }
        }
        return matches;
    }

    /**
     * Gives corresponding minot for given splitted snp
     */
    char minor(int splittedSnp) {
        int allele1 = splittedSnp % minorCount >= Utils.getMajorAllele(consensus(), al, splittedSnp / minorCount) ?
                splittedSnp % minorCount + 1 :
                splittedSnp % minorCount;
        return al.charAt(allele1);
    }

    /**
     * Transforms position + minor into splitted SNP position
     */
    int splittedPosition(int pos, char minor) {
        int r = pos * minorCount;
        int allele = al.indexOf(minor) >= Utils.getMajorAllele(consensus(), al, pos) ? al.indexOf(minor) - 1 : al.indexOf(minor);
        return r + allele;
    }

    /**
     * Calculates hits for given row of minors.
     * Gets all minors in column for a certain position for minorRow and increments appropriate hits position
     *
     * @param struct          SNV data structure
     * @param rowMinor        array for all minors for some position
     * @param referenceLength Reference length
     * @return array with all minor hits for given row
     */
    public int[] getHits(SNVStructure struct, int[] rowMinor, int referenceLength) {
        int[] hits = new int[referenceLength];

        for (int j = 0; j < rowMinor.length; j++) {
            int[] column = struct.colMinors[rowMinor[j]];
            for (int aColumn : column) {
                hits[aColumn]++;
            }
        }
        return hits;
    }

    /**
     * Transform each clique into a string of equal length. The string corresponds to clique representation in all positions
     * that are covered by any clique. If clique doesn't have some position in it, than consensus allele will be in this position
     *
     * @param struct                SNV structure
     * @param cliques               Set of cliques in form of splitted SNPs(where position and minor are encoded in a single number)
     * @param allPositionsInCliques List of all positions covered by any clique
     * @return list with strings representing given cliques
     */
    List<String> getAllCliquesCharacters(SNVStructure struct, Set<Set<Integer>> cliques, List<Integer> allPositionsInCliques) {
        List<String> allCliquesCharacters = new ArrayList<>();
        StringBuilder str = new StringBuilder();
        allPositionsInCliques.forEach(i -> str.append(consensus().charAt(i)));
        Set<Clique> clicuesSet = new HashSet<>();
        cliques.forEach(c -> clicuesSet.add(new Clique(c, consensus())));
        clicuesSet.forEach(c -> {
            StringBuilder tmp = new StringBuilder(str);
            for (int i = 0; i < c.snps.size(); i++) {
                int snpPosition = allPositionsInCliques.indexOf(c.snps.get(i));
                tmp.replace(snpPosition, snpPosition + 1, String.valueOf(c.minors.charAt(i)));
            }
            allCliquesCharacters.add(tmp.toString());
        });
        return allCliquesCharacters;
    }

    /*
     *   if for any 2 positions we have edges like  X <-> Y and X <-> Z,
     *   then we delete edge with less frequency of second allele(to avoid false positive cliques)
     */
    void removeEdgesForSecondMinors(List<Set<Integer>> adjacencyMatrix, SNVStructure struct) {
        for (int i = 0; i < adjacencyMatrix.size(); i++) {

            //convert adjacencyMatrix for each position into map (position -> all correlated alleles in this position)
            Map<Integer, Set<Integer>> columnsEdges = new HashMap<>();
            adjacencyMatrix.get(i).forEach(j -> {
                if (!columnsEdges.containsKey(j / minorCount)) {
                    columnsEdges.put(j / minorCount, new HashSet<>());
                }
                columnsEdges.get(j / minorCount).add(j);
            });

            int finalI1 = i;
            // for every position where we have more than 1 edge
            columnsEdges.entrySet().stream().filter(e -> e.getValue().size() > 1).forEach(e -> {
                final double[] max = {0};
                final int[] maxM = {0};
                e.getValue().forEach(m -> {
                    if (struct.profile[al.indexOf(minor(m))][m / minorCount] > max[0]) {
                        max[0] = struct.profile[al.indexOf(minor(m))][m / minorCount];
                        maxM[0] = m;
                    }
                });
                //for illumina we expect false edges only if frequent minor is deletion
                if (technology() == Technology.ILLUMINA && minor(maxM[0]) != '-') {
                    return;
                }
                //remove all non-maximum frequency edges
                e.getValue().forEach(m -> {
                    if (m != maxM[0]) {
                        adjacencyMatrix.get(finalI1).remove(m);
                        adjacencyMatrix.get(m).remove(finalI1);
                        log("remove " + finalI1 / minorCount + " " + m / minorCount + " " + minor(finalI1) + " " + minor(m));
                    }
                });
            });
        }
    }

    public int getAllele(int i) {
        return i % minorCount >= Utils.getMajorAllele(consensus(), al, i / minorCount) ? i % minorCount + 1 : i % minorCount;
    }

    //TODO refactor this
    Clique getSourceClique(List<Integer> allPositionsInCliques, Set<Clique> cliquesSet, String cliqueString, SNVResultContainer container) {
        Clique sourceClique = container.sourceClique;
        for (Clique clique : cliquesSet) {
            boolean fl = true;
            for (int i = 0; i < clique.snps.size(); i++) {
                if (cliqueString.charAt(allPositionsInCliques.indexOf(clique.snps.get(i))) != clique.minors.charAt(i)) {
                    fl = false;
                    break;
                }
            }
            if (fl) {
                sourceClique = clique;
                //small hack. If clique is empty than fl will be true for any source clique
                if (clique.minors.length() > 0) {
                    break;
                }
            }
        }
        return sourceClique;
    }

    /**
     * Creates clique with all given SNPs, estimates frequency based on given transition probability
     *
     * @param snps given SNPs
     * @param eps  transition probability
     * @return clique with all SNPs and frequencies
     * @see TrueFrequencyEstimator#estimateFrequency(double[], int, int, double)
     */
    Clique getVCCliqie(Set<Integer> snps, double eps) {
        List<Double> frequencies = new ArrayList<>();
        snps.stream().sorted().forEach(snp -> {
            double[] f = new double[minorCount + 1];
            double sum = 0;
            for (int i = 0; i < minorCount + 1; i++) {
                f[i] = profile()[i][snp / minorCount];
                sum += f[i];
            }
            for (int i = 0; i < minorCount + 1; i++) {
                f[i] /= sum;
            }
            long count = snps.stream().filter(e -> e / minorCount == snp / minorCount).count();
            double[] t = TrueFrequencyEstimator.estimateFrequency(f, (int) count + 1, getAllele(snp), eps);
            frequencies.add(t[getAllele(snp)]);
        });
        return new Clique(snps, consensus(), frequencies);
    }

    /**
     * Returns consensus for given input
     *
     * @return consensus
     */
    public abstract String consensus();

    /**
     * Returns profile for given input
     *
     * @return consensus
     */
    public abstract double[][] profile();

    /**
     * Specify reads' technology
     *
     * @return
     */
    protected abstract Technology technology();

    protected abstract long[] getOsImp(int i, int j);

    /**
     * Return number of common reads between two splitted positions
     */
    public abstract long getCommonReadsCount(int i, int j);

    long[] getOs(int i, int j) {
        long[] os = osCache.get(i + "_" + j);
        if (os == null) {
            os = getOsImp(i, j);
            osCache.putIfAbsent(i + "_" + j, os);
        }
        return os;
    }

    /**
     * Returns p-value between splitted allele i and splitted allele j
     * Returns 0 if O22 / (MAX_READ_ERROR * O11 + MAX_READ_ERROR * (O12 + O21)) > MAX_READ_ERROR - this cannot be explained by read errors
     *
     * @param i splitted allele
     * @param j splitted allele
     * @return
     */
    public abstract double getP(int i, int j);

    public boolean hitsFitThreshold(int hits, long coverage) {
        return hits >= MIN_O22_THRESHOLD && ((double) hits) / coverage > MIN_O22_FREQ;
    }

    /**
     * Checks if for given Oij there exist 22 edge (correlation between minors)
     *
     * @param o11
     * @param o12
     * @param o21
     * @param o22
     * @param reads
     * @param adjustment
     * @param referenceLength
     * @return true if minors have correlation, false otherwise
     */
    public boolean hasO22Edge(long o11, long o12, long o21, long o22, long reads, double adjustment, int referenceLength) {
        double p = (o12 * o21) / ((double) o11 * reads);
        if (p > 1) {
            return false;
        }
        if (p < 1E-12) {
            return true;
        } else {
            double pvalue = USE_LOG_PVALUE ? Utils.binomialLogPvalue((int) o22, p, (int) reads) : Utils.binomialPvalue((int) o22, p, (int) reads);
            boolean fl = USE_LOG_PVALUE ? pvalue > LOG_THRESHOLD : pvalue < adjustment / (referenceLength * 1. * (referenceLength - 1) / 2);
            return fl
                    || o22 / (MAX_READ_ERROR * o11 + (1 - MAX_READ_ERROR) * (o12 + o21)) > MAX_READ_ERROR;
        }
    }

    void readAnswerHaplotypes() {
        if (Start.answer != null) {
            answer = new ArrayList<>();
            for (String sequence : Start.answer.reads) {
                Set<Integer> tmp = new HashSet<>();
                for (int i = 0; i < consensus().length(); i++) {
                    if (sequence.charAt(i)
                            != consensus().charAt(i)) {
                        tmp.add(splittedPosition(i, sequence.charAt(i)));
                    }

                }
                answer.add(new Clique(tmp, consensus()));
            }
            log("Answer:");
            if (log) answer.forEach(System.out::println);
            log("");
        }
    }

    public void outputAnswerChecking(List<SNVResultContainer> snvResultContainers) {
        if (Start.answer != null && !snvResultContainers.isEmpty()) {
            try {
                for (int i = 0; i < Start.answer.reads.length; i++) {
                    String ans = Start.answer.reads[i].substring(0, snvResultContainers.get(0).haplotype.length());
                    HammingDistance hd = new HammingDistance();
                    if (snvResultContainers.stream().anyMatch(s -> hd.apply(s.haplotype, ans) == 0)) {
                        log("Found haplotype for " + answer.get(i));
                    } else {
                        log("Failed to find haplotype for " + answer.get(i));
                        int min = ans.length();
                        int minI = 0;
                        for (int j = 0; j < snvResultContainers.size(); j++) {
                            String haplotype = snvResultContainers.get(j).haplotype;
                            if (hd.apply(ans, haplotype) < min) {
                                min = hd.apply(ans, haplotype);
                                minI = j;
                            }
                        }
                        log(" closest is within " + min + " distance " + snvResultContainers.get(minI).haploClique);
                    }
                }
            } catch (IllegalArgumentException e) {
                log("Error occured while comparing answer with obtained haplotypes");
                // okaaaaaaay. Stuff happens
            }
        }
    }

    /**
     * Will find all possible merged cliques based on edges between cliques and 'no edges'
     */
    private Set<Set<Integer>> slowCliquesMerging(List<Set<Integer>> adjacencyList) {
        long st = System.currentTimeMillis();
//        List<Set<Integer>> cliques = new ArrayList<>(newCliquesFinding(adjacencyList));
        List<Set<Integer>> cliques = AlgorithmUtils.findCliquesIgnoreIsolated(adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());

        if (log && cliques.size() < 150 && !cliques.isEmpty())
            cliques.stream().map(c -> new Clique(c, consensus())).forEach(System.out::println);
        if (cliques.isEmpty()) {
            return new HashSet<>();
        }
        log("Cliques time " + (System.currentTimeMillis() - st));
        log("Cliques before merge " + cliques.size());
        double average = cliques.stream().mapToInt(Set::size).average().getAsDouble();
        log("Average clique size "+ average);
        //average < 2.001 && cliques.size() > 500
        if(cliques.size() > 500){
            log("Too much cliques to process");
            return new HashSet<>();
        }
        // adjacency list for relations between cliques - at least one edge between them and no 'no edges'
        List<Set<Integer>> realCliqueEdges = new ArrayList<>();
        for (int i = 0; i < cliques.size(); i++) {
            realCliqueEdges.add(new HashSet<>());
        }
        log("Calculating read cliques graph");
        List<Callable<List<Integer>>> tasks = new ArrayList<>();
        for (int i = 0; i < cliques.size(); i++) {
            int finalI = i;
            tasks.add(() -> {
                List<Integer> result = new ArrayList<>();
                for (int j = finalI + 1; j < cliques.size(); j++) {
                    if (mergeScore(cliques.get(finalI), cliques.get(j), adjacencyList) > 0) {
                        result.add(j);
                    }
                }
                return result;
            });
        }
        try {
            List<Future<List<Integer>>> futures = Start.executor.invokeAll(tasks);
            for (int i = 0; i < cliques.size(); i++) {
                List<Integer> edges = futures.get(i).get();
                int finalI = i;
                edges.forEach(e -> {
                    realCliqueEdges.get(finalI).add(e);
                    realCliqueEdges.get(e).add(finalI);
                });
            }
        } catch (InterruptedException | ExecutionException e) {
            e.printStackTrace();
        }
        log("");
        log("Calculating full graph");
        //full graph between cliques minus all 'no edges'
        List<Set<Integer>> cAdList = new ArrayList<>();
        for (int i = 0; i < cliques.size(); i++) {
            cAdList.add(new HashSet<>());
            for (int j = 0; j < cliques.size(); j++) {
                if (i == j) continue;
                cAdList.get(i).add(j);
            }
        }
        computeNoEdges(adjacencyList, cliques, realCliqueEdges, cAdList);
        if (!cAdList.isEmpty()) {
            double degree = cAdList.stream().mapToInt(Set::size).average().getAsDouble();
            log("Average clique degree " + degree + " start compute cliques");
//            if (degree > 200){
//                return new HashSet<>();
//            }
        }
        st = System.currentTimeMillis();
        Set<Set<Integer>> cliqueCliques = AlgorithmUtils.findCliques(cAdList);
        log("Cliques cliques time " + (System.currentTimeMillis() - st));

        int[] components = AlgorithmUtils.connectedComponents(realCliqueEdges);
        if (components.length < 200) log("Components " + Arrays.toString(components));

        Set<Set<Integer>> mergedCliques = new HashSet<>();
        log("Cliques of cliques size " + cliqueCliques.size());
//        if (cliqueCliques.size() > 500){
//            return new HashSet<>();
//        }
        boolean fl = cliqueCliques.size() < 100;
        //merge cliques but if they are from different components then merge only for components (if A,B from 1, and C,D from 2. Then it will merge only A with B and C with D)
        for (Set<Integer> setOfCliquesToMerge : cliqueCliques) {
            if (fl) log(setOfCliquesToMerge.toString());
            Map<Integer, Set<Integer>> cliquesByComponents = new HashMap<>();
            for (Integer c : setOfCliquesToMerge) {
                cliquesByComponents.putIfAbsent(components[c], new HashSet<>());
                cliquesByComponents.get(components[c]).add(c);
            }
            for (Set<Integer> clcl : cliquesByComponents.values()) {
                Set<Integer> newC = new HashSet<>();
                for (Integer cl : clcl) {
                    newC.addAll(cliques.get(cl));
                }
                mergedCliques.add(newC);
            }
        }
        log("");
        log("Finished with merging. Start filtering");
        //some cliques after merge (especially when we split by components) may entirely lie in other cliques. We remove them
        Set<Set<Integer>> toRemove = new HashSet<>();
        for (Set<Integer> clique : mergedCliques) {
            for (Set<Integer> clique2 : mergedCliques) {
                if (clique != clique2 && clique2.containsAll(clique)) {
                    toRemove.add(clique);
                }
            }
        }
        log("");
        mergedCliques.removeAll(toRemove);
        if (mergedCliques.size() > 300){
            log("Too much merged cliques "+mergedCliques.size());
            return  new HashSet<>();
        }
        return mergedCliques;
    }

    /**
     * Compute 'no edges' in parallel
     *
     * @param adjacencyList   adjacency list for snp alleles
     * @param cliques         list of all cliques in snps
     * @param realCliqueEdges adjacency list of relations between cliques
     * @param cAdList         adjusted list of relations (full - no edges)
     */
    private void computeNoEdges(List<Set<Integer>> adjacencyList, List<Set<Integer>> cliques, List<Set<Integer>> realCliqueEdges, List<Set<Integer>> cAdList) {
        //splitting all pairs of cliques into batches to run in parallel. Each thread will take equivalent number of cliques park
        List<List<Integer>> first = new ArrayList<>();
        List<List<Integer>> second = new ArrayList<>();
        int cores = Start.threadsNumber();
        for (int i = 0; i < cores; i++) {
            first.add(new ArrayList<>());
            second.add(new ArrayList<>());
        }
        int t = 0;
        for (int i = 0; i < cliques.size(); i++) {
            for (int j = i + 1; j < cliques.size(); j++) {
                first.get(t % cores).add(i);
                second.get(t % cores).add(j);
                t++;
            }
        }
        List<Callable<List<String>>> tasks = new ArrayList<>();
        log("Total no edge tasks " + t);
        for (int i = 0; i < cores; i++) {
            tasks.add(new NotEdgeParallelTask(this, cliques, adjacencyList, first.get(i), second.get(i)));
        }
        log("Start finding not edges");
        try {
            List<Future<List<String>>> futures = Start.executor.invokeAll(tasks);
            futures.forEach(future -> {
                try {
                    List<String> strings = future.get();
                    for (String s : strings) {
                        String[] split = s.split("_"); // I'm just too lazy to do it nice
                        int i = Integer.parseInt(split[0]);
                        int j = Integer.parseInt(split[1]);
                        //remove no edges from both adjacency lists
                        cAdList.get(i).remove(j);
                        cAdList.get(j).remove(i);
                        realCliqueEdges.get(i).remove(j);
                        realCliqueEdges.get(j).remove(i);
                    }
                } catch (InterruptedException | ExecutionException e) {
                    System.err.println("Error! Parallel tasks were not successful on get");
                    e.printStackTrace();
                }
            });
        } catch (InterruptedException e) {
            System.err.println("Error! Parallel tasks were not successful on invoke");
            e.printStackTrace();
        }
        log(" - DONE");
    }

    /**
     * Merge just best matching cliques. So one clique will match only one haplotype
     */
    private Set<Set<Integer>> fastCliquesMerging(List<Set<Integer>> adjacencyList) {
        long st = System.currentTimeMillis();
        Set<Set<Integer>> cliques = AlgorithmUtils.findCliquesIgnoreIsolated(adjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toSet());
        log("Cliques time " + (System.currentTimeMillis() - st));
        log("Cliques before merge " + cliques.size());
        if (log && !cliques.isEmpty())
            cliques.stream().map(c -> new Clique(c, consensus())).forEach(System.out::println);
        Set<Set<Integer>> mergedCliques = new HashSet<>();
        Set<Set<Integer>> alreadyMerged = new HashSet<>();
        boolean fl = true;
        while (fl) {
            fl = false;
            for (Set<Integer> clique : cliques) {
                int bestScore = 0;
                Set<Integer> bestClique = null;
                if (alreadyMerged.contains(clique)) {
                    continue;
                }
                //find best merging score for current clique
                for (Set<Integer> clique2 : cliques) {
                    if (clique == clique2 || alreadyMerged.contains(clique2)) {
                        continue;
                    }
                    int score = mergeScore(clique, clique2, adjacencyList);
                    if (bestClique == null || (double) score / clique2.size() > (double) bestScore / bestClique.size()) {
                        //do not merge if we are sure that there is 'no edge' for sure
                        if (noEdge(clique, bestClique, adjacencyList)) {
                            continue;
                        }
                        bestScore = score;
                        bestClique = clique2;
                    }
                }
                // if we found two cliques without 'no edge' and they have something in common
                if (bestScore > 0) {
                    HashSet<Integer> newClique = new HashSet<>(clique);
                    log("Merge " + new Clique(clique, consensus()) + " " + new Clique(bestClique, consensus()) + " " + cliques.size());
                    newClique.addAll(bestClique);
                    mergedCliques.add(newClique);
                    alreadyMerged.add(bestClique);
                    alreadyMerged.add(clique);
                    fl = true;
                } else {
                    mergedCliques.add(clique);
                    alreadyMerged.add(clique);
                }
            }
            cliques = mergedCliques;
            mergedCliques = new HashSet<>();
            alreadyMerged = new HashSet<>();
        }
        return cliques;
    }

    private Set<Set<Integer>> newCliquesFinding(List<Set<Integer>> adjacencyList) {
        List<Set<Integer>> newAList = new ArrayList<>(adjacencyList.size());
        for (Set<Integer> integers : adjacencyList) {
            newAList.add(new HashSet<>(integers));
        }
        int[] components = AlgorithmUtils.connectedComponents(newAList);
        for (int i = 0; i < newAList.size(); i++) {
            for (int j = 0; j < components.length; j++) {
                if (i != j && components[j] == components[i]) {
                    newAList.get(i).add(j);
                }
            }
        }
        int[] nonIsolated = IntStream.range(0, newAList.size()).filter(i -> !newAList.get(i).isEmpty()).toArray();
        for (int i = 0; i < nonIsolated.length; i++) {
            for (int j = i + 1; j < nonIsolated.length; j++) {
                int first = nonIsolated[i];
                int second = nonIsolated[j];
                if (getP(first, second) > 0.01 && Math.abs(first - second) / minorCount > 5) {
                    newAList.get(first).remove(second);
                    newAList.get(second).remove(first);
                    //remove transitive edges from neighborhood or first vertex that are connected to second vertex
                    Iterator<Integer> transIt = newAList.get(first).iterator();
                    while (transIt.hasNext()) {
                        Integer next = transIt.next();
                        if (!adjacencyList.get(first).contains(next) && adjacencyList.get(second).contains(next) && getP(first, next) > 1e-6) {
                            transIt.remove();
                        }
                    }
                    //same for second vertex
                    transIt = newAList.get(second).iterator();
                    while (transIt.hasNext()) {
                        Integer next = transIt.next();
                        if (!adjacencyList.get(second).contains(next) && adjacencyList.get(first).contains(next) && getP(second, next) > 1e-6) {
                            transIt.remove();
                        }
                    }
                    // remove transitive edges between neighbors of first and second
                    for (Integer firstN : adjacencyList.get(first)) {
                        for (Integer secondN : adjacencyList.get(second)) {
                            if (newAList.get(firstN).contains(secondN) && !adjacencyList.get(firstN).contains(secondN) && getP(firstN, secondN) > 5e-5) {
                                newAList.get(firstN).remove(secondN);
                                newAList.get(secondN).remove(firstN);
                            }
                        }
                    }
                }
            }
        }
        return AlgorithmUtils.findCliques(newAList).stream().filter(c -> c.size() > 1).collect(Collectors.toSet());
    }

    /**
     * Determines if there is no edge between two cliques
     *
     * @param c1 first clique
     * @param c2 second clique
     * @return true if there is 'no edge', false - otherwise
     */
    public boolean noEdge(Set<Integer> c1, Set<Integer> c2, List<Set<Integer>> adjacencyList) {
        boolean wasAlready = false;
        if (c1 == null || c2 == null) {
            return false;
        }
        for (Integer i : c1) {
            for (Integer j : c2) {
                if (i / minorCount == j / minorCount && !i.equals(j)) {
                    return true;
                }
                if (Math.abs(i - j) < CLOSE_POSITIONS_THRESHOLD * minorCount || adjacencyList.get(i).contains(j) || getCommonReadsCount(i, j) < MINIMUM_READS_FOR_NO_EDGE) {
                    continue;
                }

                double p = getP(i, j);
                if (p > NO_EDGE_THRESHOLD) { // || p > GREY_ZONE_NO_EDGE_THROSHOLD && wasAlready) {
                    return true;
                }
//                if (p > GREY_ZONE_NO_EDGE_THROSHOLD) {
//                    wasAlready = true;
//                }
            }
        }
        return false;
    }

    public boolean isLog() {
        return log;
    }

    public void setLog(boolean log) {
        this.log = log;
    }

    enum Technology {
        PACBIO, ILLUMINA
    }
}
