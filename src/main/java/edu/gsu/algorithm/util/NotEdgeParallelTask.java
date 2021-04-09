package edu.gsu.algorithm.util;

import edu.gsu.algorithm.AbstractSNV;
import edu.gsu.model.Clique;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

/**
 * Class to compute no edges between given pairs of cliques.
 * Return list of strings because stupid Java doesn't have tuples and I'm too lazy to write something more elegant. Can you believe?
 */
public class NotEdgeParallelTask implements Callable<List<String>> {
    public static AtomicInteger FP = new AtomicInteger();
    public static AtomicInteger FN = new AtomicInteger();
    private final AbstractSNV a;
    private final List<Set<Integer>> cliquesList;
    private final List<Set<Integer>> adjacencyList;
    private final List<Integer> firstCliques;
    private final List<Integer> secondCliques;

    public NotEdgeParallelTask(AbstractSNV a, List<Set<Integer>> cliquesList, List<Set<Integer>> adjacencyList, List<Integer> firstCliques, List<Integer> secondCliques) {
        this.a = a;
        this.cliquesList = cliquesList;
        this.adjacencyList = adjacencyList;
        this.firstCliques = firstCliques;
        this.secondCliques = secondCliques;
        FP = new AtomicInteger();
        FN = new AtomicInteger();
    }

    @Override
    public List<String> call() {
        List<String> notEdge = new ArrayList<>();
        for (int i = 0; i < firstCliques.size(); i++) {
            boolean nonEdge = false;
            if (a.noEdge(cliquesList.get(firstCliques.get(i)), cliquesList.get(secondCliques.get(i)), adjacencyList)) {
                notEdge.add(firstCliques.get(i) + "_" + secondCliques.get(i));
                nonEdge = true;
            }

            if (a.answer != null) {
                Set<Integer> firstClique = cliquesList.get(firstCliques.get(i));
                Set<Integer> secondClique = cliquesList.get(secondCliques.get(i));
                int firstIdx = -1, secondIdx = -1;

                for (int j = 0; j < a.answer.size(); j++) {
                    if (a.answer.get(j).splittedSnps.containsAll(firstClique)) {
                        firstIdx = j;
                    }
                    if (a.answer.get(j).splittedSnps.containsAll(secondClique)) {
                        secondIdx = j;
                    }
                }
                if (firstIdx == -1 || secondIdx == -1) {
                    continue;
                }
                List<Integer> firstSorted = firstClique.stream().sorted().collect(Collectors.toList());
                List<Integer> secondSorted = secondClique.stream().sorted().collect(Collectors.toList());
                if (firstIdx == secondIdx && nonEdge) {
                    a.log("Wrong non edge " + new Clique(firstClique, a.consensus()) + " " + new Clique(secondClique, a.consensus()));
                    a.noEdge(firstClique, secondClique, adjacencyList, true);
                    FP.incrementAndGet();
                }

                if (!nonEdge && firstIdx != secondIdx &&
                        !(firstSorted.get(0) - secondSorted.get(secondSorted.size() - 1) > a.sampleFragmentLength + 20 // we cannot find no edge if cliques don't have close SNPs
                                || secondSorted.get(0) - firstSorted.get(firstSorted.size() - 1) > a.sampleFragmentLength + 20)) {
                    a.log("FN non edge");
                    a.noEdge(firstClique, secondClique, adjacencyList);
                    FN.incrementAndGet();
                }
            }
        }
        return notEdge;
    }
}
