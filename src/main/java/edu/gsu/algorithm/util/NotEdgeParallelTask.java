package edu.gsu.algorithm.util;

import edu.gsu.algorithm.AbstractSNV;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Class to compute no edges between given pairs of cliques.
 * Return list of strings because stupid Java doesn't have tuples and I'm too lazy to write something more elegant. Can you believe?
 */
public class NotEdgeParallelTask implements Callable<List<String>> {
    private static AtomicInteger counter = new AtomicInteger();
    private AbstractSNV a;
    private List<Set<Integer>> cliquesList;
    private List<Set<Integer>> adjacencyList;
    private List<Integer> firstCliques;
    private List<Integer> secondCliques;

    public NotEdgeParallelTask(AbstractSNV a, List<Set<Integer>> cliquesList, List<Set<Integer>> adjacencyList, List<Integer> firstCliques, List<Integer> secondCliques) {
        this.a = a;
        this.cliquesList = cliquesList;
        this.adjacencyList = adjacencyList;
        this.firstCliques = firstCliques;
        this.secondCliques = secondCliques;
        counter = new AtomicInteger();
    }

    @Override
    public List<String> call() {
        List<String> notEdge = new ArrayList<>();
        for (int i = 0; i < firstCliques.size(); i++) {
            if (a.noEdge(cliquesList.get(firstCliques.get(i)), cliquesList.get(secondCliques.get(i)), adjacencyList)) {
                notEdge.add(firstCliques.get(i) + "_" + secondCliques.get(i));
                continue;
            }
            if (a.isLog()) {
                if (counter.incrementAndGet() % 1000 == 0)
                    System.out.print("\r" + counter.get());
            }
        }
        return notEdge;
    }
}
