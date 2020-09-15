package edu.gsu.algorithm.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;

public class AlgorithmUtils {

    /**
     * Binary search return index of x in sorted array or index of first element that is lower than x.
     * Returns 0 if x is less then arr[0]
     *
     * @param arr Sorted array
     * @param x   Value to find
     * @return index of x
     */
    public static int binarySearch(int arr[], int x) {
        int l = 0, r = arr.length - 1;
        while (l <= r) {
            int m = l + (r - l) / 2;

            // Check if x is present at mid
            if (arr[m] == x)
                return m;

            // If x greater, ignore left half
            if (arr[m] < x)
                l = m + 1;

                // If x is smaller, ignore right half
            else
                r = m - 1;
        }

        // if we reach here, then element was not present
        return r == -1 ? 0 : r;
    }

    /**
     * Simple implementation of Bron - Kerbosch algorithm for finding all maximum cliques (psedocode is on Wiki)
     */
    public static Set<Set<Integer>> findCliques(List<Set<Integer>> adjacencyList) {
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < adjacencyList.size(); i++) {
            p.add(i);
        }
        return findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyList);
    }

    public static Set<Set<Integer>> findCliquesIgnoreIsolated(List<Set<Integer>> adjacencyList) {
        // put only non isolated vertices, since cliques of size 1 are not interested for us
        // significantly reduce recursive calls since the graph is way smaller
        List<Set<Integer>> newAdjacencyList = new ArrayList<>();
        Map<Integer, Integer> newToOldMap = new HashMap<>();
        Map<Integer, Integer> oldToNewMap = new HashMap<>();
        int nonIsolatedCount = 0;
        for (int i = 0; i < adjacencyList.size(); i++) {
            if (adjacencyList.get(i).size() > 0) {
                newToOldMap.put(nonIsolatedCount, i);
                oldToNewMap.put(i, nonIsolatedCount);
                nonIsolatedCount++;
                newAdjacencyList.add(new HashSet<>());
            }
        }
        for (int i = 0; i < nonIsolatedCount; i++) {
            int realV = newToOldMap.get(i);
            for (Integer realU : adjacencyList.get(realV)) {
                newAdjacencyList.get(i).add(oldToNewMap.get(realU));
            }
        }
        Set<Integer> p = new HashSet<>();
        for (int i = 0; i < newAdjacencyList.size(); i++) {
            p.add(i);
        }
        Set<Set<Integer>> cliques = findCliques(new HashSet<>(), p, new HashSet<>(), newAdjacencyList);
        // revert everything back
        Set<Set<Integer>> result = new HashSet<>();
        for (Set<Integer> clique : cliques) {
            HashSet<Integer> mappedClique = new HashSet<>();
            for (Integer newU : clique) {
                mappedClique.add(newToOldMap.get(newU));
            }
            result.add(mappedClique);

        }
        return result;
    }

    private static Set<Set<Integer>> findCliques(Set<Integer> clique, Set<Integer> p, Set<Integer> x, List<Set<Integer>> adjacencyList) {
        Set<Set<Integer>> result = new HashSet<>();
        if (p.isEmpty() && x.isEmpty()) {
            result.add(clique);
        }
        Set<Integer> intersectionPX = new HashSet<>(p);
        intersectionPX.addAll(x);
        int u = -1;
        int max = 0;
        for (Integer c : intersectionPX) {
            int m = 0;
            for (Integer v : adjacencyList.get(c)) {
                if (intersectionPX.contains(v)) m++;
            }
            if (m > max) {
                max = m;
                u = c;
            }
        }
        Set<Integer> differenceP = new HashSet<>(p);
        if (u != -1) {
            differenceP.removeAll(adjacencyList.get(u));
        }
        Iterator<Integer> it = differenceP.iterator();
        while (it.hasNext()) {
            Integer v = it.next();
            Set<Integer> newCLique = new HashSet<>(clique);
            newCLique.add(v);
            Set<Integer> pIntersection = new HashSet<>(p);
            Set<Integer> xIntersection = new HashSet<>(x);
            pIntersection.retainAll(adjacencyList.get(v));
            xIntersection.retainAll(adjacencyList.get(v));
            result.addAll(findCliques(newCLique, pIntersection, xIntersection, adjacencyList));
            p.remove(v);
            x.add(v);
        }
        return result;
    }

    public static int[] connectedComponents(List<Set<Integer>> adjacencyList) {
        if (adjacencyList.size() == 0) {
            return new int[0];
        }
        int[] result = new int[adjacencyList.size()];
        int count = 0;
        int currentComponent = 1;
        Queue<Integer> q = new LinkedList<>();
        q.add(0);
        result[0] = 1;
        while (count < result.length) {
            if (q.isEmpty()) {
                for (int i = 0; i < result.length; i++) {
                    if (result[i] == 0) {
                        q.add(i);
                        currentComponent++;
                        result[i] = currentComponent;
                        break;
                    }
                }
            }
            Integer poll = q.poll();
            count++;
            for (Integer v : adjacencyList.get(poll)) {
                if (result[v] == 0) {
                    result[v] = result[poll];
                    q.add(v);
                }
            }
        }
        return result;
    }

    public static int getSortedArraysIntersectionCount(int[] x, int[] y) {
        int result = 0;
        if (x.length == 0 || y.length == 0) {
            return 0;
        }
        int si = AlgorithmUtils.binarySearch(y, x[0]);
        if(x[x.length-1] < y[0]){
            return 0;
        }
        int fi = 0;
        int firstEnd = x.length;
        int secondEnd = y.length;
        while (fi < firstEnd && si < secondEnd) {
            if (x[fi] < y[si]) {
                fi++;
            } else if (x[fi] > y[si]) {
                si++;
            } else {
                result++;
                fi++;
                si++;
            }
        }
        return result;
    }
}
