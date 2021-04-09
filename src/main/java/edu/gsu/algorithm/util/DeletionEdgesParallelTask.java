package edu.gsu.algorithm.util;

import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.Interval;
import edu.gsu.model.PairEndRead;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;

public class DeletionEdgesParallelTask implements Callable<Map<Interval, Integer>> {
    int taskNumber;
    int totalThreads;
    IlluminaSNVSample sample;
    int threshold;

    public DeletionEdgesParallelTask(int taskNumber, int totalThreads, IlluminaSNVSample sample, int threshold) {
        this.taskNumber = taskNumber;
        this.totalThreads = totalThreads;
        this.sample = sample;
        this.threshold = threshold;
    }

    @Override
    public Map<Interval, Integer> call() {
        int size = sample.reads.size();
        Map<Interval, Integer> result = new HashMap<>();
        for (int i = taskNumber; i < size; i += totalThreads) {
            PairEndRead read = sample.reads.get(i);
            List<Interval> readDeletions = getReadDeletions(read.l, read.lOffset);
            for (Interval readDeletion : readDeletions) {
                int count = result.getOrDefault(readDeletion, 0);
                result.put(readDeletion, count + 1);
            }
            readDeletions = getReadDeletions(read.r, read.rOffset);
            for (Interval readDeletion : readDeletions) {
                int count = result.getOrDefault(readDeletion, 0);
                result.put(readDeletion, count + 1);
            }
        }
        return result;
    }

    /**
     * Gets read deletions that are at least 3 bases long
     *
     * @param read   read string
     * @param offset read offset
     * @return a list with ranges of
     */
    private List<Interval> getReadDeletions(String read, int offset) {
        List<Interval> result = new ArrayList<>();
        int start = 0;
        boolean deletionRegion = false;
        for (int i = 0; i < read.length(); i++) {
            if (read.charAt(i) == '-' && !deletionRegion) {
                deletionRegion = true;
                start = i;
            }
            if (read.charAt(i) != '-' && deletionRegion) {
                deletionRegion = false;
                if (i - 1 - start >= 3) {
                    result.add(new Interval(offset + start, offset + i - 1));
                }
            }
        }
        if (deletionRegion && read.length() - start >= 3) {
            result.add(new Interval(offset + start, offset + read.length()));
        }
        return result;
    }
}
