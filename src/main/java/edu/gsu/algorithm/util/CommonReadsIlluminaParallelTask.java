package edu.gsu.algorithm.util;

import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.SNVStructure;

import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

public class CommonReadsIlluminaParallelTask implements Callable<Boolean> {
    private static AtomicInteger iter = new AtomicInteger();
    private int i;
    private int[][] commonReads;
    private SNVStructure struct;
    private IlluminaSNVSample src;
    private boolean log;
    // start end end position where we are searching for SNPs
    private int start;
    private int end;


    public CommonReadsIlluminaParallelTask(int i, int[][] commonReads, SNVStructure struct, IlluminaSNVSample src, int start, int end, boolean log) {
        this.i = i;
        this.commonReads = commonReads;
        this.struct = struct;
        this.src = src;
        iter = new AtomicInteger();
        this.log = log;
        this.start = start;
        this.end = end;
    }

    @Override
    public Boolean call() {
        if (log) {
            if (iter.incrementAndGet() % 100 == 0)
            System.out.print("\r" + iter.get());
        }
        if(i < start || i > end){
            return null;
        }
        for (int j = i; j < src.referenceLength; j++) {
            if (j < start){
                continue;
            }
            if (j> end){
                break;
            }
            if (i == j) {
                commonReads[i][j] = struct.readsAtPosition[i].length;
                continue;
            }
            commonReads[i][j] = AlgorithmUtils.getSortedArraysIntersectionCount(struct.readsAtPosition[i], struct.readsAtPosition[j]);
            commonReads[j][i] = commonReads[i][j];
        }
        return null;
    }
}
