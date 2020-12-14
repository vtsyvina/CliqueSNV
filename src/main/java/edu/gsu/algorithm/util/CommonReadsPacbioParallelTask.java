package edu.gsu.algorithm.util;

import edu.gsu.model.SNVStructure;
import edu.gsu.model.Sample;

import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;

public class CommonReadsPacbioParallelTask implements Callable<Boolean> {
    private static AtomicInteger iter = new AtomicInteger();
    private int i;
    private int[][] commonReads;
    private SNVStructure struct;
    private Sample src;
    private boolean log;
    // start and end positoins where we search for SNPs
    private int start;
    private int end;


    public CommonReadsPacbioParallelTask(int i, int[][] commonReads, SNVStructure struct, Sample src, int start, int end, boolean log) {
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
        for (int j = i; j < commonReads.length; j++) {
            if (j < start){
                continue;
            }
            if (j> end){
                break;
            }
            if (i == j) {
                commonReads[i][j] = src.reads.length - struct.rowN[i].length;
                continue;
            }
            commonReads[i][j] = src.reads.length - struct.rowN[i].length - struct.rowN[j].length + AlgorithmUtils.getSortedArraysIntersectionCount(struct.rowN[i], struct.rowN[j]);
            commonReads[j][i] = commonReads[i][j];
        }
        return true;
    }
}
