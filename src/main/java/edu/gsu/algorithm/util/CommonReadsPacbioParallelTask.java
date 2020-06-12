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


    public CommonReadsPacbioParallelTask(int i, int[][] commonReads, SNVStructure struct, Sample src, boolean log) {
        this.i = i;
        this.commonReads = commonReads;
        this.struct = struct;
        this.src = src;
        iter = new AtomicInteger();
        this.log = log;
    }

    @Override
    public Boolean call() {
        if (log) {
            if (iter.incrementAndGet() % 100 == 0)
                System.out.print("\r" + iter.get());
        }
        for (int j = i; j < commonReads.length; j++) {
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
