package edu.gsu.util.builders;

import com.carrotsearch.hppc.IntArrayList;
import edu.gsu.algorithm.AbstractSNV;
import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.SNVStructure;
import edu.gsu.model.Sample;
import edu.gsu.util.Utils;

public class SNVStructureBuilder {
    private static int minorCount = AbstractSNV.minorCount;
    private static String al  = AbstractSNV.al;

    /**
     * Builds data structure for PacBio reads
     *
     * @param src        source sample
     * @param srcProfile source profile
     * @return SNV structure with lists of minor positions
     */
    public static SNVStructure buildPacBio(Sample src, double[][] srcProfile, String consensus) {
        SNVStructure result = new SNVStructure();

        IntArrayList[] rows = initIntArrayList(src.reads[0].length() * minorCount);
        result.rowMinors = new int[src.reads[0].length() * minorCount][];

        IntArrayList[] cols = initIntArrayList(src.reads.length);
        result.colMinors = new int[src.reads.length][];

        result.rowN = new int[src.reads[0].length()][];
        IntArrayList[] rowN = initIntArrayList(src.reads[0].length());

        result.majorsInRow = new int[src.reads[0].length() * minorCount];
        for (int i = 0; i < src.reads.length; i++) {
            fillRowN(src.reads[i], 0, i, rowN);
            fillMajorsCount(src.reads[i], 0, consensus, result.majorsInRow);
        }

        for (int i = 0; i < src.reads.length; i++) {
            fillRowsAndCols(src.reads[i], 0, i, consensus, cols, rows);
        }
        result.rowN = copyFromIntListToArray(rowN);
        result.colMinors = copyFromIntListToArray(cols);
        result.rowMinors = copyFromIntListToArray(rows);
        result.profile = srcProfile;
        return result;
    }

    public static SNVStructure buildIllumina(IlluminaSNVSample src, String consensus, double[][] srcProfile) {
        SNVStructure result = new SNVStructure();
        IntArrayList[] rows = initIntArrayList(src.referenceLength * minorCount);
        result.rowMinors = new int[src.referenceLength * minorCount][];

        IntArrayList[] cols = initIntArrayList(src.reads.size());
        result.colMinors = new int[src.reads.size()][];

        result.rowN = new int[src.referenceLength][];
        IntArrayList[] rowN = initIntArrayList(src.referenceLength);

        result.majorsInRow = new int[src.referenceLength];
        IntArrayList[] readsAtPositions = initIntArrayList(src.referenceLength);
        for (int i = 0; i < src.reads.size(); i++) {
            PairEndRead read = src.reads.get(i);
            fillMajorsCount(read.l, read.lOffset, consensus, result.majorsInRow);
            fillMajorsCount(read.r, read.rOffset, consensus, result.majorsInRow);
            fillRowN(read.l, read.lOffset, i, rowN);
            fillRowN(read.r, read.rOffset, i, rowN);
            fillReadsAtPosition(read.l, read.lOffset, i, readsAtPositions);
            fillReadsAtPosition(read.r, read.rOffset, i, readsAtPositions);
        }
        for (int i = 0; i < src.reads.size(); i++) {
            PairEndRead read = src.reads.get(i);
            fillRowsAndCols(read.l, read.lOffset, i, consensus, cols, rows);
            fillRowsAndCols(read.r, read.rOffset, i, consensus, cols, rows);
        }

        result.rowN = copyFromIntListToArray(rowN);
        result.colMinors = copyFromIntListToArray(cols);
        result.rowMinors = copyFromIntListToArray(rows);
        result.readsAtPosition = copyFromIntListToArray(readsAtPositions);
        result.profile = srcProfile;
        return result;
    }

    private static void fillRowsAndCols(String read, int offset, int readNumber, String consensus, IntArrayList[] cols, IntArrayList[] rows) {
        for (int j = 0; j < read.length(); j++) {
            if (read.charAt(j) != consensus.charAt(offset+j) && read.charAt(j) != 'N') {
                int position = Utils.splittedPosition(offset + j, read.charAt(j), minorCount, al, consensus);
                cols[readNumber].add(position);
                rows[position].add(readNumber);
            }
        }
    }

    private static void fillRowN(String read, int offset, int readNumber, IntArrayList[] rowN) {
        for (int j = 0; j < read.length(); j++) {
            if (read.charAt(j) == 'N') {
                rowN[j+offset].add(readNumber);
            }
        }
    }

    private static void fillReadsAtPosition(String read, int offset, int readNumber, IntArrayList[] readsAtPosition) {
        for (int j = 0; j < read.length(); j++) {
            if (read.charAt(j) != 'N') {
                readsAtPosition[j + offset].add(readNumber);
            }
        }
    }

    private static IntArrayList[] initIntArrayList(int l){
        IntArrayList[] result = new IntArrayList[l];
        for (int i = 0; i < l; i++) {
            result[i] = new IntArrayList();
        }
        return result;
    }

    private static int[][] copyFromIntListToArray(IntArrayList[] src){
        int[][] result = new int[src.length][];
        for (int i = 0; i < src.length; i++) {
            result[i] = src[i].toArray();
        }
        return result;
    }

    private static void fillMajorsCount(String read, int offset, String consensus, int[] majorCount){
        for (int i = 0; i < read.length(); i++) {
            if (read.charAt(i) == consensus.charAt(i+offset)){
                majorCount[i+offset]++;
            }
        }
    }
}
