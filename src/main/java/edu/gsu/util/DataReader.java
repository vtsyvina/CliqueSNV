package edu.gsu.util;

import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.PairEndRead;
import edu.gsu.model.Sample;
import edu.gsu.start.Start;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;

/**
 * Class to read input data from standard input files
 */
public class DataReader {

    private static final Pattern begin = Pattern.compile("^-*");
    private static final Pattern end = Pattern.compile("-*$");

    public static Sample readSample(File file) throws IOException {
        return readSample(file.toPath());
    }

    public static Sample readSample(Path filePath) throws IOException {
        return readList(filePath);
    }

    public static Sample readList(Path filePath) throws IOException {
        List<String> raw = Files.readAllLines(filePath);
        List<String> result = new ArrayList<>();
        StringBuilder seq = new StringBuilder();
        for (int j = 0; j < raw.size(); j++) {
            String str = raw.get(j);
            if (str.startsWith(">") && seq.length() > 0 || str.length() == 0) {
                result.add(seq.toString());
                seq.setLength(0);
            } else if (!str.startsWith(">")) {
                seq.append(str);
                //workaround for for cases when file doesn't contain last empty string
                if (j + 1 == raw.size()) {
                    result.add(seq.toString());
                }
            }
        }
        return new Sample(filePath.toFile().getName().replaceFirst("[.][^.]+$", ""), result.toArray(new String[0]));
    }

    /**
     * For imputation data
     */
    public static Sample readLineByLine(Path filePath) throws IOException{
        return new Sample(filePath.toFile().getName().replaceFirst("[.][^.]+$", ""), Files.readAllLines(filePath).toArray(new String[0]));
    }


    /**
     * Read reads from .fas or .sam file. For fas file all reads should be of equal size: 'N' should be added in the beginning and in the end if needed. Example:
     * >read1
     * ACTTTAAA-ACG
     * >read2
     * NNNTTACAANNN
     */
    public static Sample getPacBioReads(File file) throws IOException {
        String[] split = file.getName().split("\\.");
        if (split.length < 2) {
            return null;
        }
        String extension = split[split.length - 1];
        if (extension.equals("fas")) {
            return readSample(file);
        }
        SamReader open = SamReaderFactory.make().open(file);
        SAM4WebLogo sam4WebLogo = new SAM4WebLogo(open);
        List<String> tmp = new ArrayList<>();
        List<String> readNames = new ArrayList<>();
        for (SAMRecord anOpen : open) {
            StringBuilder str = new StringBuilder(sam4WebLogo.printRead(anOpen, false));
            int i = 0;
            while (i < str.length() && str.charAt(i) == '-') {
                str.setCharAt(i++, 'N');
            }
            i = str.length() - 1;
            while (i >= 0 && str.charAt(i) == '-') {
                str.setCharAt(i--, 'N');
            }
            //TODO think about it. getAlignmentStart is 1 for some reason everytime
            str.delete(0, 1);
            tmp.add(str.toString());
            readNames.add(anOpen.getReadName());
        }
        return new Sample(split[split.length - 2], tmp.toArray(new String[0]), readNames.toArray(new String[0]));
    }

    public static IlluminaSNVSample getIlluminaPairedReads(File file) {
        Map<String, List<SAMRecord>> readsSet = new HashMap<>();
        System.out.println("Start read sam");
        SamReader open = SamReaderFactory.make().open(file);
        int limit = Start.tryParseInt(Start.settings.get("-limit"), 1_000_000);
        for (SAMRecord anOpen : open) {

            if (readsSet.size() % 100_000 == 0) {
                System.out.print("\r" + readsSet.size());
            }
            if (anOpen.getReadUnmappedFlag()) {
                continue;
            }
            if (anOpen.getAlignmentStart() > limit){
                continue;
            }
            if (!readsSet.containsKey(anOpen.getReadName())) {
                readsSet.put(anOpen.getReadName(), new ArrayList<>());
            }
            readsSet.get(anOpen.getReadName()).add(anOpen);
        }
        System.out.println(" DONE");
        SAM4WebLogo sam4WebLogo = new SAM4WebLogo(open);
        System.out.println("Start convert");
        List<PairEndRead> pairedReads = new ArrayList<>(readsSet.size());
        int cores = Start.threadsNumber();
        List<Map<String, List<SAMRecord>>> parts = new ArrayList<>();
        for (int i = 0; i < cores; i++) {
            parts.add(new HashMap<>());
        }
        int[] index = {0};
        readsSet.forEach((key1, value1) -> parts.get(index[0]++ % cores).put(key1, value1));
        List<Callable<List<PairEndRead>>> tasks = new ArrayList<>();
        AtomicInteger counter = new AtomicInteger();
        for (int i = 0; i < cores; i++) {
            tasks.add(new ConvertTask(parts.get(i), sam4WebLogo, counter));
        }
        try {
            List<Future<List<PairEndRead>>> futures = Start.executor.invokeAll(tasks);
            futures.forEach(future -> {
                try {
                    pairedReads.addAll(future.get());
                } catch (InterruptedException | ExecutionException e) {
                    System.err.println("Error! Parallel tasks were not successful on get");
                    e.printStackTrace();
                }
            });
        } catch (InterruptedException e) {
            System.err.println("Error! Parallel tasks were not successful on invoke");
            e.printStackTrace();
        }
        System.out.println(" DONE");
        System.out.println("Total reads number: " + pairedReads.size());
        return new IlluminaSNVSample(file.getName().replaceFirst("[.][^.]+$", ""), pairedReads,
                pairedReads.parallelStream().mapToInt(r -> Math.max(r.lOffset + r.l.length(), r.rOffset + r.r.length())).max().orElse(0));
    }

    private static String getSplittedRead(String read, int offset, String consensus, String alphabet) {
        StringBuilder str = new StringBuilder();
        for (int j = 0; j < read.length(); j++) {
            int major = Utils.getMajorAllele(consensus, alphabet, offset + j);
            int minor = 0;
            for (int k = 0; k < alphabet.length() - 1; k++, minor++) {
                if (minor == major) {
                    minor++;
                }
                int allele = alphabet.indexOf(read.charAt(j));
                str.append(allele == minor ? "2" : "1");
            }
        }
        return str.toString();
    }

    private static class ConvertTask implements Callable<List<PairEndRead>> {
        Map<String, List<SAMRecord>> readsSet;
        SAM4WebLogo sam4WebLogo;
        AtomicInteger counter;

        ConvertTask(Map<String, List<SAMRecord>> readsSet, SAM4WebLogo sam4WebLogo, AtomicInteger counter) {
            this.readsSet = readsSet;
            this.sam4WebLogo = sam4WebLogo;
            this.counter = counter;
        }

        @Override
        public List<PairEndRead> call() {
            List<PairEndRead> pairedReads = new ArrayList<>(readsSet.size());
            readsSet.forEach((key, value) -> {
                if (value.size() == 1) {
                    String s = sam4WebLogo.printRead(value.get(0), true);
                    //s = cutRead(s, value.get(0).getAlignmentStart());//end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                    pairedReads.add(new PairEndRead(s,
                            "",
                            value.get(0).getAlignmentStart() - 1,
                            -1, key)
                    );
                } else {
                    for (int i = 0; i < value.size(); i++) {
                        //if read is unpaired
                        if (!value.get(i).getReadPairedFlag()) {
                            String s = sam4WebLogo.printRead(value.get(0), true);
                            //s = cutRead(s, value.get(0).getAlignmentStart());//end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                            pairedReads.add(new PairEndRead(s,
                                    "",
                                    value.get(0).getAlignmentStart() - 1,
                                    -1, key)
                            );
                            continue;
                        }
                        for (int j = 0; j < value.size(); j++) {
                            if (value.get(i).getMateAlignmentStart() == value.get(j).getAlignmentStart()
                                    && value.get(i).getAlignmentStart() <= value.get(j).getAlignmentStart()) {
                                if (i == j) {
                                    String s = sam4WebLogo.printRead(value.get(i), true);
                                    //s = cutRead(s, value.get(i).getAlignmentStart());//end.matcher(begin.matcher(s).replaceAll("")).replaceAll("");
                                    pairedReads.add(new PairEndRead(s,
                                            "",
                                            value.get(i).getAlignmentStart() - 1,
                                            -1, key)
                                    );
                                } else {
                                    String s = sam4WebLogo.printRead(value.get(i), true);
                                    //s =  cutRead(s, value.get(i).getAlignmentStart());
                                    String s2 = sam4WebLogo.printRead(value.get(j), true);
                                    //s2 = cutRead(s2, value.get(j).getAlignmentStart());//end.matcher(begin.matcher(s2).replaceAll("")).replaceAll("");
                                    pairedReads.add(new PairEndRead(s,
                                            s2,
                                            value.get(i).getAlignmentStart() - 1,
                                            value.get(j).getAlignmentStart() - 1,
                                            key)
                                    );
                                }
                            }
                        }
                    }
                }
                if (pairedReads.size() % 100_000 == 0) {
                    System.out.print("\r" + counter.addAndGet(100000));
                }
            });
            return pairedReads;
        }

        /**
         * Read comes with leading and ending '-', we need to cut them
         */
//        private String cutRead(String read, int alignmentStart){
//            read = read.substring(alignmentStart);
//            int end = read.length()-1;
//            for (int i = read.length()-1; i > 0; i--) {
//                if (read.charAt(i) != '-'){
//                    end = i+1;
//                    break;
//                }
//            }
//            return read.substring(0,end);
//        }

    }
}
