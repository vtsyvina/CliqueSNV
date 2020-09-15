package edu.gsu.start;

import edu.gsu.algorithm.SNVIlluminaMethod;
import edu.gsu.algorithm.SNVPacBioMethod;
import edu.gsu.algorithm.imputation.Imputation;
import edu.gsu.model.Clique;
import edu.gsu.model.IlluminaSNVSample;
import edu.gsu.model.SNVResultContainer;
import edu.gsu.model.Sample;
import edu.gsu.util.DataReader;
import edu.gsu.util.VCFWriter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;

public class Start {
    private static boolean log;
    public static Map<String, String> settings = new HashMap<>();
    public static Sample answer;
    public static ExecutorService executor;
    public static Sample coronaHaplotypes;

    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        printVersion();
        parseArgs(args);
        System.out.println("Settings: " + settings);
        executor = Executors.newFixedThreadPool(threadsNumber());
        if (settings.containsKey("-answer")) {
            answer = DataReader.readSample(new File(settings.get("-answer")));
        }
        if (settings.get("-help") != null) {
            helpOutput(null, false);
        }
        log = settings.containsKey("-log");
        try {
            if (settings.get("-m") != null) {
                switch (settings.get("-m")) {
                    case "snv-pacbio":
                        pacBioSNV(false);
                        break;
                    case "snv-pacbio-vc":
                        pacBioSNVVC();
                        break;
                    case "snv-illumina":
                        illumina2SNV(false);
                        break;
                    case "snv-illumina-vc":
                        illumina2SNVVC();
                        break;
                    case "imputation":
                        imputation(false);
                        break;
                    case "test":
                        multiTest();
                        break;
                    case "consensus-illumina":
                        consensus("illumina");
                        break;
                    case "consensus-pacbio":
                        consensus("pacbio");
                        break;
                    default:
                        helpOutput(settings.get("-m"), true);
                }
            } else {
                helpOutput("", false);
            }
        } finally {
            executor.shutdown();
        }
    }

    private static void helpOutput(String arg, boolean error) {
        if (error) {
            System.out.println("Error! Wrong argument value: " + arg);
        }
        System.out.println("How to set arguments:");
        System.out.println("-m snv-illumina -- run one of predefined methods. Methods are: snv-pacbio, snv-illumina, snv-pacbio-vc, snv-illumina-vc");
        System.out.println("-in /usr/name/tmp/reads.sam -- input file");
        System.out.println("-outDir /usr/name/tmp/ -- folder with output.");
        System.out.println("-threads 4 -- how many threads to use in parallel. Usually just the number of cores is the best choice");
        System.out.println("Final command can look as follows:");
        System.out.println("java -jar clique-snv.jar -m snv-illumina -in data/flu.sam -outDir output/");
        System.out.println("More on parameters read at  https://github.com/vtsyvina/CliqueSNV");
        System.exit(1);
    }

    private static void illumina2SNV(boolean test) throws IOException {
        long start = System.currentTimeMillis();
        String coronaPath = settings.getOrDefault("-corona", "");
        if (coronaPath.length() > 0) {
            coronaHaplotypes = DataReader.getPacBioReads(new File(coronaPath));
            System.out.println();
        }
        String pathname = settings.getOrDefault("-in", "data/Illumina_reads/reads.sam");
        IlluminaSNVSample sample = DataReader.getIlluminaPairedReads(new File(pathname));
        System.out.println("read " + (System.currentTimeMillis() - start));
        SNVIlluminaMethod snvIlluminaMethod = new SNVIlluminaMethod(sample,
                tryParseInt(settings.get("-t"), 10),
                tryParseDouble(settings.get("-tf"), 0.05),
                log);
        List<SNVResultContainer> haplotypes = snvIlluminaMethod.getHaplotypes();
        writeSNVResultsToFile(sample.name, haplotypes, null, false);
        System.out.printf("SNV got %d haplotypes\n%n", haplotypes.size());
        if (!test) System.out.println(haplotypes);
        if (test) {
            snvIlluminaMethod.setLog(true);
            snvIlluminaMethod.outputAnswerChecking(haplotypes);
        }
        System.out.println("time,ms " + (System.currentTimeMillis() - start));
    }

    private static void illumina2SNVVC() throws IOException {
        long start = System.currentTimeMillis();
        String pathname = settings.getOrDefault("-in", "data/Illumina_reads/reads.sam");
        File input = new File(pathname);
        if (!input.exists()) {
            System.out.printf("Input file %s does not exists%n", input.getCanonicalPath());
            return;
        }
        IlluminaSNVSample sample = DataReader.getIlluminaPairedReads(input);
        System.out.println("read " + (System.currentTimeMillis() - start));
        SNVIlluminaMethod snvIlluminaMethod = new SNVIlluminaMethod(sample,
                tryParseInt(settings.get("-t"), 10),
                tryParseDouble(settings.get("-tf"), 0.05),
                log);
        Clique vc = snvIlluminaMethod.getVC();
        writeVCFResults(input, sample.name, snvIlluminaMethod.consensus(), vc);
        System.out.println("time,ms " + (System.currentTimeMillis() - start));
    }

    private static void pacBioSNV(boolean test) throws IOException, ExecutionException, InterruptedException {
        File file = new File(settings.getOrDefault("-in", "data/PacBio_reads/reads.sam"));
        if (!file.exists()) {
            System.out.printf("Input file %s does not exists%n", file.getCanonicalPath());
            return;
        }
        Sample sample = DataReader.getPacBioReads(file);
        System.out.println("Reads number " + sample.reads.length);
        long start;
        start = System.currentTimeMillis();
        SNVPacBioMethod snvPacBioMethod = new SNVPacBioMethod(sample,
                tryParseInt(settings.get("-t"), 10),
                tryParseDouble(settings.get("-tf"), 0.05),
                log);
        List<SNVResultContainer> haplotypes = snvPacBioMethod.getHaplotypes();
        System.out.printf("SNV got %d haplotypes\n%n", haplotypes.size());
        if (!test) System.out.println(haplotypes);
        if (test) {
            snvPacBioMethod.setLog(true);
            snvPacBioMethod.outputAnswerChecking(haplotypes);
        }
        writeSNVResultsToFile(sample.name, haplotypes, sample, false);
        System.out.println("time,ms " + (System.currentTimeMillis() - start));
    }

    private static void consensus(String method) throws IOException {
        File file = new File(settings.getOrDefault("-in", "data/PacBio_reads/reads.sam"));
        if (!file.exists()) {
            System.out.printf("Input file %s does not exists%n", file.getCanonicalPath());
            return;
        }
        String consensus = "";
        String name = "";
        if (method.equals("pacbio")) {
            Sample sample = DataReader.getPacBioReads(file);
            SNVPacBioMethod snvPacBioMethod = new SNVPacBioMethod(sample, log);
            consensus = snvPacBioMethod.consensus();
            name = sample.name;
        } else {
            IlluminaSNVSample sample = DataReader.getIlluminaPairedReads(file);
            SNVIlluminaMethod snvIlluminaMethod = new SNVIlluminaMethod(sample, log);
            consensus = snvIlluminaMethod.consensus();
            name = sample.name;
        }
        List<SNVResultContainer> c = new ArrayList<>();
        c.add(new SNVResultContainer("", null, null, consensus));
        writeSNVResultsToFile(name + "_consensus", c, null, true);
    }

    private static void imputation(boolean test) throws IOException, InterruptedException {
        long start;
        start = System.currentTimeMillis();
        File file = new File(settings.getOrDefault("-in", "data/50plus50/50snps.1.txt"));
        if (!file.exists()) {
            System.out.printf("Input file %s does not exists%n", file.getCanonicalPath());
            return;
        }
        if (file.isDirectory()) {
            if (settings.get("-outDir") == null) {
                System.out.println(" Error! -outDir is not specified for multiple files input. This parameter is mandatory for multiple file input");
                return;
            }
            ImputationTask.length = file.listFiles().length;
            List<ImputationTask> tasks = new ArrayList<>();
            for (File f : file.listFiles()) {
                tasks.add(new ImputationTask(f));
            }
            ExecutorService executor = Executors.newFixedThreadPool(threadsNumber());
            List<Future<Integer>> futures = executor.invokeAll(tasks);
            executor.shutdown();
            futures.forEach(f -> {
                try {
                    f.get();
                } catch (InterruptedException | ExecutionException e) {
                    System.err.println("Error! Parallel tasks were not successful on get");
                    e.printStackTrace();
                }
            });
            System.out.println();
            System.out.println("Results are available at " + settings.get("-outDir"));
        } else {
            Sample sample = DataReader.readLineByLine(file.toPath());
            System.out.println("Reads number " + sample.reads.length);
            Imputation imputation = new Imputation(sample,
                    tryParseInt(settings.get("-t"), 10),
                    tryParseDouble(settings.get("-tf"), 0.05),
                    log);
            Sample imputedHaplotypes = imputation.getImputedHaplotypes();
            writeImputationResults(imputedHaplotypes, false);
        }
        System.out.println("time,ms " + (System.currentTimeMillis() - start));
    }

    private static void pacBioSNVVC() throws IOException, ExecutionException, InterruptedException {
        File input = new File(settings.getOrDefault("-in", "data/PacBio_reads/reads.sam"));
        if (!input.exists()) {
            System.out.println(String.format("Input file %s does not exists", input.getCanonicalPath()));
            return;
        }
        Sample sample = DataReader.getPacBioReads(input);
        long start;
        start = System.currentTimeMillis();
        SNVPacBioMethod snvPacBioMethod = new SNVPacBioMethod(sample,
                tryParseInt(settings.get("-t"), 10),
                tryParseDouble(settings.get("-tf"), 0.05),
                log);
        Clique vc = snvPacBioMethod.getVC();
        writeVCFResults(input, sample.name, snvPacBioMethod.consensus(), vc);
        System.out.println("time,ms " + (System.currentTimeMillis() - start));
    }

    /**
     * Method for debug only. Run file with all given parameters
     */
    private static void multiTest() throws IOException, ExecutionException, InterruptedException {
        File input = new File(settings.getOrDefault("-in", "data/tests.in"));
        List<String> strings = Files.readAllLines(input.toPath());
        for (String args : strings) {
            System.out.println("Start test with arguments: " + args);
            parseArgs(args.split(" "));
            switch (settings.get("-m")) {
                case "snv-pacbio":
                    pacBioSNV(true);
                    break;
                case "snv-pacbio-vc":
                    pacBioSNVVC();
                    break;
                case "snv-illumina":
                    illumina2SNV(true);
                    break;
                case "snv-illumina-vc":
                    illumina2SNVVC();
                    break;
                case "imputation":
                    imputation(true);
                    break;
                default:
                    helpOutput(settings.get("-m"), true);
            }
            System.out.println();
        }
    }

    private static void writeSNVResultsToFile(String name, List<SNVResultContainer> haplotypes, Sample sample, boolean fastaOnly) throws IOException {
        String snvOutput = settings.getOrDefault("-outDir", "snv_output/");
        boolean writeReadNames = settings.containsKey("-rn");
        List<SNVResultContainer> haplotypesCopy = new ArrayList<>();
        if (haplotypes.size() == 0) {
            System.out.println("CliqueSNV didn't find any haplotypes (too low coverage)");
            return;
        }
        haplotypes.forEach(h -> haplotypesCopy.add(h.copy()));
        if (!snvOutput.endsWith("/")) {
            snvOutput += "/";
        }
        Path path;
        int start = Integer.parseInt(Start.settings.getOrDefault("-os", "0"));
        int end = Integer.parseInt(Start.settings.getOrDefault("-oe", String.valueOf(haplotypes.get(0).haplotype.length())));
        haplotypesCopy.forEach(h -> h.haplotype = h.haplotype.substring(start, end));
        if (!fastaOnly) {
            path = preparePath(snvOutput + (name == null ? "snv_output.txt" : name + ".txt"));
            System.out.println("Results are available in: " + path.toFile().getCanonicalPath());
            Files.write(path, ("Used parameters:" + Start.settings.toString() + "\n").getBytes(), StandardOpenOption.WRITE);
            Files.write(path, String.format("SNV got %d haplotypes\n", haplotypesCopy.size()).getBytes(), StandardOpenOption.APPEND);
            Files.write(path, haplotypesCopy.toString().getBytes(), StandardOpenOption.APPEND);
        }
        path = preparePath(snvOutput + (name == null ? "snv_output.fasta" : name + ".fasta"));
        System.out.println("Results are available in: " + path.toFile().getCanonicalPath());
        Path finalPath = path;
        int[] i = {1};
        String format = settings.getOrDefault("-fdf", "short"); //extended
        boolean defaultFormat = format.equals("short");
        int precision = 2;
        if (!defaultFormat && !format.equals("extended")) {
            precision = Integer.parseInt(format.substring(8));
        }
        int finalPrecision = precision;
        Path readPath = writeReadNames ? preparePath(snvOutput + name + "_read_names.txt") : null;
        haplotypesCopy.forEach(h -> {
            try {
                String haplotypeName = ">" + (i[0]++) + "_fr_" + h.frequency;
                if (!defaultFormat) {
                    haplotypeName = ">" + name + "_" + (i[0]++) + "_" + String.format("%." + finalPrecision + "f", h.frequency);
                }
                Files.write(finalPath, (haplotypeName + "\n").getBytes(), StandardOpenOption.APPEND);
                Files.write(finalPath, h.haplotype.getBytes(), StandardOpenOption.APPEND);
                Files.write(finalPath, "\n\n".getBytes(), StandardOpenOption.APPEND);
                if (writeReadNames) {
                    writeReadNames(readPath, h, sample, haplotypeName);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

        });
    }

    private static void writeReadNames(Path readPath, SNVResultContainer h, Sample sample, String haplotypeName) throws IOException {
        Files.write(readPath, (haplotypeName + "\n").getBytes(), StandardOpenOption.APPEND);
        StringBuilder str = new StringBuilder();
        if (h.illuminaCluster != null) {
            for (int j = 0; j < h.illuminaCluster.size(); j++) {
                str.append(h.illuminaCluster.get(j).name).append("\n");
                if (j % 5000 == 0) {
                    Files.write(readPath, str.toString().getBytes(), StandardOpenOption.APPEND);
                    str.setLength(0);
                }
            }
            if (str.length() > 0) {
                Files.write(readPath, str.toString().getBytes(), StandardOpenOption.APPEND);
            }
        }
        if (h.pacBioCluster != null) {
            int j = 0;
            for (Map.Entry<Integer, String> entry : h.pacBioCluster.entrySet()) {
                str.append(sample.readNames[entry.getKey()]).append("\n");
                if (j % 5000 == 0) {
                    Files.write(readPath, str.toString().getBytes(), StandardOpenOption.APPEND);
                    str.setLength(0);
                }
                j++;
            }
            if (str.length() > 0) {
                Files.write(readPath, str.toString().getBytes(), StandardOpenOption.APPEND);
            }
        }
    }

    private static Path preparePath(String filePath) throws IOException {
        Path path = Paths.get(filePath);
        Files.deleteIfExists(path);
        if (path.toFile().getParentFile() != null) {
            path.toFile().getParentFile().mkdirs();
        }
        Files.createFile(path);
        return path;
    }

    private static void parseArgs(String[] args) throws IOException {
        settings = new HashMap<>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("-") && i + 1 < args.length && !args[i + 1].startsWith("-")) {
                settings.put(args[i], args[i + 1]);
                i++;
            } else {
                settings.put(args[i], "true");
            }
        }
        if (settings.containsKey("-answer")) {
            answer = DataReader.readSample(new File(settings.get("-answer")));
        }
        if (settings.get("-help") != null) {
            helpOutput(null, false);
        }
        log = settings.containsKey("-log");
    }

    private static void writeVCFResults(File input, String sampleName, String consensus, Clique vcClique) throws IOException {
        String snvOutput = settings.getOrDefault("-outDir", "snv_output/");
        if (!snvOutput.endsWith("/")) {
            snvOutput += "/";
        }
        Path path = preparePath(snvOutput + (sampleName == null ? "snv_output.vcf" : sampleName + ".vcf"));
        VCFWriter.create(input, path.toFile(), vcClique, consensus);
        System.out.println("Results are available in: " + path.toFile().getCanonicalPath());
    }

    private static void writeImputationResults(Sample sample, boolean multiple) throws IOException {
        String snvOutput = settings.getOrDefault("-outDir", "snv_output/");
        if (!snvOutput.endsWith("/")) {
            snvOutput += "/";
        }
        Path path = preparePath(snvOutput + (sample.name == null ? "imputation_output.txt" : sample.name + ".txt"));
        try (FileWriter writer = new FileWriter(path.toFile())) {
            for (String read : sample.reads) {
                writer.append(read).append("\n");
            }
        }
        if (!multiple) System.out.println("Results are available in: " + path.toFile().getCanonicalPath());
    }

    public static String getConsensusFromParameter() {
        if (settings.get("-ref") == null) {
            return null;
        }
        try {
            Sample sample = DataReader.readSample(new File(settings.get("-ref")));
            return sample.reads[0];
        } catch (IOException e) {
            try {
                System.out.println("File " + new File(settings.get("-ref")).getCanonicalPath() + " for reference doesn't not exist");
            } catch (IOException e1) {
                e1.printStackTrace();
            }
        }
        return null;
    }

    public static int threadsNumber() {
        int processors = Runtime.getRuntime().availableProcessors();
        String def = String.valueOf(processors);
        return tryParseInt(settings.getOrDefault("-threads", def), processors);
    }

    public static int tryParseInt(String value, int def) {
        if (value == null) {
            return def;
        }
        try {
            return Integer.parseInt(value);
        } catch (NumberFormatException nfe) {
            System.out.println("Couldn't parse " + value);
            return def;
        }
    }

    private static double tryParseDouble(String value, double def) {
        if (value == null) {
            return def;
        }
        try {
            return Double.parseDouble(value);
        } catch (NumberFormatException nfe) {
            System.out.println("Couldn't parse " + value);
            return def;
        }
    }

    private static void printVersion() throws IOException {
        final Properties properties = new Properties();
        properties.load(Start.class.getClassLoader().getResourceAsStream("project.properties"));
        System.out.println("CliqueSNV version: " + properties.getProperty("version"));
    }

    static class ImputationTask implements Callable<Integer> {
        static private AtomicInteger counter;
        static public int length;
        private File file;

        ImputationTask(File file) {
            counter = new AtomicInteger();
            this.file = file;
        }

        @Override
        public Integer call() throws Exception {
            Sample sample = DataReader.readLineByLine(file.toPath());
            Imputation imputation = new Imputation(sample,
                    tryParseInt(settings.get("-t"), 10),
                    tryParseDouble(settings.get("-tf"), 0.05),
                    log);
            Sample imputedHaplotypes = imputation.getImputedHaplotypes();
            writeImputationResults(imputedHaplotypes, true);
            int i = counter.incrementAndGet();
            System.out.print("\rProcessed " + i + " out of " + length);
            return i;
        }
    }
}
