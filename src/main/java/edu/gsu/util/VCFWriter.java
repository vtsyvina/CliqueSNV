package edu.gsu.util;

import edu.gsu.model.Clique;
import edu.gsu.start.Start;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

/**
 * Class is responsible for output for Variant Calling
 */
public class VCFWriter {

    /**
     * Writes SNP info into vcf format. Skips deletions SNPs(not supported by format).
     *
     * @param input      Input sam file
     * @param output     Output vcf file
     * @param variants   Clique with SNPs that were got for input file after applying CliqueSNV
     * @param consensus  Consensus string for given input file
     */
    public static void create(File input, File output, Clique variants, String consensus) {
        SamReader in = SamReaderFactory.make().open(input);
        VariantContextWriter writer = new VariantContextWriterBuilder()
                .setOutputFile(output)
                .setOutputFileType(VariantContextWriterBuilder.OutputType.VCF)
                .setReferenceDictionary(in.getFileHeader().getSequenceDictionary())
                .build();
        VariantContextBuilder builder = new VariantContextBuilder();
        //merge all SNPs by position
        Map<Integer, Set<Container>> snps = new HashMap<>();
        for (int i = 0; i < variants.snps.size(); i++) {
            snps.putIfAbsent(variants.snps.get(i), new HashSet<>());
            snps.get(variants.snps.get(i)).add(new Container(variants.minors.substring(i, i + 1), variants.frequencies.get(i)));
        }
        VCFHeader vcfHeader = new VCFHeader();

        VCFHeaderLine vcfHeaderLine = new VCFInfoHeaderLine("AF", VCFHeaderLineCount.A, VCFHeaderLineType.Float, "Allele frequency");
        vcfHeader.addMetaDataLine(vcfHeaderLine);
        writer.writeHeader(vcfHeader);
        int start = Integer.parseInt(Start.settings.getOrDefault("-os","0"));
        int end = Integer.parseInt(Start.settings.getOrDefault("-oe","10000000"));
        //add SNPs as Alleles to writer
        snps.entrySet().stream().sorted(Comparator.comparingInt(Map.Entry::getKey)).forEach(e -> {
            Integer pos = e.getKey();
            if (pos < start || pos > end){
                return;
            }
            if (e.getValue().stream().anyMatch(s -> !s.allele.equals("-")) && consensus.charAt(pos) != '-') {
                List<Allele> alleles = new ArrayList<>();
                alleles.add(Allele.create(String.valueOf(consensus.charAt(pos)), true));
                StringBuilder freq = new StringBuilder();
                e.getValue().stream().filter(s -> !s.allele.equals("-")).forEach(a -> {
                    alleles.add(Allele.create(a.allele, false));
                    freq.append(Utils.smartDoubleToString(a.frequency)).append(",");
                });
                //remove last comma
                freq.setLength(freq.length()-1);
                VariantContext variantContext = builder.alleles(alleles)
                        .chr("ref")
                        .start(pos + 1)
                        .stop(pos + 1)
                        .attribute("AF", freq.toString())
                        .make();
                writer.add(variantContext);
            }
        });
        writer.close();
    }

    private static class Container {
        String allele;
        double frequency;

        Container(String allele, double frequency) {
            this.allele = allele;
            this.frequency = frequency;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Container container = (Container) o;
            return Double.compare(container.frequency, frequency) == 0 &&
                    Objects.equals(allele, container.allele);
        }

        @Override
        public int hashCode() {
            return Objects.hash(allele, frequency);
        }
    }
}
