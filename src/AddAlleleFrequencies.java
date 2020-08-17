import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;

public class AddAlleleFrequencies {
	
	// The vcf file with merged variants
	static String vcfFn = "";
	
	// The mpileup file from the Illumina data
	static String illuminaMpileupFn = "";
	
	// The mpileup file form the ONT data
	static String ontMpileupFn = "";
	
	// File to write updated variants to
	static String ofn = "";
	
	// Max genome length allowed
	static int maxLen = 31000;
	
	static void usage()
	{
		System.out.println("Usage: java -cp src AddAlleleFrequencies [args]");
		System.out.println("  Example: java -cp src AddAlleleFrequencies vcf_file=merged.vcf illumina_mpileup=illumina.mpileup ont_mpileup=ont.mpileup out_file=updated.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  vcf_file         (String) - a VCF file with the merged variants");
		System.out.println("  ont_mpileup      (String) - mpileup from the Oxford Nanopore read alignments");
		System.out.println("  out_file         (String) - file to write updated variants to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  illumina_mpileup (String) - mpileup from the Illumina read alignments");
		System.out.println();
	}
	
	static void parseArgs(String[] args)
	{
		for(String s : args)
		{
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				
			}
			else
			{
				String key = s.substring(0, equalsIdx);
				String val = s.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("vcf_file")) { vcfFn = val; }
				else if(key.equalsIgnoreCase("illumina_mpileup")) { illuminaMpileupFn = val; } 
				else if(key.equalsIgnoreCase("ont_mpileup")) { ontMpileupFn = val; }
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 

			}
		}
		
		if(vcfFn.length() == 0 || ofn.length() == 0 || (ontMpileupFn.length() == 0 && illuminaMpileupFn.length() == 0))
		{
			usage();
			System.exit(1);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		Mpileup ontMpileup = null;

        if(ontMpileupFn.length() > 0)
		{
			ontMpileup = new Mpileup(ontMpileupFn);
		}
		
		Mpileup illuminaMpileup = null;
		
		if(illuminaMpileupFn.length() > 0)
		{
			illuminaMpileup = new Mpileup(illuminaMpileupFn);
		}
		
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			if(line.startsWith("#"))
			{
				out.println(line);
				continue;
			}
			VcfEntry entry = new VcfEntry(line);
			
			if(illuminaMpileup == null)
			{
				addInfoFieldsSingle(entry, ontMpileup, false);
			}
            else if(ontMpileup == null)
			{
				addInfoFieldsSingle(entry, illuminaMpileup, true);
			}
			else
			{
				addInfoFields(entry, illuminaMpileup, ontMpileup);
			}
			out.println(entry);
		}
		
		input.close();
		out.close();
		
	}
	
	static void addInfoFieldsSingle(VcfEntry entry, Mpileup ontMpileup, boolean illumina) throws Exception
	{
		String chrName = entry.getChromosome();
		int position = entry.getPos() - 1;
		
		// Get the mpileup data for this specific position
		int[][] ontCovArray = ontMpileup.allFrequencies.get(chrName)[position];
		
		String ref = entry.getRef();
		String alt = entry.getAlt();
		
		// Ignore indels - variants with len(ref) = len(alt) were already split by merging
		if(ref.length() == 1 && alt.length() == 1)
		{
			// Convert the ALT allele to an index for mpileup lookups
			int altVal = CallVariants.charToInt(alt.charAt(0));
			
			// Add up total depths for each technology and strand
			int[] ontTotals = new int[ontCovArray.length];
			
			// Compute ONT depths
			for(int i = 0; i<ontCovArray.length; i++)
			{
				for(int x : ontCovArray[i])
				{
					ontTotals[i] += x;
				}
				ontTotals[i] -= ontCovArray[i][5];
			}
			
			// ONT alt allele frequency - handle case with zero depth
			double ontAfValue = 0;
			if(ontTotals[0] != 0)
			{
				ontAfValue = 1.0 * ontCovArray[0][altVal] / ontTotals[0]; 
			}
			
			// Set ONT alt allele frequency INFO field
			String ontAf = String.format("%.6f", ontAfValue);
			entry.setInfo(illumina ? "ILLUMINA_AF" : "AF", ontAf);
			
			// Set ONT strand bias INFO field
			String ontStrandBias = String.format("%d,%d,%d,%d", 
					ontCovArray[1][altVal], ontTotals[1], 
					ontCovArray[2][altVal], ontTotals[2]);
			entry.setInfo(illumina ? "ILLUMINA_STRANDAF" : "STRANDAF", ontStrandBias);
			
			// Set fields for all alleles on each strand of ONT
			entry.setInfo(illumina ? "ILLUMINA_POSITIVE_STRAND_FREQUENCIES" : "POSITIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					ontCovArray[1][0], ontCovArray[1][1], 
					ontCovArray[1][2], ontCovArray[1][3], 
					ontCovArray[1][4], ontCovArray[1][5]));
			entry.setInfo(illumina ? "ILLUMINA_NEGATIVE_STRAND_FREQUENCIES" : "NEGATIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					ontCovArray[2][0], ontCovArray[2][1], 
					ontCovArray[2][2], ontCovArray[2][3], 
					ontCovArray[2][4], ontCovArray[2][5]));
			
		}
		
		else
		{
			// For indels set everything to 0
			entry.setInfo("ILLUMINA_AF", "0");
			entry.setInfo("ILLUMINA_STRANDAF", "0,0,0,0");
			entry.setInfo("ILLUMINA_POSITIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
			entry.setInfo("ILLUMINA_NEGATIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
		}
	}
	
	/*
	 * Adds allele-frequency INFO fields
	 */
	static void addInfoFields(VcfEntry entry, Mpileup illuminaMpileup, Mpileup ontMpileup) throws Exception
	{
		String chrName = entry.getChromosome();
		int position = entry.getPos() - 1;
		
		boolean inIllumina = illuminaMpileup.allFrequencies.containsKey(chrName);
		boolean inOnt = ontMpileup.allFrequencies.containsKey(chrName);
		if(!inIllumina && !inOnt)
		{
			return;
		}
		else if(!inIllumina)
		{
			int[][][] ontData = ontMpileup.allFrequencies.get(chrName);
			illuminaMpileup.allFrequencies.put(chrName, new int[ontData.length][ontData[0].length][ontData[0][0].length]);
		}
		else if(!inOnt)
		{
			int[][][] illuminaData = illuminaMpileup.allFrequencies.get(chrName);
			ontMpileup.allFrequencies.put(chrName, new int[illuminaData.length][illuminaData[0].length][illuminaData[0][0].length]);
		}
		
		// Get the mpileup data for this specific position
		int[][] illuminaCovArray = illuminaMpileup.allFrequencies.get(chrName)[position];
		int[][] ontCovArray = ontMpileup.allFrequencies.get(chrName)[position];
		
		String ref = entry.getRef();
		String alt = entry.getAlt();
		
		// Ignore indels - variants with len(ref) = len(alt) were already split by merging
		if(ref.length() == 1 && alt.length() == 1)
		{
			// Convert the ALT allele to an index for mpileup lookups
			int altVal = CallVariants.charToInt(alt.charAt(0));
			
			// Add up total depths for each technology and strand
			int[] illuminaTotals = new int[illuminaCovArray.length];
			int[] ontTotals = new int[ontCovArray.length];
			
			// Compute ONT depths
			for(int i = 0; i<ontCovArray.length; i++)
			{
				for(int x : ontCovArray[i])
				{
					ontTotals[i] += x;
				}
				ontTotals[i] -= ontCovArray[i][5];
			}
			
			// Compute Illumina depths
			for(int i = 0; i<illuminaCovArray.length; i++)
			{
				for(int x : illuminaCovArray[i])
				{
					illuminaTotals[i] += x;
				}
				illuminaTotals[i] -= illuminaCovArray[i][5];
			}
			
			// ONT alt allele frequency - handle case with zero depth
			double ontAfValue = 0;
			if(ontTotals[0] != 0)
			{
				ontAfValue = 1.0 * ontCovArray[0][altVal] / ontTotals[0]; 
			}
			
			// Set ONT alt allele frequency INFO field
			String ontAf = String.format("%.6f", ontAfValue);
			entry.setInfo("AF", ontAf);
			
			// Set ONT strand bias INFO field
			String ontStrandBias = String.format("%d,%d,%d,%d", 
					ontCovArray[1][altVal], ontTotals[1], 
					ontCovArray[2][altVal], ontTotals[2]);
			entry.setInfo("STRANDAF", ontStrandBias);
			
			// Illumina alt allele frequency - handle case with zero depth
			double illuminaAfValue = 0;
			if(illuminaTotals[0] != 0)
			{
				illuminaAfValue = 1.0 * illuminaCovArray[0][altVal] / illuminaTotals[0];
			}
			
			// Set Illumina alt allele frequency INFO field
			String illuminaAf = String.format("%.6f", illuminaAfValue);
			entry.setInfo("ILLUMINA_AF", illuminaAf);
			
			// Set Illumina strand bias INFO field
			String illuminaStrandBias = String.format("%d,%d,%d,%d", 
					illuminaCovArray[1][altVal], illuminaTotals[1], 
					illuminaCovArray[2][altVal], illuminaTotals[2]);
			entry.setInfo("ILLUMINA_STRANDAF", illuminaStrandBias);
			
			// Set fields for all alleles on each strand of ONT
			entry.setInfo("POSITIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					ontCovArray[1][0], ontCovArray[1][1], 
					ontCovArray[1][2], ontCovArray[1][3], 
					ontCovArray[1][4], ontCovArray[1][5]));
			entry.setInfo("NEGATIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					ontCovArray[2][0], ontCovArray[2][1], 
					ontCovArray[2][2], ontCovArray[2][3], 
					ontCovArray[2][4], ontCovArray[2][5]));
			
			// Set fields for all alleles on each strand of Illumina
			entry.setInfo("ILLUMINA_POSITIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					illuminaCovArray[1][0], illuminaCovArray[1][1], 
					illuminaCovArray[1][2], illuminaCovArray[1][3], 
					illuminaCovArray[1][4], illuminaCovArray[1][5]));
			entry.setInfo("ILLUMINA_NEGATIVE_STRAND_FREQUENCIES", String.format("%d,%d,%d,%d,%d,%d",
					illuminaCovArray[2][0], illuminaCovArray[2][1], 
					illuminaCovArray[2][2], illuminaCovArray[2][3], 
					illuminaCovArray[2][4], illuminaCovArray[2][5]));
		}
		
		else
		{
			// For indels set everything to 0
			entry.setInfo("AF", "0");
			entry.setInfo("STRANDAF", "0,0,0,0");
			entry.setInfo("ILLUMINA_AF", "0");
			entry.setInfo("ILLUMINA_STRANDAF", "0,0,0,0");
			entry.setInfo("POSITIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
			entry.setInfo("NEGATIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
			entry.setInfo("ILLUMINA_POSITIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
			entry.setInfo("ILLUMINA_NEGATIVE_STRAND_FREQUENCIES", "0,0,0,0,0,0");
		}

	}
	
	static class Mpileup
	{
		// Map chromosome name to an array of frequencies indexed by (position, strand, base)
		HashMap<String, int[][][]> allFrequencies;
		
		// The reference characters
		HashMap<String, char[]> genome;
		
		/*
		 * Take in an mpileup file and store the allele frequencies at each position
		 */
		Mpileup(String fn) throws Exception
		{
			Scanner input = new Scanner(new FileInputStream(new File(fn)));
			allFrequencies = new HashMap<String, int[][][]>();
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("@"))
				{
					continue;
				}
				
				String[] tokens = line.split("\t");
				
				// Get chromosome, position, and ref allele
				String chrName = tokens[0];
				int refPos = Integer.parseInt(tokens[1]) - 1;
				char refChar = tokens[2].charAt(0);
				
				if(!allFrequencies.containsKey(chrName))
				{
					allFrequencies.put(chrName, new int[maxLen][3][6]);
				}
				
				// Fill the frequency array at this position
				int[][][] covArray = allFrequencies.get(chrName);
				covArray[refPos] = CallVariants.getAlleleFreqs(refChar, tokens[4]);
			}
			input.close();
		}
	}

}
