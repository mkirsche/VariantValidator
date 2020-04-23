import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeSet;

public class CheckVariants {
	static int maxLen = 31000;
	static String samFn = "", genomeFn = "", vcfFn = "";
	static int covThreshold = 20;
	static double fnThreshold = .6;
	static double fpThreshold = .4;
	
	/*
	 * Prints out usage instructions
	 */
	static void usage()
	{
		System.out.println("Usage: java -cp src CheckVariants [args]");
		System.out.println("  Example: java -cp src CheckVariants sam_file=jhu004.sam vcf_file=jhu004.vcf genome_file=ref.fa");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  sam_file    (String) - a SAM file with the alignments of the reads");
		System.out.println("  vcf_file    (String) - a VCF file with the variant calls");
		System.out.println("  genome_file (String) - a FASTA file with the reference genome");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  coverage_threshold  (int)    [20]    - the min coverage needed to possibly flag a region");
		System.out.println("  genome_max_len      (int)    [31000] - an upper bound on the genome length");
		System.out.println("  missed_variant_freq (float)  [0.6]   - call a possible missed variant if ref allele frequency < this value (or wrong alt allele if its frequency < this value)");
		System.out.println("  fp_freq             (float)  [0.4]   - call a false positive if ref allele frequency > this value");

		System.out.println();
	}
	
	/*
	 * Parses command line arguments
	 */
	static void parseArgs(String[] args)
	{
		for(String str : args)
		{
			String s = str;
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
			}
			else
			{
				String key = s.substring(0, equalsIdx).toLowerCase();
				String val = s.substring(1 + equalsIdx);
				
				if(key.equals("sam_file"))
				{
					samFn = val;
				}
				else if(key.equals("vcf_file"))
				{
					vcfFn = val;
				}
				else if(key.equals("genome_file"))
				{
					genomeFn = val;
				}
				else if(key.equals("coverage_threshold"))
				{
					covThreshold = Integer.parseInt(val);
				}
				else if(key.equals("genome_max_len"))
				{
					maxLen = Integer.parseInt(val);
				}
				else if(key.equals("missed_variant_freq"))
				{
					fnThreshold = Double.parseDouble(val);
				}
				else if(key.equals("fp_freq"))
				{
					fpThreshold = Double.parseDouble(val);
				}
			}
		}
		if(samFn.length() == 0 || vcfFn.length() == 0 || genomeFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
public static void main(String[] args) throws Exception
{
	parseArgs(args);
	
	// Read in variants 
	System.err.println("Reading variants");
	TreeSet<Variant> vars = new TreeSet<Variant>();
	Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0 || line.startsWith("#"))
		{
			continue;
		}
		
		Variant v = new Variant(line);
		if(vars.contains(v))
		{
			System.out.println("Multiple variants at position: " + v.chr + ":" + v.pos);
		}
		else
		{
			vars.add(v);
		}
	}
	input.close();
	
	// Read in genome
	System.err.println("Reading genome");

	HashMap<String, String> genome = new HashMap<String, String>();
	input = new Scanner(new FileInputStream(new File(genomeFn)));
	
	StringBuilder seq = new StringBuilder("");
	String refName = "";
	
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith(">"))
		{
			String name = line.split(" ")[0].substring(1);
			if(refName.length() > 0)
			{
				// add last contig
				genome.put(refName, seq.toString());
				seq = new StringBuilder("");
			}
			refName = name;
		}
		else
		{
			seq.append(line);
		}
	}
	if(refName.length() > 0)
	{
		// add last contig
		genome.put(refName, seq.toString());
	}
	input.close();
	
	// Parse the cigar strings of read alignments and count up allele frequencies 
	System.err.println("Counting coverage from alignments");
	input = new Scanner(new FileInputStream(new File(samFn)));
	HashMap<String, int[][]> cov = new HashMap<String, int[][]>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith("@"))
		{
			continue;
		}
		String[] tokens = line.split("\t");
		
		String chrName = tokens[2];
		if(!cov.containsKey(chrName))
		{
			cov.put(chrName, new int[maxLen][7]);
		}
		int[][] covArray = cov.get(chrName);
		
		String cigar = tokens[5];
		int n = cigar.length();
		
		String readSeq = tokens[9];
		
		int refPos = Integer.parseInt(tokens[3]);
		int queryPos = 0;
		
		int lenSoFar = 0;
		for(int i = 0; i<n; i++)
		{
			char c = cigar.charAt(i);
			if(c >= '0' && c <= '9')
			{
				lenSoFar = lenSoFar * 10 + (c - '0');
			}
			else
			{
				// deletion
				if(consumesReference(c) && !consumesQuery(c))
				{
					for(int j = 0; j<lenSoFar; j++)
					{
						covArray[refPos + j - 1][6]++;
					}
				}
				
				// insertion
				else if(consumesQuery(c) && !consumesReference(c))
				{
					covArray[refPos - 1][5]++;
				}
				else if(consumesQuery(c) && consumesReference(c))
				{
					for(int j = 0; j<lenSoFar; j++)
					{
						int charVal = charToInt(readSeq.charAt(queryPos + j));
						covArray[refPos + j - 1][charVal]++;
					}
				}
				
				if(consumesReference(c)) refPos += lenSoFar;
				if(consumesQuery(c)) queryPos += lenSoFar;
				lenSoFar = 0;
			}
		}
	}
	input.close();
	
	// Now go through every position and check consistency between allele frequencies and variant presence/absence
	
	System.err.println("Validating variants");
	
	// Loop over every ref contig
	for(String s : cov.keySet())
	{
		// Get an array of allele frequencies for this contig
		int[][] covArray = cov.get(s);
		for(int i = 0; i<covArray.length; i++)
		{
			// Total coverage over this position only counting matches/mismatches
			int totalCov = 0;
			for(int j = 0; j<5; j++) totalCov += covArray[i][j];
			
			if(totalCov < covThreshold) continue;
			
			int refChar = charToInt(genome.get(s).charAt(i));
			int refCov = covArray[i][refChar];
			
			boolean hasVar = vars.contains(new Variant(s, i));
			
			double refProp = 1.0 * refCov / totalCov;
			
			if(!hasVar && refProp < fnThreshold)
			{
				System.out.println("Possible missed variant at " + s + ":" + (i+1) + "; Ref allele = " + genome.get(s).charAt(i) + "; Ref proportion = " + String.format("%.3f", refProp) + "; Allele freqs = " + covToString(covArray[i]));
			}
			
			if(hasVar)
			{
				char alt = vars.ceiling(new Variant(s, i)).alt;
				int altChar = charToInt(alt);
				if(refProp > fpThreshold) 
				{
					System.out.println("Possible false positive at " + s + ":" + (i+1) + "; Ref allele = " + genome.get(s).charAt(i) + "; Alt allele = " + alt + "; Ref proportion = " + String.format("%.3f", refProp) + "; Allele freqs = " + covToString(covArray[i]));
				}
				else if(altChar < 4 && 1.0 * covArray[i][altChar] / totalCov < fnThreshold)
				{
					System.out.println("Possible wrong ALT at " + s + ":" + (i+1) + "; Alt allele = " + alt + "; Alt proportion = " + String.format("%.3f", 1.0 * covArray[i][altChar] / totalCov) + "; Allele freqs = " + covToString(covArray[i]));
				}
			}
		}
	}
}

/*
 * Prints the coverage array of a position in a human-readable format
 */
static String covToString(int[] cov)
{
	return String.format("A:%d, C:%d, G:%d, T:%d, N:%d, INS:%d, DEL:%d", cov[0], cov[1], cov[2], cov[3], cov[4], cov[5], cov[6]);
}

/*
 * Converts a basepair charater to an integer index
 */
static int charToInt(char c)
{
	if(c == 'a' || c == 'A') return 0;
	else if(c == 'c' || c == 'C') return 1;
	else if(c == 'g' || c == 'G') return 2;
	else if(c == 't' || c == 'T') return 3;
	return 4;
}

/*
 * Whether or not a CIGAR character consumes the reference
 */
static boolean consumesReference(char c)
{
	return c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X';
}

/*
 * Whether or not a CIGAR character consumes the reference
 */
static boolean consumesQuery(char c)
{
	return c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X';
}

/*
 * Stores the chr/pos/ref/alt of a variant
 */
static class Variant implements Comparable<Variant>
{
	String chr;
	int pos;
	char ref, alt;
	
	Variant(String chr, int pos)
	{
		this.chr = chr;
		this.pos = pos;
	}
	
	// Fill the variant from a vcf line
	Variant(String line)
	{
		String[] tokens = line.split("\t");
		chr = tokens[0];
		pos = Integer.parseInt(tokens[1])-1;
		ref = tokens[3].charAt(0);
		alt = tokens[4].charAt(0);
	}

	// Sort order is increasing by chromosome then position
	public int compareTo(Variant o)
	{
		if(!chr.equals(o.chr))
		{
			return chr.compareTo(o.chr);
		}
		return pos - o.pos;
	}
}
}