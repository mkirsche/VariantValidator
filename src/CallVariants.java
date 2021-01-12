/*
 * Code for calling SNPs based on simple thresholds with samtools mpileup
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class CallVariants {
	static int maxLen = 31000;
	static String pileupFn = "", ofn = "";
	static int covThreshold = 20;
	static double refThreshold = .6;
	static double altThreshold = .15;
	static double indelThreshold = .4;
	static String flagPrefix = "";
	
	/*
	 * Prints out usage instructions
	 */
	static void usage()
	{
		System.out.println("Usage: java -cp src CallVariants [args]");
		System.out.println("  Example: java -cp src CallVariants pileup_file=jhu004.mpileup out_file=calls.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  pileup_file (String) - the output of samtools mpileup with the read alignments and reference");
		System.out.println("  out_file    (String) - the output file to call SNPs");

		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  coverage_threshold  (int)    [20]    - the min coverage needed to possibly flag a region");
		System.out.println("  genome_max_len      (int)    [31000] - an upper bound on the genome length");
		System.out.println("  alt_threshold       (float)  [0.15]  - call a variant if any alt allele frequency > this value");
		System.out.println("  ref_threshold       (float)  [0.60]  - call an N even if no alt allele frequency is high enough if ref allele frequency < this value");
		System.out.println("  indel_threshold       (float)  [0.40]  - call a variant an indel if no other variant is called there and indel frequency > this value");
		System.out.println("  flag_prefix         (String) []      - add this to AF and STRANDAF flag names");

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
				
				if(key.equals("pileup_file"))
				{
					pileupFn = val;
				}
				else if(key.equals("out_file"))
				{
					ofn = val;
				}
				else if(key.equals("alt_threshold"))
				{
					altThreshold = Double.parseDouble(val);
				}
				else if(key.equals("ref_threshold"))
				{
					refThreshold = Double.parseDouble(val);
				}
				else if(key.equals("indel_threshold"))
				{
					indelThreshold = Double.parseDouble(val);
				}
				else if(key.equals("genome_max_len"))
				{
					maxLen = Integer.parseInt(val);
				}
				else if(key.equals("coverage_threshold"))
				{
					covThreshold = Integer.parseInt(val);
				}
				else if(key.equals("flag_prefix"))
				{
					flagPrefix = val;
				}
			}
		}
		if(pileupFn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
public static void main(String[] args) throws Exception
{
	parseArgs(args);
	
	// Parse the cigar strings of read alignments and count up allele frequencies 
	System.err.println("Counting coverage from alignments");
	Scanner input = new Scanner(new FileInputStream(new File(pileupFn)));
	HashMap<String, int[][][]> cov = new HashMap<String, int[][][]>();
	HashMap<String, String[]> lines = new HashMap<String, String[]>();
	HashMap<String, char[]> genome = new HashMap<String, char[]>();
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.startsWith("@"))
		{
			continue;
		}
		//System.out.println(line);
		String[] tokens = line.split("\t");
		
		String chrName = tokens[0];
		int refPos = Integer.parseInt(tokens[1]) - 1;

		char refChar = tokens[2].charAt(0);
		
		if(!cov.containsKey(chrName))
		{
			cov.put(chrName, new int[maxLen][3][6]);
			lines.put(chrName, new String[maxLen]);
			genome.put(chrName, new char[maxLen]);
		}
		
		genome.get(chrName)[refPos] = refChar;
		
		int[][][] covArray = cov.get(chrName);
		
		covArray[refPos] = getAlleleFreqs(refChar, tokens[4]);
		
		lines.get(chrName)[refPos] = line;
	}
	input.close();
	
	// Now go through every position and output a variant if the allele frequencies indicate a variant
	System.err.println("Calling variants");
	
	PrintWriter out = new PrintWriter(ofn);
	
	// Loop over every ref contig
	int varId = 0;
	for(String s : cov.keySet())
	{
		// Get an array of allele frequencies for this contig
		int[][][] covArray = cov.get(s);
		HashSet<Integer> deleted = new HashSet<Integer>();
		for(int i = 0; i<covArray.length; i++)
		{
			// Total coverage over this position only counting matches/mismatches
			int totalCov = 0;
			for(int j = 0; j<covArray[i][0].length; j++) totalCov += covArray[i][0][j];
			
			if(totalCov < covThreshold) continue;
			
			// The character in the reference at this position
			int refChar = charToInt(genome.get(s)[i]);
			
			// Check for possible ALT alleles with high enough frequency
			int alt = -1;
			for(int j = 0; j<4; j++)
			{
				if(j != refChar && covArray[i][0][j] >= totalCov * altThreshold)
				{
					// Just in case multiple alts qualify, take the more frequent one
					if(alt != -1 && covArray[i][0][alt] > covArray[i][0][j])
					{
						continue;
					}
					alt = j;
				}
			}
			
			String indelSeq = "";
			
			if(alt == -1 && (covArray[i][0][5] >= totalCov * indelThreshold || covArray[i][0][6] >= totalCov * indelThreshold))
			{
				String x = getIndelSeq(lines.get(s)[i]);
				if(x.length() > 0)
				{
					if(covArray[i][0][5] >= totalCov * indelThreshold)
					{
						alt = 5;
					}
					else if(covArray[i][0][6] >= totalCov * indelThreshold)
					{
						for(int j = i; j<i+x.length(); j++)
						{
							deleted.add(j+1);
						}
						alt = 6;
					}
					System.out.println("Indel at position " + i + ": "+Arrays.toString(covArray[i][0])+" "+totalCov);
					indelSeq = x;
				}
			}
			
			if(alt == -1 && covArray[i][0][refChar] < (totalCov) * refThreshold && !deleted.contains(i))
			{
				System.out.println("Calling N at " + i + " " + Arrays.toString(covArray[i][0]) + " " + refChar);
				alt = 4;
			}
			
			if(alt != -1)
			{
				int totalPositive = 0, totalNegative = 0;
				for(int j = 0; j<covArray[i][1].length; j++)
				{
					totalPositive += covArray[i][1][j];
					totalNegative += covArray[i][2][j];
				}
				String refString = genome.get(s)[i] + "";
				String altString = intToChar(alt) + "";
				if(alt == 6)
				{
					altString = refString;
					refString = refString + indelSeq;
				}
				else if(alt == 5)
				{
					altString = refString + indelSeq;
				}
				out.printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;%s\n",
						s,
						i+1,
						"var" + varId,
						refString,
						altString,
						".",
						".",
						flagPrefix + "AF=" + String.format("%.6f", 1.0 * covArray[i][0][alt] / totalCov),
						flagPrefix + "STRANDAF=" + String.format("%d,%d,%d,%d", 
								covArray[i][1][alt], totalPositive, covArray[i][2][alt], totalNegative));
				varId++;
			}
		}
	}
	
	out.close();
}

static String getIndelSeq(String pileup)
{
	ArrayList<String> seqs = new ArrayList<String>();
	for(int i = 0; i<pileup.length(); i++)
	{
		char c = pileup.charAt(i);
		
		if(c == '+' || c == '-')
		{
			int end = i;
			int length = 0;
			while(end+1 < pileup.length() && pileup.charAt(end+1) >= '0' && pileup.charAt(end+1)<= '9')
			{
				end++;
				length = length * 10 + pileup.charAt(end) - '0';
			}
			String seq = pileup.substring(end+1, end+1+length);
			seqs.add(seq.toUpperCase());
		}
	}
	HashMap<String, Integer> seqFreq = new HashMap<String, Integer>();
	for(String x : seqs)
	{
		seqFreq.put(x, 1 + seqFreq.getOrDefault(x, 0));
	}
	for(String x : seqFreq.keySet())
	{
		if(seqFreq.get(x) > .8 * seqs.size())
		{
			return x;
		}
	}
	return "";
}

/*
 * Gets the number of A/C/G/T/N's covering a position from an mpileup string
 */
static int[][] getAlleleFreqs(char refChar, String pileup)
{
	int[][] res = new int[3][7];
	for(int i = 0; i<res.length; i++)
	{
		res[i] = new int[7];
	}
	for(int i = 0; i<pileup.length(); i++)
	{
		char c = pileup.charAt(i);
		
		// Exact match so use ref character
		if(c == '.' || c == ',')
		{
			res[0][charToInt(refChar)]++;
			if(c == '.')
			{
				res[1][charToInt(refChar)]++;
			}
			else
			{
				res[2][charToInt(refChar)]++;
			}
		}
		
		// Insertion or deletion after this base so ignore
		else if(c == '+' || c == '-')
		{
			int end = i;
			int length = 0;
			while(end+1 < pileup.length() && pileup.charAt(end+1) >= '0' && pileup.charAt(end+1)<= '9')
			{
				end++;
				length = length * 10 + pileup.charAt(end) - '0';
			}
			boolean capital = (pileup.charAt(end+1) >= 'A' && pileup.charAt(end+1) <= 'Z') || pileup.charAt(end+1) == '*';
			i = end + length;
			
			int idx = 5;
			if(c == '-') idx = 6;
			res[0][idx]++;
			if(capital)
			{
				res[1][idx]++;
			}
			else
			{
				res[2][idx]++;
			}
		}
		
		else if(c == '*')
		{
			res[0][5]++;
			res[1][5]++;
		}
		
		else if(c == '#')
		{
			res[0][5]++;
			res[2][5]++;
		}
		
		// Last character indicator - ignore
		else if(c == '$')
		{
			continue;
		}
		
		// First character indicator - ignore
		else if(c == '^')
		{
			i++;
			continue;
		}
		
		// Mismatch or N so count this character after converting it to an integer
		else
		{
			int val = charToInt(c);
			if(val != -1)
			{
				res[0][charToInt(c)]++;
				
				if(Character.isUpperCase(c) || c == '>')
				{
					res[1][charToInt(c)]++;
				}
				if(Character.isLowerCase(c) || c == '<')
				{
					res[2][charToInt(c)]++;
				}
			}
		}
	}
	return res;
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
	else if(c == '>' || c == '<' || c == 'n' || c == 'N') return 4;
	else return -1;
}

static char intToChar(int val)
{
	if(val == 0) return 'A';
	else if(val == 1) return 'C';
	else if(val == 2) return 'G';
	else if(val == 3) return 'T';
	else if(val == 5) return '.';
	else return 'N';
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
