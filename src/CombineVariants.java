import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeSet;

public class CombineVariants
{
	static String vcfFn = "", ofn = "";
	static String gffFn = "", genomeFn = "";
	static boolean usingGenes = false;
	
	// These are used when incorporating gene annotations
	static TreeSet<Integer> orfStarts;
	static HashMap<String, String> genome;
	
	static void usage()
	{
		System.out.println("Usage: java -cp src CombineVariants [args]");
		System.out.println("  Example: java -cp src CombineVariants vcf_file=merged.vcf out_file=merged_combined.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  vcf_file    (String) - vcf file containing the variants after merging across samples");
		System.out.println("  out_file    (String) - file to output variants after combining adjacent positions");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  gene_file   (String) - gff file containing genes: only groups variants together if they are in the same CDS reading frame");
		System.out.println("  genome_file (String) - path to genome, required if using a gene file");
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
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 
				else if(key.equalsIgnoreCase("gene_file")) { gffFn = val; usingGenes = true; } 
				else if(key.equals("genome_file")) { genomeFn = val; }
			}
		}
		
		if(vcfFn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(1);
		}
		
		if(usingGenes && genomeFn.length() == 0)
		{
			usage();
			System.exit(1);;
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		
		ArrayList<VcfEntry> allEntries = new ArrayList<VcfEntry>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				out.println(line);
				continue;
			}
			
			VcfEntry entry = new VcfEntry(line);
			
			allEntries.add(entry);
		}
		
		Collections.sort(allEntries);
		
		orfStarts = new TreeSet<Integer>();
		genome = new HashMap<String, String>();

		if(usingGenes)
		{
			Scanner geneInput = new Scanner(new FileInputStream(new File(gffFn)));
			while(geneInput.hasNext())
			{
				String line = geneInput.nextLine();
				if(line.startsWith("#"))
				{
					continue;
				}
				String[] tokens = line.split("\t");
				if(tokens.length < 5 || !tokens[2].equals("CDS"))
				{
					continue;
				}
				int start = Integer.parseInt(tokens[3]);
				int end = Integer.parseInt(tokens[4]);
				for(int i = start; i+2 <= end; i+=3)
				{
					orfStarts.add(i);
				}
			}
			geneInput.close();
			
			// Read in genome
			Scanner genomeInput = new Scanner(new FileInputStream(new File(genomeFn)));
			
			StringBuilder seq = new StringBuilder("");
			String refName = "";
			
			while(genomeInput.hasNext())
			{
				String line = genomeInput.nextLine();
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
			genomeInput.close();
		}
				
		// The current run of adjacent SNPs
		ArrayList<VcfEntry> currentEntries = new ArrayList<VcfEntry>();
		
		// Indels which occurred somewhere in the middle of the current run and need to be output afterwards
		ArrayList<VcfEntry> pendingEntries = new ArrayList<VcfEntry>();
		for(VcfEntry entry : allEntries)
		{
			boolean isSnp = entry.getRef().length() == 1 && entry.getAlt().length() == 1;
			
			boolean isAdjacent = currentEntries.size() > 0 && currentEntries.get(currentEntries.size()-1).getPos() + 1 >= entry.getPos();
			
			if(usingGenes)
			{
				isAdjacent = true;
				if(currentEntries.size() == 0)
				{
					isAdjacent = false;
				}
				else
				{
					Integer currentOrfStart = orfStarts.floor(currentEntries.get(0).getPos());
					Integer newOrfStart = orfStarts.floor(entry.getPos());
					
					if(currentOrfStart == null || newOrfStart == null)
					{
						isAdjacent = false;
					}
					else if(currentOrfStart + 2 < orfStarts.floor(currentEntries.get(0).getPos()))
					{
						isAdjacent = false;
					}
					else if(newOrfStart + 2 < entry.getPos())
					{
						isAdjacent = false;
					}
					else if(currentOrfStart != newOrfStart)
					{
						isAdjacent = false;
					}
				}
			}
			
			if(isSnp && (currentEntries.size() == 0 || isAdjacent))
			{
				currentEntries.add(entry);
			}
			
			else
			{
				if(currentEntries.size() > 0 && !isAdjacent)
				{
					processAdjacentVariants(currentEntries, out);
					currentEntries = new ArrayList<VcfEntry>();
					
					for(VcfEntry e : pendingEntries)
					{
						out.println(e);
					}
					pendingEntries = new ArrayList<VcfEntry>();
					
					if(isSnp)
					{
						currentEntries.add(entry);
					}
					else
					{
						out.println(entry);
					}
				}
				else if(currentEntries.size() > 0)
				{
					pendingEntries.add(entry);
				}
				else if(isSnp)
				{
					currentEntries.add(entry);
				}
				else
				{
					out.println(entry);
				}
			}
		}
		
		if(currentEntries.size() > 0)
		{
			processAdjacentVariants(currentEntries, out);
		}
		input.close();
		out.close();
	}
	
	/*
	 * Processes a list of adjacent variant calls
	 */
	static void processAdjacentVariants(ArrayList<VcfEntry> entries, PrintWriter out) throws Exception
	{
		int minPos = entries.get(0).getPos();
		int maxPos = entries.get(entries.size() - 1).getPos();
		
		// Get the reference character at each position
		char[] refs = new char[maxPos - minPos + 1];
		Arrays.fill(refs, '.');
		for(VcfEntry entry : entries)
		{
			refs[entry.getPos() - minPos] = entry.getRef().charAt(0);
		}
		
		if(genome != null && genome.size() > 0)
		{
			String chrName = entries.get(0).getChromosome();
			for(int i = 0; i<refs.length; i++)
			{
				if(refs[i] == '.')
				{
					refs[i] = genome.get(chrName).charAt(i + minPos - 1);
				}
				else if(refs[i] != genome.get(chrName).charAt(i + minPos - 1));
			}
		}
		
		int numSamples = entries.get(0).getInfo("SUPP_VEC").length();
		
		if(numSamples == 0)
		{
			char[] alt = new char[maxPos - minPos + 1];
			for(int i = 0; i<alt.length; i++)
			{
				alt[i] = refs[i];
			}
			for(VcfEntry entry : entries)
			{
				
				for(int i = 0; i<entry.getAlt().length(); i++)
				{
					alt[i + entry.getPos() - minPos] = entry.getAlt().charAt(i);
				}
			}
			
			VcfEntry copy = new VcfEntry(entries.get(0).originalLine);
			copy.setRef(new String(refs));
			copy.setAlt(new String(alt));
			out.println(copy);
		}
		
		// Create an alt sequence for each sample
		char[][] alts = new char[numSamples][maxPos - minPos + 1];
		
		// Initialize alt sequences to the ref sequence
		for(int i = 0; i<numSamples; i++)
			for(int j = 0; j<alts[i].length; j++)
			{
				alts[i][j] = refs[j];
			}
		
		// Go through the variant calls and update the ALT sequences
		for(VcfEntry entry : entries)
		{
			String suppVec = entry.getInfo("SUPP_VEC");
			
			for(int i = 0; i<numSamples; i++)
			{
				if(suppVec.charAt(i) == '0')
				{
					continue;
				}
				for(int j = 0; j<entry.getAlt().length(); j++)
				{
					alts[i][j + entry.getPos() - minPos] = entry.getAlt().charAt(j);
				}
			}
		}
		
		// Now create a set of ALT sequences being mapped to which samples they occur in
		HashMap<String, ArrayList<Integer>> altMap = new HashMap<String, ArrayList<Integer>>();
		for(int i = 0; i<numSamples; i++)
		{
			String altSequence = new String(alts[i]);
			if(!altMap.containsKey(altSequence))
			{
				altMap.put(altSequence, new ArrayList<Integer>());
			}
			altMap.get(altSequence).add(i);
			
		}
		
		// Create new entries
		for(String s : altMap.keySet())
		{
			VcfEntry copy = new VcfEntry(entries.get(0).originalLine);
			copy.setRef(new String(refs));
			copy.setAlt(s);
			
			// Ignore samples with no variants
			if(copy.getRef().equals(copy.getAlt()))
			{
				continue;
			}
			
			char[] suppVec = new char[numSamples];
			Arrays.fill(suppVec, '0');
			for(int sampleId : altMap.get(s))
			{
				suppVec[sampleId] = '1';
			}
			copy.setInfo("SUPP_VEC", new String(suppVec));
			out.println(copy);
		}
	}
}
