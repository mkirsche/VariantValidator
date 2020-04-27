/*
 * Converts a post-filtering TSV to two separate VCFs - one with all variants and one with only consensus variants
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class TableToVcf
{
	static String tableFn = "", consensusFn = "", allFn = "";
	
	/*
	 * Prints usage message
	 */
	static void usage()
	{
		System.out.println("Usage: java -cp src TableToVcf [args]");
		System.out.println("  Example: java -cp src TableToVcf table_file=merged.vcf consensus_file=consensus.vcf all_file=all.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  table_file     (String) - vcf file containing the variants after merging across samples");
		System.out.println("  consensus_file (String) - file to output consensus variants to");
		System.out.println("  all_file       (String) - file to output all variants to");
		System.out.println();
	}
	
	/*
	 * Parse command line arguments
	 */
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
				if(key.equalsIgnoreCase("table_file")) { tableFn = val; }
				else if(key.equalsIgnoreCase("consensus_file")) { consensusFn = val; } 
				else if(key.equalsIgnoreCase("all_file")) { allFn = val; } 

			}
		}
		
		if(tableFn.length() == 0 || consensusFn.length() == 0 || allFn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		Scanner input = new Scanner(new FileInputStream(new File(tableFn)));
		VariantTable variantTable = new VariantTable(input.nextLine());
		while(input.hasNext())
		{
			variantTable.addRow(input.nextLine());
		}
		input.close();
		
		PrintWriter consensusOut = new PrintWriter(new File(consensusFn));
		PrintWriter allOut = new PrintWriter(new File(allFn));
		
		// Print VCF headers
		printVcfHeader(consensusOut);
		printVcfHeader(allOut);
		
		// Print VCF entries based on the table
		for(int i = 0; i<variantTable.rows.size(); i++)
		{
			variantTable.printRowVcf(i, consensusOut, true);
			variantTable.printRowVcf(i, allOut, false);
		}
		
		consensusOut.close();
		allOut.close();
	}
	
	/*
	 * Prints a standard VCF header to a given PrintWriter
	 */
	static void printVcfHeader(PrintWriter out)
	{
		out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
	}
	
	/*
	 * A table of variants parsed from a post-filtering TSV
	 */
	static class VariantTable
	{
		// Map of category names to which column they correspond to
		HashMap<String, Integer> categoryToIndex;
		
		// A list of data entries (rows in the table)
		ArrayList<String[]> rows;
		
		/*
		 * Parses the header line and initializes a table with those fields
		 */
		VariantTable(String headerLine)
		{
			categoryToIndex = new HashMap<String, Integer>();
			String[] categories = headerLine.split("\t");
			for(int i = 0; i<categories.length; i++)
			{
				categoryToIndex.put(categories[i].toLowerCase(), i);
			}
			
			rows = new ArrayList<String[]>();
		}
		
		/*
		 * Adds a row to the table based on a TSV row
		 */
		void addRow(String line)
		{
			// Ignore empty lines
			if(line.length() == 0)
			{
				return;
			}
			String[] tokens = line.split("\t");
			
			rows.add(tokens);
		}
		
		/*
		 * Gets a particular field's value in a given row
		 */
		String getValue(int rowIndex, String category)
		{
			return rows.get(rowIndex)[categoryToIndex.get(category.toLowerCase())];
		}
		
		/*
		 * Prints the row in VCF format to the given PrintWriter
		 * If the consensusOnly flag is true, don't print the variant if it's non-consensus
		 */
		void printRowVcf(int rowIndex, PrintWriter out, boolean consensusOnly)
		{
			boolean inConsensus = getValue(rowIndex, "in_consensus").equalsIgnoreCase("true");
			if(!inConsensus && consensusOnly)
			{
				return;
			}
			String chr = getValue(rowIndex, "chrom");
			int pos = Integer.parseInt(getValue(rowIndex, "pos"));
			String ref = getValue(rowIndex, "ref");
			String alt = getValue(rowIndex, "alt");
			String id = ".";
			String filter = ".";
			String qual = ".";
			String info = ".";
			out.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", chr, pos, id, ref, alt, qual, filter, info);
		}
	}
}
