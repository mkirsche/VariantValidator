/*
 * Converts an iVar TSV to a VCF
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

public class IvarToVcf
{
	static String tableFn = "", ofn = "";
	
	/*
	 * Prints usage message
	 */
	static void usage()
	{
		System.out.println("Usage: java -cp src IvarToVcf [args]");
		System.out.println("  Example: java -cp src IvarToVcf table_file=ivar.tsv out_file=all.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  table_file     (String) - tsv output from ivar variant caller");
		System.out.println("  out_file       (String) - vcf file to output all variants to");
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
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 

			}
		}
		
		if(tableFn.length() == 0 || ofn.length() == 0)
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
		
		PrintWriter out = new PrintWriter(new File(ofn));
		
		// Print VCF headers
		printVcfHeader(out);
		
		// Print VCF entries based on the table
		for(int i = 0; i<variantTable.rows.size(); i++)
		{
			variantTable.printRowVcf(i, out, false);
		}
		
		out.close();
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
			String chr = getValue(rowIndex, "region");
			int pos = Integer.parseInt(getValue(rowIndex, "pos"));
			String ref = getValue(rowIndex, "ref");
			String alt = getValue(rowIndex, "alt");
			String id = ".";
			String filter = ".";
			String qual = ".";
			String info = String.format("IVAR_REF_DP=%s;IVAR_REF_RV=%s;IVAR_REF_QUAL=%s;IVAR_ALT_DP=%s;IVAR_ALT_RV=%s;"
					+ "IVAR_ALT_QUAL=%s;IVAR_ALT_FREQ=%s;IVAR_TOTAL_DP=%s;IVAR_PVAL=%s",
					getValue(rowIndex, "ref_dp"),
					getValue(rowIndex, "ref_rv"),
					getValue(rowIndex, "ref_qual"),
					getValue(rowIndex, "alt_dp"),
					getValue(rowIndex, "alt_rv"),
					getValue(rowIndex, "alt_qual"),
					getValue(rowIndex, "alt_freq"),
					getValue(rowIndex, "total_dp"),
					getValue(rowIndex, "pval")
			);
			out.printf("%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", chr, pos, id, ref, alt, qual, filter, info);
		}
	}
}
