import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.TreeSet;

public class MergeVariants
{
	// A file containing the absolute path to each VCF to merge
	static String fileList = "";
	
	// File to print merged variants to
	static String ofn = "";
	
	static void usage()
	{
		System.out.println("Usage: java -cp src MergeVariants [args]");
		System.out.println("  Example: java -cp src MergeVariants filelist=vcflist.txt out_file=merged.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list    (String) - a txt file containing absolute paths to VCF files, one on each line");
		System.out.println("  out_file    (String) - file to write merged variants to");

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
				if(key.equalsIgnoreCase("file_list")) { fileList = val; }
				else if(key.equalsIgnoreCase("out_file")) { ofn = val; } 
			}
		}
		
		if(fileList.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(1);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		String[] vcfs = getFilesFromList();
		TreeSet<VcfEntry> vars = new TreeSet<VcfEntry>();
		for(int i = 0; i<vcfs.length; i++)
		{
			String vcf = vcfs[i];
			Scanner input = new Scanner(new FileInputStream(new File(vcf)));
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.startsWith("#"))
				{
					continue;
				}
				VcfEntry v = new VcfEntry(line);
				v.support.add(i);
				
				if(vars.contains(v))
				{
					vars.floor(v).merge(v);
				}
				else
				{
					vars.add(v);
				}
			}
			input.close();
		}
		
		int numSamples = vcfs.length;
		
		PrintWriter out = new PrintWriter(new File(ofn));
		
		for(VcfEntry entry : vars)
		{
			char[] suppVec = new char[numSamples];
			Arrays.fill(suppVec, '0');
			for(int x : entry.support)
			{
				suppVec[x] = '1';
			}
			entry.setInfo("SUPP_VEC", new String(suppVec));
			entry.setInfo("SUPP", entry.support.size() + "");
			
			out.println(entry);
		}
		
		out.close();
	}
	
	static String[] getFilesFromList() throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(fileList)));
		ArrayList<String> res = new ArrayList<String>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() > 0)
			{
				res.add(line);
			}
		}
		String[] array = new String[res.size()];
		for(int i = 0; i<res.size(); i++)
		{
			array[i] = res.get(i);
		}
		input.close();
		return array;
	}
}
