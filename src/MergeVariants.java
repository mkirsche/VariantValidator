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
	
	// Bam file with illumina reads
	static String illuminaBam = "";
	
	static void usage()
	{
		System.out.println("Usage: java -cp src MergeVariants [args]");
		System.out.println("  Example: java -cp src MergeVariants file_list=vcflist.txt out_file=merged.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list    (String) - a txt file containing absolute paths to VCF files, one on each line");
		System.out.println("  out_file     (String) - file to write merged variants to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  illumina_bam (String) - a filename which will be added to the header as ILLUMINABAM");
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
				else if(key.equalsIgnoreCase("illumina_bam")) { illuminaBam = val; } 

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
				VcfEntry entry = new VcfEntry(line);
				entry.tabTokens[7] = entry.tabTokens[7].replaceAll(";;", ";");
				
				if(entry.getAlt().startsWith("-"))
				{
					String oldRef = entry.getRef();
					entry.setRef(entry.getRef() + entry.getAlt().substring(1));
					entry.setAlt(oldRef);
				}
				else if(entry.getAlt().startsWith("+"))
				{
					entry.setAlt(entry.getRef() + entry.getAlt().substring(1));
				}
								
				while(entry.getRef().length() > 1 && entry.getAlt().length() > 1)
				{
					int refLength= entry.getRef().length();
					int altLength = entry.getAlt().length();
					if(entry.getRef().charAt(refLength - 1) == entry.getAlt().charAt(altLength - 1))
					{
						entry.setRef(entry.getRef().substring(0, refLength - 1));
						entry.setAlt(entry.getAlt().substring(0, altLength - 1));
					}
					else if(entry.getRef().substring(0, 1).equals(entry.getAlt().substring(0, 1)))
					{
						entry.setPos(1 + entry.getPos());
						entry.setRef(entry.getRef().substring(1));
						entry.setAlt(entry.getAlt().substring(1));
					}
					else
					{
						break;
					}
				}
				
				entry.setKey();
				
				VcfEntry[] splitEntries = split(entry);
				
				for(VcfEntry v : splitEntries)
				{
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
			}
			input.close();
		}
		
		int numSamples = vcfs.length;
		
		PrintWriter out = new PrintWriter(new File(ofn));
		
		String files = "";
		for(int i = 0; i<vcfs.length; i++)
		{
			files += vcfs[i];
			if(i<vcfs.length - 1) files += ",";
		}
		out.println("##filelist=" + files);
		if(illuminaBam.length() > 0)
		{
			out.println("##ILLUMINABAM=" + illuminaBam);
		}
		out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
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
			
			String[] tokens = entry.tabTokens;
			for(int i = 0; i<8; i++)
			{
				out.print(tokens[i] + (i == 7 ? "\n" : "\t"));
			}
		}
		
		out.close();
	}
	
	/*
	 * Splits an entry if it has reflen and altlen equal and greater than 1
	 */
	static VcfEntry[] split(VcfEntry entry) throws Exception
	{
		if(entry.getRef().length() == 1 || entry.getAlt().length() == 1)
		{
			return new VcfEntry[] {entry};
		}
		
		ArrayList<VcfEntry> res = new ArrayList<VcfEntry>();
		for(int i = 0; i<entry.getRef().length() && i < entry.getAlt().length(); i++)
		{
			VcfEntry cur = new VcfEntry(entry.toString());
			cur.setRef(entry.getRef().charAt(i) + "");
			cur.setAlt(entry.getAlt().charAt(i) + "");
			if(i == entry.getRef().length() - 1 
					&& entry.getAlt().length() > entry.getRef().length())
			{
				cur.setAlt(entry.getAlt().substring(i));
			}
			if(i == entry.getAlt().length() - 1 
					&& entry.getRef().length() > entry.getAlt().length())
			{
				cur.setRef(entry.getRef().substring(i));
			}
			cur.setPos(entry.getPos() + i);
			cur.setId(entry.getId());
			if(!entry.getId().equals("."))
			{
				cur.setId(entry.getId() + "_" + i);
			}
			cur.setKey();
			if(cur.getRef().equalsIgnoreCase(cur.getAlt()))
			{
				continue;
			}
			res.add(cur);
		}
		VcfEntry[] array = new VcfEntry[res.size()];
		for(int i = 0; i<res.size(); i++)
		{
			array[i] = res.get(i);
		}
		return array;
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
