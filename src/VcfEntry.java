/*
 * Methods for storing VCF v4.2 entries for structural variants
 * It allows easy access to main variant fields, plus some methods
 * to do more parsing and error-checking for things like type and length.
 */

import java.util.Arrays;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VcfEntry implements Comparable<VcfEntry>
{

	String originalLine;
	String[] tabTokens;
	String oldId;
	String key;
	
	HashSet<Integer> support;
	
	public VcfEntry(String line) throws Exception
	{
		support = new HashSet<Integer>();
		originalLine = line;
		tabTokens = line.split("\t");
		if(tabTokens.length < 8)
		{
			throw new Exception("VCF line had too few entries: "
					+ Arrays.toString(tabTokens));
		}
		setKey();
	}
	
	void setKey() throws Exception
	{
		key = getChromosome() + "_" + String.format("%08d", getPos()) + "_" + getRef() + "_" + getAlt();
	}
	
	/*
	 * Reconstruct the VCF line by concatenating and tab-separating the fields
	 */
	public String toString()
	{
		StringBuilder sb = new StringBuilder("");
		for(int i = 0; i<tabTokens.length; i++)
		{
			sb.append(tabTokens[i]);
			if(i < tabTokens.length - 1)
			{
				sb.append("\t");
			}
		}
		return sb.toString();
	}
	
	/*
	 * Get the chromosome field
	 */
	public String getChromosome()
	{
		return tabTokens[0];
	}
	
	/*
	 * Set the chromosome field
	 */
	public void setChromosome(String s)
	{
		tabTokens[0] = s;
	}
	
	/*
	 * Get the POS field
	 */
	public int getPos() throws Exception
	{
		return Integer.parseInt(tabTokens[1]);
	}
	
	/*
	 * Set the POS field
	 */
	public void setPos(int val)
	{
		tabTokens[1] = val+"";
	}
	
	/*
	 * Get the variant ID field
	 */
	public String getId()
	{
		return tabTokens[2];
	}
	
	/*
	 * Set the variant ID field
	 */
	public void setId(String s)
	{
		tabTokens[2] = s;
	}
	
	/*
	 * Get the REF sequence field
	 */
	public String getRef()
	{
		return tabTokens[3];
	}
	
	/*
	 * Set the REF sequence field
	 */
	public void setRef(String s)
	{
		tabTokens[3] = s;
	}
	
	/*
	 * Get the ALT sequence field
	 */
	public String getAlt()
	{
		return tabTokens[4];
	}
	
	/*
	 * Set the ALT sequence field
	 */
	public void setAlt(String s)
	{
		tabTokens[4] = s;
	}
	
	/*
	 * Set a particular VCF INFO field, adding the field if it doesn't already exist
	 */
	public void setInfo(String field, String val) throws Exception
	{
		if(tabTokens[7].equals("."))
		{
			tabTokens[7] = field + "=" + val;
			return;
		}
		String[] infoFields = tabTokens[7].split(";");
		for(String semitoken : infoFields)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				String updatedToken = key + "=" + val;
				
				// Special case if this is the first INFO field
				if(tabTokens[7].startsWith(semitoken))
				{
					tabTokens[7] = tabTokens[7].replaceFirst(Pattern.quote(semitoken), Matcher.quoteReplacement(updatedToken));
				}
				else
				{
					tabTokens[7] = tabTokens[7].replace(";" + semitoken, ";" + updatedToken);
				}
				return;
			}
		}
		
		// Field not found, so add it!
		tabTokens[7] += ";" + field + "=" + val;
	}
	
	/*
	 * Get the value of a particular INFO field
	 */
	public String getInfo(String field) throws Exception
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				return semitoken.substring(1 + equalIndex);
			}
		}
		return "";
	}
	
	/*
	 * Whether or not the line has a particular INFO field
	 */
	public boolean hasInfoField(String fieldName)
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(fieldName))
			{
				return true;
			}
		}
		return false;
	}
	
	// Merges v into this variant
	void merge(VcfEntry v) throws Exception
	{
		for(int sample: v.support)
		{
			support.add(sample);
		}
		
		String[] infoFields = v.tabTokens[7].split(";");
		for(String semitoken : infoFields)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			String val = semitoken.substring(1 + equalIndex);
			
			if(!hasInfoField(key))
			{
				setInfo(key, val);
			}
		}
	}

	public int compareTo(VcfEntry o)
	{
		return key.compareTo(o.key);
	}
	
}
