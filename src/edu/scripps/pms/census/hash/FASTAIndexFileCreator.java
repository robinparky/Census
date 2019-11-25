/*
 * FASTAIndexFileCreator.java
 *
 * Created on March 25, 2005, 10:36 AM
 */

package edu.scripps.pms.census.hash;

import gnu.trove.TIntLongHashMap;

/**
 * @author Robin Park
 * @version $Id: FASTAIndexFileCreator.java,v 1.1 2013/01/15 22:44:20 rpark Exp $
 */
import java.io.*;
import edu.scripps.pms.census.util.RelExFileFilter;

public class FASTAIndexFileCreator {
    
    private final static String MS1_FILE = "ms1";
    private final static String MS2_FILE = "ms2";
    private final static String MSZM = "mszm";
  
    private String filePath;
    private static java.text.NumberFormat n = java.text.NumberFormat.getInstance();

    static 
    {
	n.setMaximumFractionDigits(3);
    }
       

    private static void printError()
    {
        System.out.println("Usage : java FASTAIndexFileCreator [-option]");
        System.out.println("-a\t Generate index files for all ms files in the same folder");
        System.out.println("-f\t ms file name");

    }


    public static void main(String args[]) throws IOException
    {
      
	if(args.length<1 || ("-f".equals(args[0]) && args.length<2) )
	{
            printError();

            return;
	}
        
	System.out.println("indexing files...");

        if("-a".equals(args[0]))
        {
            File f = new File(".");
            String[] list = f.list();

            for(int i=0;i<list.length;i++)
            {
                if(null != list[i] && (list[i].endsWith(MS1_FILE) || list[i].endsWith(MS2_FILE) || list[i].endsWith(MSZM)) )
                {
                    createIndexFile(list[i]);
                    System.out.println(list[i] + " was indexed.");
                }

            }

            return;
        } else if("-f".equals(args[0]))
        {
            createIndexFile(args[1]);
            System.out.println(args[1] + " was indexed.");
        }
        else
        {
           printError();
           return;
        }
            
	System.out.println(" completed");

    }


    public FASTAIndexFileCreator(String filePath)
    {
	this.filePath = filePath;
    }

    public static void createIndexFile(String fileName) throws IOException
    {
	try {
        if( fileName.endsWith("fasta"))
            createFASTAIndexFile(fileName);
        else if( fileName.endsWith("ms2") )
            createMS2IndexFile(fileName);

	} catch(Exception e)
	{
	    e.printStackTrace();
	    throw new IOException("Failed to create index file. " + e);
	}

    }

    private static void createFASTAIndexFile(String fileName) throws IOException
    {
        File file = new File(fileName);
        InputStream fisIn = new FileInputStream(fileName);

        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String lastLine;

        FileOutputStream out = new FileOutputStream(fileName + ".index");
        PrintStream p = new PrintStream( out );
        String temp;

        int size = (int)file.length();

        byte[] byteBuffer = new byte[size];

        fisIn.read(byteBuffer, 0, size);
        StringBuffer sb = new StringBuffer(); 

        long pos;
        
        for(int i=0;i<byteBuffer.length;i++)
        {
//		System.out.print("\tc = " + byteBuffer[i]);
            if( (char)byteBuffer[i] == '>' && ((i==0) || (char)byteBuffer[i-1] == '\n'))
            {
//		System.out.println("\tfound > byteIndex = "+i+"\n");
                pos = i;
                int j = i+1; //skip S char and space
                char ch = (char)byteBuffer[j];
                int tabCount=0;
                 
                while( ch != '\n')
                {
                    if(ch == '\t')
                        tabCount++;

         //           if(tabCount==2)
         //           {
                        sb.append(ch);
         //           }

                    ch = (char)byteBuffer[++j];

                }
//		System.out.println("message from fastaIndexing: sb = [" + sb.toString() +"]");
		String accession = extractAccession(sb.toString());
//		System.out.println("message from fastaIndexing: accession = " + accession);
                p.print(accession);
            //    System.out.println(sb.toString());
                p.print("\t");
                p.println(pos);
            //    System.out.println(pos);
                
                sb.delete(0, sb.length());
            }
        }

        fisIn.close();
        p.close();
        out.close();
    }

    public static String extractAccession(String accession)
    {
        //NCBI, IPI, or others such as UNIT_PROT, SGD, NCI
//        accession = getDefline().substring( getDefline().indexOf('>')+1 );
        //accession = getDefline();

        //There are many corruptted sqt file.  Ignore it.
        try
        {
            if( accession.startsWith("gi") && accession.contains("|") ) //NCBI
            {
                String[] arr = accession.split("\\|");

                if( arr.length>=4 && ("gb".equals(arr[2]) || "ref".equals(arr[2]) || "emb".equals(arr[2]) || "dbj".equals(arr[2]) || "prf".equals(arr[2]) ||"sp".equals(arr[2])) || "tpd".equals(arr[2]) ||"tpg".equals(arr[2]) ||"tpe".equals(arr[2]) )
                    accession = arr[3];
                else
                {
                    arr = accession.split(" ");
                    accession = arr[1];
                }

                //Accession # should end with digit.  If accession # does not end with digit,
                //grap next string (We assume this next one ends with digit.)
                /*
                if( pattern.matcher(arr[3]).matches() )
                    accession = arr[3];
                else
                    accession = arr[4].substring(0, arr[4].indexOf(" "));
                */

            }
            else if( accession.startsWith("IPI") ) //IPI
            {
                String arr[] = accession.split("\\|");
                String subArr[] = arr[0].split(":");

                if(subArr.length>1)
                    accession = subArr[1];
                else
                    accession = subArr[0];
            }
            else if( accession.startsWith("Re") || accession.startsWith("contam") || accession.startsWith("Contam")) //Reverse database
            {
                int space = accession.indexOf(" ");
                int tab = accession.indexOf("\t");

                if(space<0) space = 40;
                if(tab<0) tab = 40;

                int index = (tab>space)?space:tab;

                int end;

                if(index<=0 || index>=40) //no space
                {
                    int length = accession.length();
                    end = (length>40)?40:length;
                }
                else  //cut by the first space
                    end = index;

                accession = accession.substring(0, end);
            }
            else //UNIT_PROT, NCI or SGD

            {
                int spaceIndex = accession.indexOf(" ");
                int tabIndex;

                if(spaceIndex>0)
                {
                    tabIndex = accession.indexOf("\t");

                    if(tabIndex>0 && spaceIndex>tabIndex)
                        accession = accession.substring(0, tabIndex);
                    else
                        accession = accession.substring(0, spaceIndex);
                }
            }
        }
        catch(Exception e)
        {
            //System.out.println("No Correct Accession found, but this will be handled by MSP system." + accession + " " +  e);

            int i = accession.indexOf(" ");
            if(i<0)
                return accession;
            else
                return accession.substring(0, i);

        }

        return accession;
    }

    private static void createMS2IndexFile(String fileName) throws IOException
    {
        File file = new File(fileName);
        InputStream fisIn = new FileInputStream(fileName);

        BufferedReader br = new BufferedReader(new FileReader(fileName));
	
        String lastLine = br.readLine();
	String[] strArr = lastLine.split("\t");

	boolean isOldVersion = false;

	if(strArr.length>3)
	{
	    System.out.println("===" + lastLine);
	    lastLine = strArr[2].substring( strArr[2].lastIndexOf("/")+1 );
	    lastLine = lastLine.substring( 0, lastLine.indexOf(" ") );
	    int year = Integer.parseInt(lastLine);

	    isOldVersion = (year<2006)?true:false;
	}

        FileOutputStream out = new FileOutputStream(fileName + ".index");
        PrintStream p = new PrintStream( out );

        String temp;

        int size = (int)file.length();

        byte[] byteBuffer = new byte[size];

        fisIn.read(byteBuffer, 0, size);
        StringBuffer sb = new StringBuffer(); 
        StringBuffer precurSb = new StringBuffer(); 

        long pos;
	double prevRet=-1;
	double inc=0.000d;

        for(int i=0;i<byteBuffer.length;i++)
        {
            if( (char)byteBuffer[i] == 'S' && (char)byteBuffer[i-1] == '\n')
            {
                pos = i;
                int j = i; //skip S char and space
                char ch = (char)byteBuffer[j];
                int tabCount=0;
                 
                while( ch != '\n')
                {
                    if(ch == '\t')
                        tabCount++;

                    if(tabCount==2)
                    {
                        sb.append(ch);
                    }
                    else if(tabCount==3)
                    {
                       precurSb.append(ch); 
                    }

                    ch = (char)byteBuffer[++j];

                }

                p.print(Integer.parseInt(sb.toString().trim()));
                p.print("\t");
                p.print(pos);
                p.print("\t");
                p.print( Float.parseFloat(precurSb.toString().trim()) );

		if(isOldVersion)
		    p.println("");
		else
		    p.print("\t");

                
                sb.delete(0, sb.length());
                precurSb.delete(0, precurSb.length());

            }

            if( !isOldVersion && (char)byteBuffer[i] == 'I' && 
		(char)byteBuffer[i+1] == '\t' &&
		(char)byteBuffer[i+2] == 'R' &&
		(char)byteBuffer[i+3] == 'e' &&
		(char)byteBuffer[i+4] == 't' &&
		(char)byteBuffer[i+5] == 'T' )
            {

                pos = i;
                int j = i;
                char ch = (char)byteBuffer[j];
                int tabCount=0;
                 
                while( ch != '\n')
                {
                    if(ch == '\t')
                        tabCount++;

                    if(tabCount==2)
                    {
                        sb.append(ch);
                    }

                    ch = (char)byteBuffer[++j];

                }

		double curRet = Double.parseDouble(sb.toString().trim());
		if(curRet == prevRet)
		{
		    inc += 0.001d;
		    curRet += inc;
		}
		else
		{
		    inc = 0.000d;
		    prevRet = curRet;
		}

		p.println( n.format(curRet) );
                
                sb.delete(0, sb.length());
            }

        }

        p.close();
        out.close();
        fisIn.close();
    }
}
