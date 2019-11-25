import org.jdom.*;
import org.jdom.output.XMLOutputter;
import org.jdom.input.*;
import java.math.BigInteger;
import java.io.*;
import java.util.*;


public class FilterLowQualityChro
{
    public static void main(String args[]) throws Exception
    {
	if(args.length<2)
	{
	    System.out.println("Usage: FilterLowQuanityChro census_chro.xml low_quality_filtered_peptides.txt");
	    ///home/gcantin/gcantin_on_data/lee/wei_yi/SILAC_phos_expers_112006/L-EGF_H_112206/light_parc

	    System.exit(0);
	}

	BufferedReader br = new BufferedReader(new FileReader(args[1]));
	String eachLine = "";
	String proLine = "";

	Set set = new HashSet();
	while( null != (eachLine= br.readLine()) )
	{
	    if(eachLine.startsWith("P"))
	    {
		String[] pArr = eachLine.split("\t");
		proLine = pArr[1];
	    }

	    if(eachLine.startsWith("&S"))
	    {
		String[] pArr = eachLine.split("\t");
		set.add(proLine + "_" + pArr[2]);
	    }
	}
	SAXBuilder builder = new SAXBuilder();

	Document doc = builder.build(new File(args[0]));
	Element rootEle = doc.getRootElement();

	ArrayList<Element> proList = new ArrayList<Element>();
	int offset=0;

	//for(Iterator<Element> itr=rootEle.getChildren("protein").iterator(); itr.hasNext(); )
	List<Element> proteinList = rootEle.getChildren("protein");
	for(int i=0;i<proteinList.size();i++)
	{
	    Element proEle = proteinList.get(i);
	    //Element proEle = itr.next();

	    proLine = proEle.getAttributeValue("locus");

	    int removePepCount=0;
	    int pepSize=0;

	    if(proEle.getChildren("peptide").size()<=0)
		proList.add(proEle);

	    for(Iterator<Element> pepitr=proEle.getChildren("peptide").iterator(); pepitr.hasNext(); )
	    {
		pepSize++;
		Element pepEle = pepitr.next();

		//System.out.println(pepEle + " " + pepEle.getAttributeValue("seq"));
		if( !set.contains(proLine + "_" + pepEle.getAttributeValue("seq")) )
		{
		    pepitr.remove();
		    removePepCount++;
		}
	    }

/*
	    if(proLine.equals("IPI00193918.1"))
	    {
		System.out.println(proLine + " " + pepSize + " " + removePepCount);	
	    }
*/	    
	    if(removePepCount==pepSize && pepSize>0)
	    {
		int remSize = proList.size();
		int remIndex = i;
//		System.out.println("--00000000000" + " " + i + " " + remSize + " " + offset + " " + proList  + " " + remIndex);

		for(int j=0;j<remSize+1;j++)
		{
		    proteinList.remove(i); 
		    i--;
		}

		offset = offset + remSize + 1;

/*
		for(Iterator<Element> tmpProItr=proList.iterator(); tmpProItr.hasNext(); )
		{
		    tmpProItr.next();


		    proteinList.remove(i - proList.size()); // - count-1);
		    offset++;

//		    Element e = tmpProItr.next();

//		    tmpProItr.remove();
		}
*/

		proList = new ArrayList<Element>();	
	
	    }
	    else if(removePepCount!=pepSize && pepSize>0)
		proList = new ArrayList<Element>();	

	    
	}

	try {
	    OutputStream out = new FileOutputStream("census_chro_low_quality.xml");
	    org.jdom.output.XMLOutputter serializer = new org.jdom.output.XMLOutputter();
//	    serializer.setNewlines(true);
	    serializer.output(doc, out);

	    out.close();
	}
	catch (IOException e) {
	    System.err.println(e);
	}

    }
}
