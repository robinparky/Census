package scripts;

import java.sql.*;
import java.util.*;

public class CensusDB
{
    public static void main(String[] args) throws Exception
    {

String sql =
"select * from download_users where software_id=10";

    Connection con;
	//Class.forName("com.mysql.jdbc.Driver");
	Class.forName("org.gjt.mm.mysql.Driver");
	//innerConnection = DriverManager.getConnection("jdbc:mysql://localhost/pmsdb", "pmsuser", "myd474");
	con = DriverManager.getConnection("jdbc:mysql://homer.scripps.edu/newpmsdb", "pmsuser", "myd474");
	//con = DriverManager.getConnection("jdbc:mysql://skinner.scripps.edu/pmsdb", "pmsuser", "myd474");

	Statement stmt = con.createStatement();
	ResultSet rs = stmt.executeQuery(sql);

	Hashtable<String, Integer> ht = new Hashtable<String, Integer>();
	while(rs.next())
	{
	    
//	    System.out.println( rs.getString("email") );

	    String downloadDate = rs.getString("input_date").split(" ")[0];

	    Integer count = ht.get(downloadDate);
	    if(null == count)
	    {
		ht.put(downloadDate, new Integer(1));
	    } else {
		int countT = count.intValue();
		countT++;
		ht.put(downloadDate, countT);
	    }

	    
	}

	rs.close();

	stmt.close();
	con.close();

	Object[] arr = ht.keySet().toArray(); 
	Arrays.sort( arr );

	
	for(Object each : arr)
	    System.out.println( each + "\t" +ht.get(each) );

	System.out.println( "" );
	int total=0;
	int currentSeg = 1;
	int prevSeg = 1;


	String prev = "";

	for(Object eachObj : arr)
	{
	    String each = eachObj.toString();

	    String date = each.substring(each.length()-2);

	    int dateNum = Integer.parseInt(date);

	    if(dateNum<=10)
		currentSeg=1;
//		System.out.println("1\t" + dateNum);
	    else if(dateNum<=20)
		currentSeg=2;
	    else
		currentSeg=3;

	    if(prevSeg != currentSeg)
	    {
		System.out.println( prev + ".\t" + total);
		total=0;
	    }

	    total += ht.get(eachObj);

	    prevSeg = currentSeg;
	    prev = each;
	    
	}

	System.out.println( arr[arr.length-1] + ".\t" + total);
	
    }
}
