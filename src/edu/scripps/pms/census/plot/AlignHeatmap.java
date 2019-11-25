package edu.scripps.pms.census.plot;

import java.awt.*;
import javax.swing.*;
import java.awt.event.*;
import java.io.*;

public class AlignHeatmap extends JFrame 
{

    //JPanel  pan;

    JScrollPane pan;

    public AlignHeatmap(String fileName) throws Exception
    {
	
	long start;
	start = System.currentTimeMillis();

	JPanel topPanel = new JPanel();
	topPanel.setLayout( new BorderLayout() );
	getContentPane().add( topPanel );
	
	pan  = new myPan(fileName, start);

	Icon image = new ImageIcon( "dummy.jpg" );
	JLabel label = new JLabel( image );
	pan.getViewport().add( label );


//	setBounds(10,10,100,500);
	setSize(950, 950);
//	getContentPane().add(pan,BorderLayout.CENTER);
	topPanel.add( pan, BorderLayout.CENTER );
	setVisible(true);
    }

    public class myPan extends JScrollPane 
    {
	private String fileName;
	private long start;
	boolean runIt = false;

	public myPan(String fileName, long start)
	{
	    super();
	    this.fileName = fileName;
	    this.start = start;
	}

	public void repaint(long tm,
		int x,
		int y,
		int width,
		int height)
	{

	}

	public void paint(Graphics g)
	{
	/*
	    if(runIt)
		return;

	    runIt = true;
	  */  
	    super.paint(g);

	    
	    int min=-1;
	    int max=1;
	    double mid = (min+max)/2;

	    
	    try 
	    {

		BufferedReader br = new BufferedReader(new FileReader(fileName));

		String line;
		int index=0;
                int  tick=5;

		int y=0;

		System.out.println("start");

		while ((line = br.readLine()) != null) {
		    index++;
		    String[] arr = line.split("\t");

		    if(index%tick != 0)
			continue;

		    y++;
		    System.out.print(".");
		    int x=0;
		    for(int i=0;i<arr.length;i++)
		    {

			if(i%tick != 0)
			    continue;

			x++;
			double input;
			if("".equals(arr[i]) || null == arr[i])
			{
			    input = 0;
                            System.out.print("arr[i] = null\ti = " + i);
			    continue;
			}
	
			//input = (int)Math.log10(Double.parseDouble(arr[i]));
			//input = Math.log10(Double.parseDouble(arr[i]));
			input = Double.parseDouble(arr[i]);
                   
			//if(input>-1)
                            //System.out.println(input);
			if(input<=min)
			{
			    g.setColor(new Color(0, 255, 0));//green
			    g.fillRect(x, y, 1, 1);
			    continue;
			}
			
			if(input>=max)
			{
			    g.setColor(new Color(255, 0, 0));//red
			    g.fillRect(x, y, 1, 1);
			    continue;
			}
			
			if(input<=mid) //shade of green color
			{
			    g.setColor(new Color(0, (int)((mid-input)/(max-mid)*255), 0));
			    g.fillRect(x, y, 1, 1);
			}
			else //shade of Red color
			{
			    g.setColor(new Color((int)((input-mid)/(max-mid)*255), 0, 0));
			    g.fillRect(x, y, 1, 1);
			}

		    }
		}

	    } catch (Exception e) {System.out.println(e); }
	System.out.println( System.currentTimeMillis() - start );
	}

    }

    public static void main (String[] args) throws Exception
    {
	new AlignHeatmap(args[0]);

    }
}




