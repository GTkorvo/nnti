// @HEADER
// ***********************************************************************
// 
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

package Jpetra;

/*
 * VecView.java
 * Visualization class for Jpetra Vectors
 *
 * Created on February 3, 2002
 * Last modified on March 19, 2002
 *
 */

/**
 * @author  kahansen
 * @version
 */

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.lang.*;

public class VecView {
    
    public final static double BORDERFACTOR = 1.05;
    
    private double[] vector;
    private double minx, maxx, miny, maxy, xrange, yrange;
    private java.util.Vector data = new java.util.Vector();
    private Jpetra.Vector JpetraVector;
    String title = "";
    private Color histColor = Color.red;
    private Color xyColor = Color.blue;
    private int widthOfFrame = 400;
    private int heightOfFrame = 400;
    
    /** Creates new VecView object using values of double[] 
    
        @param vals the array of doubles
    */
    public VecView(double[] vals) {
	vector = vals;
	initialize(vector);
    }
    
    
    public VecView(int[] vals) {
        vector = new double[vals.length];
    	for (int k=0; k<vals.length; k++) {
		vector[k] = (double)vals[k];
	}	
	initialize(vector);
    }
	
    /** Creates new VecView object using values of double[] 
    
        @param vals the array of doubles
	@param vectorName the name of the vector to be displayed 
                          on the title bar
    */
    public VecView(double[] vals, String vectorName) {
	vector = vals;
	title = vectorName;
	initialize(vector);
    }
    
    /** Creates new VecView object using values of double[] 
    
    	@param vals the name of the Jpetra.Vector
    */
    public VecView(Jpetra.Vector vals) {
	vector = vals.extractView()[0];
	initialize(vector);
    }
    
    
    /** Creates new VecView object using values of double[] 
    
    	@param vals the name of the JpetraVector
	@param vectorName the name of the vector to be displayed 
                          on the title bar
    */
    public VecView(Jpetra.Vector vals, String vectorName) {
	vector = vals.extractView()[0];
	title = vectorName;
	initialize(vector);
    }
    
    
    /** Draws a histogram of the Jpetra vector */
    public void histogram() {
	
	//Computes X and Y range, scales frame and plots data
	class GraphPane extends JPanel {
	    public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Dimension s = getSize();
		
		// Compute scale factors,
		// Leave space for axis lines and labels
		final float xfact = (s.width-45) / (float)xrange;
		final float yfact = (s.height-25) / (float)(maxy*BORDERFACTOR);
		
		// Scale and plot the data
		for (int i=0; i<vector.length; i++) {
		    Apoint d = (Apoint)data.elementAt(i);
		    float x = (float)(d.x)*(float)xfact;
		    float y = (float)(d.y)*(float)yfact;
		    
		    g.setColor(histColor);
		    
		    // Draw a rectangle using (x, y), width is scaled 
		    //   depending on the number of values
		    // Fill rectangle down to the x-axis line
		    g.fillRect(((int)x+45), s.height-(int)y-25, 
		                    (int)((.6*s.width)/maxx), (int)y);
		    
		    // Draw x axis line, label every 5th number
		    g.setColor(Color.black);
		    Font axisLabel = new Font("Serif", Font.BOLD, 8);
		    g.setFont(axisLabel);
		    g.drawLine((int)minx+40, s.height-25, 
		                           s.width ,s.height-25);
		    for (int j=0; j<vector.length; j=j+5) {
			if (i == j)
			    g.drawString(""+i+"", (int)x+45, s.height-7);
		    }
		}
		// Draw and label y axis line
		yAxisHistogram(g, s);
	    }
	}

	final GraphJFrame f = new GraphJFrame();
	f.setTitle(title);
	GraphPane graph = new GraphPane();
	f.getContentPane().add(graph);
	f.setVisible(true);
	
    }
    
    /** Draw an xy graph of the Jpetra vector */
    public void xygraph() {
	
	// integer arrays for input into 
	// the drawPolyline function of Graphpane
	final int[] int_y = new int[vector.length];
	final int[] int_x = new int[vector.length];
		
	//Computes X and Y range, scales frame and plots data
	class GraphPane extends JPanel {
	    public void paintComponent(Graphics g) {
		super.paintComponent(g);
		Dimension s = getSize();
		
		// Compute scale factors
		// Leave space for axis lines and labels
		final float xfact = (s.width-45) / (float)xrange;
		final float yfact = (s.height-25) / (float)yrange;
		
		// Scale and plot the data
		for (int i=0; i<vector.length; i++) {
		    Apoint d = (Apoint)data.elementAt(i);
		    float x = (float)(d.x)*(float)xfact;
		    float y = (float)(d.y-(miny))*(float)yfact;
		    
		    g.setColor(xyColor);
		    
		    // Draw a 3-pixel rectangle centered, 
		    //   so -1 both x and y.
		    // AWT numbers y from 0 down, so invert:
		    g.fillRect(((int)x)-1+45, s.height-1-(int)y-25, 3, 3);
		   
		    // fill the integer arrays, x_int & y_int
		    int_x[i] = (int)x+45;
		    int_y[i] = s.height-(int)y-25;
		    
		    // Draw x axis line, label every 5th value
		    g.setColor(Color.black);
		    Font axisLabel = new Font("Serif", Font.BOLD, 8);
		    g.setFont(axisLabel);
		    g.drawLine((int)minx+40, s.height-25, 
		                           s.width ,s.height-25);
		    for (int j=0; j<vector.length; j=j+5) {
			if (i == j)
			    g.drawString(""+i+"", (int)x+45, s.height-7);
		    }
		}
		
		// Draw and label y axis line
		yAxis(g, s);
		
		// Draw a line between the points on the graph
		g.setColor(xyColor);
		g.drawPolyline(int_x, int_y, vector.length);		
	    }
	}
	
	final GraphJFrame f = new GraphJFrame();
	f.setTitle(title);
	GraphPane graph = new GraphPane();
	f.getContentPane().add(graph);
	f.setVisible(true);
	
    }
    
    public void setColor(Color col) {
        histColor = col;
	xyColor = col;
    }
    
    public void setSize(int width, int height) {
        widthOfFrame = width;
	heightOfFrame = height;
    }	
    
    
    /** Creates the skeleton frame of any VecView graph */
    class GraphJFrame extends JFrame implements ActionListener {
	
	JMenuBar GraphMenuBar;
	JMenu FileMenu, HelpMenu;
	JMenuItem StatisticsItem,CloseVecVisItem,OpenHelpItem;
	
	public void main(String [] Args) {
	    GraphJFrame G=new GraphJFrame();
	    G.show();
	}
	
	public GraphJFrame() {
	    // Set up the menu bar
	    GraphMenuBar=new JMenuBar();
	    FileMenu=new JMenu("File");
	    HelpMenu=new JMenu("Help");
	    GraphMenuBar.add(FileMenu);
	    GraphMenuBar.add(HelpMenu);
	    setJMenuBar(GraphMenuBar);
	    StatisticsItem=new JMenuItem("Statistics");
	    StatisticsItem.addActionListener(this);
	    CloseVecVisItem=new JMenuItem("Close VecVis");
	    CloseVecVisItem.addActionListener(this);
	    FileMenu.add(StatisticsItem);
	    FileMenu.add(CloseVecVisItem);
	    OpenHelpItem=new JMenuItem("Open Help File");
	    OpenHelpItem.addActionListener(this);
	    HelpMenu.add(OpenHelpItem);
	    
	    Container Pane=getContentPane();
	    
	    setSize(widthOfFrame,heightOfFrame);
	    addWindowListener(new WindowAdapter() {
		public void windowClosing(WindowEvent e)
		{  System.exit(0);
		}
	    } );
	}
	
	
	public void actionPerformed(ActionEvent e) {
	    Object source=e.getSource();
	    
	    // Show Statistics Box
	    if(source==StatisticsItem) {
		Font statsFont = new Font("Dialog", Font.PLAIN, 12);
		String output = ("Vector length: "+vector.length+"\n");
		output += ("Minimum Value: " + miny + "\n");
		output += ("Maximum Value: " + maxy + "\n");
		
		// Access methods of Jpetra MultiVector class for 
		// 1-, 2-, and infinity-norms, also median value
		
		double [] norm1 = new double[1];
		double [] norm2 = new double[1];
		double [] normI = new double[1];
		double [] mean = new double[1];
		JpetraVector.norm1(norm1);
		JpetraVector.norm2(norm2);
		JpetraVector.normInf(normI);
		JpetraVector.meanValue(mean);
		output += ("Norm: " + norm1[0] + "\n");
		output += ("2-norm: " + norm2[0] + "\n");
		output += ("Infinity norm: " + normI[0] + "\n");
		output += ("Mean value: " + mean[0] + "\n");
		JTextArea stats = new JTextArea(output, 8, 8);
		stats.setEditable(false);
		stats.setFont(statsFont);
		JOptionPane.showMessageDialog(null, stats, 
		         "Vector Statistics", JOptionPane.PLAIN_MESSAGE);
	    }
	    
	    // Exit VecVis
	    else if(source==CloseVecVisItem) {
		System.exit(0);
	    }
	    
	    // Open Help File
	    else if(source==OpenHelpItem) {
		Font helpFont = new Font("Monospaced", Font.PLAIN, 12);
		String output = ("Using VecVis" + "\n" + "\n");
		output += ("Use this tool to visualize Jpetra vectors " +
		           "that you are working with." + "\n");
		output += ("Find statistics of the Jpetra vector using " +
		           "File --> Statistics." + "\n");
		output += ("\n" + "\n" + "\n");
		output += ("Created by: Kara Hansen");
		JTextArea help = new JTextArea(output, 8, 12);
		help.setEditable(false);
		help.setFont(helpFont);
		JOptionPane.showMessageDialog(null, help, "Help", 
					      JOptionPane.PLAIN_MESSAGE);
	    }
	}
    }
    
    class Apoint {
	float x;
	float y;
	public String toString() {
	    return "Apoint("+x+", "+y+")";
	}
    }
    
    /** Draws and labels the y-axis line for an xy graph */
    private void yAxis(Graphics g, Dimension s) {
	g.drawLine((int)minx+40, s.height-25, (int)minx+40, 0);
	
	// draw evenly spaced text labels, rounded to 4 decimals
	for (double d=0.0; d<1.0; d=d+.05) {
	    drawYValue(d, g, s);
	}
    }
    
    /** Draws and labels the y-axis line for a histogram */
    private void yAxisHistogram(Graphics g, Dimension s) {
	g.drawLine((int)minx+40, s.height-25, (int)minx+40, 0);
	
	// draw evenly spaced text labels, rounded to 4 decimals
	for (double d=0.0; d<1.0; d=d+.05) {
	    drawYValueHistogram(d, g, s);
	}
    }
    
    /** Returns the smallest value of the Jpetra vector */
    private double smallest(double[] vals) {
	double miny = vector[1];
	for (int i=0; i<vector.length; i++)
	    if (vector[i] < miny) miny = vector[i];
	return miny;
    }
    
    /** Returns the largest value of the Jpetra vector */
    private double largest(double[] vals) {
	double maxy = vector[1];
	for (int i=0; i<vector.length; i++)
	    if (vector[i] > maxy) maxy = vector[i];
	return maxy;
    }
    
    /** Initializes values for the constructor */
    private void initialize(double[] vector) {
	minx = 0;
	maxx = vector.length;
	miny = smallest(vector);
	maxy = largest(vector);
	xrange = (maxx - minx);
	yrange = (maxy - miny) * BORDERFACTOR;
	
	// Read the vector of doubles, convert into a vector of points
	for (int i=0; i<vector.length; i++) {
	    Apoint d = new Apoint();
	    d.x = i;
	    d.y = (float)vector[i];
	    data.add(d);
	}
	
	// Make the array of doubles into a Jpetra.Vector
	try {
		Jpetra.SerialComm comm = new Jpetra.SerialComm();
		Jpetra.Map map = new Jpetra.Map(1, 0, comm);
     		double [][] JpetraVecArray = new double [1][]; // Jpetra vector
							       //constructor needs a double [][]
    		JpetraVecArray[0] = vector; // Set first and only array of JpetraVec to
					       // point to vecViewVec
     		Jpetra.Vector JpetraVector = new Jpetra.Vector("Copy", map, JpetraVecArray);
	} catch (Jpetra.JpetraException f) {
			System.out.println(f);
      	}
    }
    
    private void drawYValue(double percent, Graphics g, Dimension s) {
	g.drawString(""+(double)((int)((percent*yrange+miny)*10000))/10000, 
	            (int)minx+5, (int)((1-percent)*(s.height-25)));
    }
    
    private void drawYValueHistogram(double percent, Graphics g, Dimension s) {
	g.drawString(""+(double)((int)((percent*(maxy*BORDERFACTOR))*10000))/10000, 
	            (int)minx+5, (int)((1-percent)*(s.height-25)));
    }
}
