/*
 * ColorTable.java
 *
 * Created on May 4, 2006, 11:25 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
 */

package edu.scripps.pms.census;

import javax.swing.JTable;
import javax.swing.SwingUtilities;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.table.TableCellRenderer;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.Component;
import java.awt.Color;

/**
 *
 * @author rpark
 */
public class ColorTable extends JTable {
    
    /** Creates a new instance of ColorTable */
    public ColorTable() {        
       super();
    }
    
    public Component prepareRenderer(
                            TableCellRenderer renderer,
                            int row, int col) {
       Component c = super.prepareRenderer(renderer,
                                  row, col);
       if (row % 2 == 1 && !isCellSelected(row,col)) {
         c.setBackground(new Color(220, 255, 213));         
       } else {
         c.setBackground(getBackground());
         //c.setBackground(Color.WHITE);
       }
       return c;
     }    
}
