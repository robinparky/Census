/*
 * MRMFragPanel.java
 *
 * Created on January 2, 2007, 10:37 AM
 */

package edu.scripps.pms.census;

import javax.swing.table.*;
import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import org.jdom.*;
import edu.scripps.pms.census.io.*;

import edu.scripps.pms.census.model.ChroProtein;
import edu.scripps.pms.census.model.ChroPeptide;
import edu.scripps.pms.census.model.ChroData;

import java.text.DecimalFormat;

/**
 *
 * @author  rpark
 */
public class MRMFragPanel extends javax.swing.JPanel {
    
    private DecimalFormat format = new DecimalFormat("0.00");
    /** Creates new form MRMFragPanel */
    public MRMFragPanel() {
        initComponents();
        initAdditional();
    }
    
    public MRMFragPanel(String chroFileName) throws IOException, JDOMException, Exception {
        this();
        
        ChroXmlReader cr = new ChroXmlReader(chroFileName, false);          
        
        java.util.List list = cr.getMrmCrvProteinList();

        Object[] simpleProteinArr = new Object[8];                        
        for(Iterator<ChroProtein> itr = list.iterator(); itr.hasNext(); )
        {
            ChroProtein protein = itr.next();
            //proteinTableModel.addRow(protein.getProteinData());

            simpleProteinArr[0] = protein.getLocus();
            simpleProteinArr[1] = "";            
            simpleProteinArr[2] = "";            
            simpleProteinArr[3] = "";            
            simpleProteinArr[4] = "";            
            simpleProteinArr[5] = ""; 
            simpleProteinArr[6] = ""; 
            simpleProteinArr[7] = protein.getDescription();

            proteinMRMTableModel.addRow( simpleProteinArr );
            
            for(Enumeration<TableColumn> enu = proteinSimpleTable.getColumnModel().getColumns(); enu.hasMoreElements(); )
            {
                TableColumn col = enu.nextElement();
                col.setCellRenderer(new FontRenderer());                
            }
            
            java.util.List pepList = protein.getPeptideList();
            
            for(Iterator<ChroPeptide> pepItr = pepList.iterator(); pepItr.hasNext(); )
            {
                ChroPeptide peptide = pepItr.next();
                
                ChroData cData = (ChroData)peptide.getDataList().get(0);
                
                long[] bionSamArr = cData.getBsIntensity();
                long[] yionSamArr = cData.getYsIntensity();

                peptide.getSequence();
                Vector tmpVec = new Vector();
                tmpVec.add("  " + peptide.getSequence());
                
                
                //peptide.get

                int charge = peptide.getChargeState();
                double d = Double.parseDouble( peptide.getMhPlus() );                
                double preCursor = (d+(charge-1)*CensusConstants.PROTON_MASS)/charge;
                        
                //System.out.println("--- " + peptide.getMhPlus() + " " + peptide.getChargeState() + " " + preCursor);
                
                boolean isPBion=false;
                boolean isMaxBion=false;
                
                int afterPIndex = -1; //ions bigger m/z than precursor
                int maxIonIndex = -1;
                double[] bStartMass = cData.getBsStartMass();
                double[] yStartMass = cData.getYsStartMass();
                long maxValue = 0;
                long maxPValue = 0;
                
                
                //System.out.println(peptide.getSequence() + " " + preCursor + " " + charge + " " + peptide.getMhPlus() + " " + maxValue + " " + maxPValue);
                
                for(int i=0;i<bionSamArr.length;i++)
                {
                    //System.out.println((i+1) + "bb " + bStartMass[i] + " " + preCursor + " " + bionSamArr[i] + " " + maxValue + " " + (bionSamArr[i]>maxValue) + " " + maxPValue);
                    
                    if(bionSamArr[i]>maxValue)
                    {                       
                        maxValue = bionSamArr[i];
                        maxIonIndex = (i+1);
                        //System.out.println(maxIonIndex);
                        isMaxBion = true;
                        if(bStartMass[i]>preCursor && bStartMass[i]>maxPValue)
                        {                            
                            maxPValue = bionSamArr[i];
                            isPBion = true;
                            afterPIndex = i;
                            
                      //      System.out.println("-->>" + maxPValue);
                        }
                    }
                }
                
                for(int i=0;i<yionSamArr.length;i++)
                {
                    //System.out.println( (yionSamArr.length-i) + "yy " + yStartMass[i] + " " + preCursor + " " + yionSamArr[i] + " " + maxValue + " " + (yionSamArr[i]>maxValue) + " " + maxPValue);
                    if(yionSamArr[i]>maxValue)
                    {                        
                        maxValue = yionSamArr[i];
                        maxIonIndex = (yionSamArr.length-i);
                        
                        //System.out.println(maxIonIndex);
                        
                        isMaxBion = false;
                        
                        //System.out.println( yStartMass[i] + " " + preCursor + " " +  yStartMass[i] + " " + maxPValue );
                        if(yStartMass[i]>preCursor && yionSamArr[i]>maxPValue)
                        {
                            maxPValue = yionSamArr[i];                    
                            isPBion = false;
                            afterPIndex = i;    
                            
                          //  System.out.println("-->>" + maxPValue);
                        }                        
                    }
                }
                
                double[] precurStartMass = null;
                double[] precurEndMass = null;
                double[] maxStartMass = null;
                double[] maxEndMass = null;
                
                if(isPBion)
                {
                    precurStartMass = cData.getBsStartMass();
                    precurEndMass = cData.getBsEndMass();
                }
                else
                {
                    precurStartMass = cData.getYsStartMass();
                    precurEndMass = cData.getYsEndMass();
                }
                
                
                if(isMaxBion)
                {
                    maxStartMass = cData.getBsStartMass();
                    maxEndMass = cData.getBsEndMass();
                }
                else
                {
                    maxStartMass = cData.getYsStartMass();
                    maxEndMass = cData.getYsEndMass();
                }

                String fragPIon = (isPBion)?" b" + (afterPIndex+1):" y" + (yionSamArr.length-afterPIndex);
                String fragMaxIon = (isMaxBion)?" b" + (maxIonIndex+1):" y" + maxIonIndex;
                    
                //System.out.println(maxValue + " " + index + " " + isBion);
                if(afterPIndex >= 0)  //found ions bigger than precursor ions
                {                    
                    
                    //tmpVec.add( ((isPBion)?" b":" y") + (afterPIndex+1) );                                        
                    tmpVec.add( fragPIon );                                        
                    tmpVec.add( " " + format.format(precurStartMass[afterPIndex]));                                                            
                    //tmpVec.add( ((isMaxBion)?" b":" y") + (maxIonIndex+1) );                    
                    tmpVec.add( fragMaxIon );                    
                    tmpVec.add( " " + format.format(maxStartMass[maxIonIndex-1])); //+endMass[maxIonIndex])/2) );                                
                    tmpVec.add( " " + charge );
                    tmpVec.add( " " + format.format(preCursor));
                    //tmpVec.add( " " + format.format(precurStartMass[afterPIndex]) + " - " + format.format(precurEndMass[afterPIndex]) );
                    tmpVec.add( " " );
                }
                else
                {
                    //Protein/Peptide", "Frag Ion", "Max mass after Precursor", "charge state", "Precursor", "Description", };
                    tmpVec.add(" N/A");                                        
                    tmpVec.add(" N/A");                    
                    //tmpVec.add( ((isMaxBion)?" b":" y") + (maxIonIndex+1) );                    
                    tmpVec.add( fragMaxIon );                    
                    tmpVec.add( " " + format.format(maxStartMass[maxIonIndex-1])); //+endMass[maxIonIndex])/2) );
                    tmpVec.add( " " + charge );
                    tmpVec.add( " " + format.format(preCursor));
                    //tmpVec.add( " " + format.format(maxStartMass[maxIonIndex]) + " - " + format.format(maxEndMass[maxIonIndex]) );
                    tmpVec.add( " " );                    
                }
                                
                proteinMRMTableModel.addRow(tmpVec);                
            }
        }
    }
    

                
    private javax.swing.JTable proteinSimpleTable;
    //private JScrollPane proteinSimplePane = new JScrollPane();
    
    private DefaultTableModel proteinMRMTableModel = 
        new DefaultTableModel( CensusConstants.PROTEIN_MRM_FRAG_COLUMNS, 0 ) 
    {
        public boolean isCellEditable(int row, int column)
        {
            return false;
        }        
    };
    
    public void initAdditional()
    {
        
        
        proteinSimpleTable = new JTable(proteinMRMTableModel) {
            
            /*
            public Component prepareRenderer(TableCellRenderer renderer, int rowIndex, int columnIndex) {
                /*
                Component c = super.prepareRenderer(renderer, rowIndex, columnIndex);
                if (rowIndex % 2 == 1 && !isCellSelected(rowIndex, columnIndex)) {
                        c.setFont(getFont());
                        c.setBackground(new Color(255, 255, 204)); // a light yellow
                } else if (isCellSelected(rowIndex, columnIndex)) {
                        // do nothing to allow operating system to handle shading
                } else {
                        c.setFont(getFont());
                        c.setBackground(new Color(255, 255, 255)); // a light yellow
                }
                return c;
            }  
            */
        };    
        
        this.proteinMRMPane.setViewportView(proteinSimpleTable);
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        proteinMRMPane = new javax.swing.JScrollPane();

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(proteinMRMPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 357, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(proteinMRMPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 257, Short.MAX_VALUE)
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents
    
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JScrollPane proteinMRMPane;
    // End of variables declaration//GEN-END:variables
    
}

class FontRenderer extends DefaultTableCellRenderer implements TableCellRenderer
    {
        //change them to whatever you like
        java.awt.Font font1 = new java.awt.Font("Arial",Font.BOLD,10); 
        java.awt.Font font2 = new java.awt.Font("Arial",Font.PLAIN,10);
        
                
        public FontRenderer()
        {
            //setFont(jTable1.getFont());
           // setFont(font1);
        }
        
        public java.awt.Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) 
        {
              java.awt.Component c = super.getTableCellRendererComponent(table,value,isSelected,hasFocus,  row, column);
 
              
              if(null == value)
                  return c;
              
              if(!value.toString().startsWith(" "))
              {
            //check for First Column
            //if(table.convertIndexToModel(column) == 0) 
                 c.setFont(font1); //or whatever font you want                                  
                 c.setBackground(new Color(182, 208, 255)); 
                         
              }
              else
              {
                  c.setFont(font2);
                  c.setBackground(new Color(255, 255, 255));
              }

            if (isSelected)
                    setBackground( new Color(255, 127, 127) );
              
            //else
              //   this.setFont(font2); //or this.setFont(table.getFont());
 
            return c;
        }
        

        
    };

