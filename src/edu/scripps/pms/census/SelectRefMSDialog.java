/*
 * SelectRefMSDialog.java
 *
 * Created on July 24, 2006, 10:17 AM
 */

package edu.scripps.pms.census;

import javax.swing.*;
import java.io.*;
import org.jdom.*;
import org.jdom.output.*;

import java.util.*;
import edu.scripps.pms.census.model.*;

/**
 *
 * @author  rpark
 */
public class SelectRefMSDialog extends javax.swing.JDialog {
    
    private Element rootEle;
    private Vector<SelectFileModel> fileList;
    
    /** Creates new form SelectRefMSDialog */
    public SelectRefMSDialog(java.awt.Frame parent, boolean modal, Element rootEle, Vector<SelectFileModel> fileList) {
        super(parent, modal);
        this.rootEle = rootEle;
        
        initComponents();
        this.fileList = fileList;                
        
        fileBox.setModel(new javax.swing.DefaultComboBoxModel(fileList));
        
        //fileBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 4" }));
        
        
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        jLabel1 = new javax.swing.JLabel();
        okBtn = new javax.swing.JButton();
        cancelBtn = new javax.swing.JButton();
        fileBox = new javax.swing.JComboBox();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        jLabel1.setText("Select Reference Spectra File");

        okBtn.setText("OK");
        okBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okBtnActionPerformed(evt);
            }
        });

        cancelBtn.setText("Cancel");
        cancelBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelBtnActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                            .add(fileBox, 0, 423, Short.MAX_VALUE)
                            .addContainerGap())
                        .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                            .add(jLabel1)
                            .add(135, 135, 135)))
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                        .add(okBtn)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(cancelBtn)
                        .add(151, 151, 151))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(jLabel1)
                .add(22, 22, 22)
                .add(fileBox, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED, 24, Short.MAX_VALUE)
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(okBtn)
                    .add(cancelBtn))
                .addContainerGap())
        );
        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void okBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okBtnActionPerformed
// TODO add your handling code here:
        setVisible(false);
        dispose();
        
        SelectFileModel fModel = (SelectFileModel)fileBox.getSelectedItem();
        Element refEle = new Element("ref");
        refEle.addContent( new Element("sample_name").addContent(fModel.getSampleName()) );
        refEle.addContent( new Element("file_name").addContent(fModel.getSpectraFileName()) );
        rootEle.addContent(refEle);
        
        
        JFileChooser choose = new JFileChooser();
        choose.setMultiSelectionEnabled(false);
        choose.setDialogTitle("Save Configuration file");
        choose.setApproveButtonText("Create");
               
        //choose.addChoosableFileFilter( new SimpleFileNameFilter("txt", "CenSus Report File (*.txt)") );
                
        File f = new File("census_config.xml");
        choose.setSelectedFile(f);        
        
        
        int returnVal = choose.showOpenDialog(this);        

        if(returnVal == choose.CANCEL_OPTION)
            return;
                
        File file = choose.getSelectedFile();       
        
        
        
        try
        {
            Document doc = new Document(rootEle);
            OutputStream os = new FileOutputStream(file); //(filePath + "census_chro.xml");             
            XMLOutputter outputter = new XMLOutputter();
            outputter.setFormat(Format.getPrettyFormat());             
            outputter.output(doc, os);
            os.close();         
        } 
        catch (IOException e)
        {
            JOptionPane.showMessageDialog(this, "Error: " + e, "Error", JOptionPane.ERROR_MESSAGE);
        }

//        JOptionPane.showMessageDialog(this, "configuration file was successfully saved", "Saved", JOptionPane.INFORMATION_MESSAGE);

    }//GEN-LAST:event_okBtnActionPerformed

    private void cancelBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelBtnActionPerformed
// TODO add your handling code here:
        setVisible(false);
        dispose();
    }//GEN-LAST:event_cancelBtnActionPerformed
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new SelectRefMSDialog(new javax.swing.JFrame(), true, null, null).setVisible(true);
            }
        });
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton cancelBtn;
    private javax.swing.JComboBox fileBox;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JButton okBtn;
    // End of variables declaration//GEN-END:variables
    
}