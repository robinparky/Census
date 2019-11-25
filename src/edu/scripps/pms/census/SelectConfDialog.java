/*
 * SelectConfDialog.java
 *
 * Created on May 2, 2006, 5:05 PM
 */

package edu.scripps.pms.census;

/**
 *
 * @author  rpark
 */
public class SelectConfDialog extends javax.swing.JDialog {
    
    private RelExMainFrame mFrame;
    
    /** Creates new form SelectConfDialog */
    public SelectConfDialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
        
        this.mFrame = (RelExMainFrame)parent;
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    // <editor-fold defaultstate="collapsed" desc=" Generated Code ">//GEN-BEGIN:initComponents
    private void initComponents() {
        btnGroup = new javax.swing.ButtonGroup();
        jPanel1 = new javax.swing.JPanel();
        labelRadioBtn = new javax.swing.JRadioButton();
        labelFreeRadBtn = new javax.swing.JRadioButton();
        selectBtn = new javax.swing.JButton();
        cancelBtn = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle("Select Scheme");
        jPanel1.setBorder(javax.swing.BorderFactory.createTitledBorder("Select Scheme"));
        btnGroup.add(labelRadioBtn);
        labelRadioBtn.setSelected(true);
        labelRadioBtn.setText("Isotope Labeling Configuration");
        labelRadioBtn.setActionCommand("label");
        labelRadioBtn.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        labelRadioBtn.setMargin(new java.awt.Insets(0, 0, 0, 0));

        btnGroup.add(labelFreeRadBtn);
        labelFreeRadBtn.setText("Label Free Configuration");
        labelFreeRadBtn.setActionCommand("label free");
        labelFreeRadBtn.setBorder(javax.swing.BorderFactory.createEmptyBorder(0, 0, 0, 0));
        labelFreeRadBtn.setMargin(new java.awt.Insets(0, 0, 0, 0));

        org.jdesktop.layout.GroupLayout jPanel1Layout = new org.jdesktop.layout.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel1Layout.createSequentialGroup()
                .add(20, 20, 20)
                .add(jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(labelFreeRadBtn)
                    .add(labelRadioBtn))
                .addContainerGap(24, Short.MAX_VALUE))
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .add(labelRadioBtn)
                .add(14, 14, 14)
                .add(labelFreeRadBtn)
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        selectBtn.setText("Select");
        selectBtn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                selectBtnActionPerformed(evt);
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
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(layout.createSequentialGroup()
                        .addContainerGap()
                        .add(jPanel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(layout.createSequentialGroup()
                        .add(85, 85, 85)
                        .add(selectBtn)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(cancelBtn)))
                .addContainerGap(42, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(layout.createSequentialGroup()
                .addContainerGap()
                .add(jPanel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(selectBtn)
                    .add(cancelBtn))
                .addContainerGap(15, Short.MAX_VALUE))
        );
        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void selectBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_selectBtnActionPerformed
// TODO add your handling code here:

        cancelBtnActionPerformed(evt);
        
        if( "label".equals(this.btnGroup.getSelection().getActionCommand()) )
        {
            LabelingConfDialog dialog = new LabelingConfDialog(this.mFrame,true);                  
            dialog.pack();
            dialog.setLocationRelativeTo(this);
            dialog.setVisible(true);
            dialog.setResizable(false);                
        }
        else //label free.
        {
            LabelFreeConfDialog dialog = new LabelFreeConfDialog(this.mFrame,true);                  
            dialog.pack();
            dialog.setLocationRelativeTo(this);
            dialog.setVisible(true);
            dialog.setResizable(false);    
            
        }
               
    }//GEN-LAST:event_selectBtnActionPerformed

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
                new SelectConfDialog(new javax.swing.JFrame(), true).setVisible(true);
            }
        });
    }
    
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup btnGroup;
    private javax.swing.JButton cancelBtn;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JRadioButton labelFreeRadBtn;
    private javax.swing.JRadioButton labelRadioBtn;
    private javax.swing.JButton selectBtn;
    // End of variables declaration//GEN-END:variables
    
}