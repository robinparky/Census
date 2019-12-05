/*
 * RelaxMainFrame.java
 *
 * Created on March 18, 2005, 3:52 PM
 */
package edu.scripps.pms.census;

import edu.scripps.pms.census.dialog.*;
import javax.swing.*;
import java.io.*;
import java.net.URI;

import java.awt.Container;
import java.awt.Color;
import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;

import java.awt.event.*;
import javax.swing.table.*;

import java.util.*;
//import org.apache.commons.math.distribution.*;
import org.apache.commons.math3.distribution.*;

import edu.scripps.pms.census.util.SimpleFileNameFilter;
import edu.scripps.pms.census.io.parse.ChroXMLParser;
import edu.scripps.pms.census.io.*;

import edu.scripps.pms.census.hash.*;
import edu.scripps.pms.census.util.*;

import edu.scripps.pms.census.util.io.DTASelectFilterReader;

import edu.scripps.pms.census.model.*;

import javax.swing.table.DefaultTableModel;

import ptolemy.plot.*;
import org.jdom.*;
import org.jdom.input.SAXBuilder;

import ptolemy.plot.plotml.PlotBoxMLParser;
import ptolemy.plot.plotml.PlotMLParser;
import edu.scripps.pms.census.plot.*;

import edu.scripps.pms.census.util.LinearRegression;

import edu.scripps.pms.census.conf.PostOptions;
import edu.scripps.pms.census.conf.Configuration;

import edu.scripps.pms.census.chroalign.*;
import gnu.trove.TLongArrayList;
import gnu.trove.TDoubleArrayList;

import edu.scripps.pms.stats.STDev;

import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.stat.regression.*;
//import org.apache.commons.math3.MathException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import rpark.statistics.CommonStat;
import rpark.statistics.model.GaussianPeakModel;
//import rpark.statistics.Outlier;

/**
 *
 * @author Robin Park
 * @version $Id: RelExMainFrame.java,v 1.91 2014/09/09 19:26:32 rpark Exp $
 */
public class RelExMainFrame extends javax.swing.JFrame implements java.beans.PropertyChangeListener {

  protected static String version;
  protected boolean chroFileOpen = false;
  private static BufferedWriter bw = null;
  private static double MIN_RATIO_VALUE = 0.001;
  private static double MAX_RATIO_VALUE = 1000.0;
  private static double LOWER_BOUND_RATIO_THRESHOLD = 0.033;
  private static double UPPER_BOUND_RATIO_THRESHOLD = 30;

  /**
   * Creates new form RelaxMainFrame
   */
  public RelExMainFrame() {

    Configuration conf = Configuration.getInstance();
    version = conf.getVersion();
    initComponents();
    initAdditional();
    //       initTest();
  }

  public static void printHeader(PrintStream p) {
    Configuration conf = Configuration.getInstance();
    p.println("H\tCensus version " + conf.getVersion());
//	p.println("H\tCensus a merged file");
    p.println("H\tCreate by ");
    p.println("H\tRobin, Sung Kyu Park rpark@scripps.edu");
    p.println("H\tJohn Venable jvenable@gnf.org");
//	p.println("H\tMichael J. MacCoss");
    p.println("H\tThe Scripps Research Institute, La Jolla, CA");
    p.println("H\tOL stands for outlier");
    p.print("H\tcreated date\t");
    p.println(new Date());
  }

  private void initAdditional() {
    welcomeLabel.setText("Welcome to Census version " + conf.getVersion());

    proteinSimpleTable = new JTable(proteinSimpleTableModel) {
      public Component prepareRenderer(TableCellRenderer renderer, int rowIndex, int columnIndex) {
        Component c = super.prepareRenderer(renderer, rowIndex, columnIndex);
                /*
                 System.out.println("===>>" + columnIndex);
                 if (rowIndex == 0) {
                 c.setFont(new Font("Arial", Font.BOLD, 18));
                 } */
        if (rowIndex % 2 == 1 && !isCellSelected(rowIndex, columnIndex)) {
          c.setFont(getFont());
          c.setBackground(new Color(182, 208, 255));
        } else if (isCellSelected(rowIndex, columnIndex)) {
          // do nothing to allow operating system to handle shading
        } else {
          c.setFont(getFont());
          c.setBackground(new Color(255, 255, 255)); // a light yellow
        }
        return c;
      }
    };

    this.proteinSimplePane.setViewportView(proteinSimpleTable);

    peptideListTable = new javax.swing.JTable(peptideTableModel) {
      public Component prepareRenderer(TableCellRenderer renderer, int rowIndex, int columnIndex) {
        Component c = super.prepareRenderer(renderer, rowIndex, columnIndex);

        if (rowIndex % 2 == 1 && !isCellSelected(rowIndex, columnIndex)) {
          c.setFont(getFont());
          c.setBackground(new Color(182, 208, 255));
        } else if (isCellSelected(rowIndex, columnIndex)) {
          // do nothing to allow operating system to handle shading
        } else {
          c.setFont(getFont());
          c.setBackground(new Color(255, 255, 255)); // a light yellow
          // if not shaded, match the table's background
          //c.setBackground(getBackground());
          //c.setFont(getFont());
          //c.setForeground(getForeground());
        }
        return c;
      }
    };

    peptideList.setViewportView(peptideListTable);

    peptideListTable.addKeyListener(new java.awt.event.KeyAdapter() {
      public void keyPressed(java.awt.event.KeyEvent evt) {
        peptideListTableKeyPressed(evt);
      }
    });
    peptideListTable.addMouseListener(new java.awt.event.MouseAdapter() {
      public void mouseClicked(java.awt.event.MouseEvent evt) {
        peptideListTableMouseClicked(evt);
      }
    });

    proteinSimpleTable.addMouseListener(new java.awt.event.MouseAdapter() {
      public void mouseClicked(java.awt.event.MouseEvent evt) {
        proteinSimpleTableMouseClicked(evt);
      }
    });

    proteinSimpleTable.addKeyListener(new java.awt.event.KeyAdapter() {
      public void keyPressed(java.awt.event.KeyEvent evt) {
        proteinSimpleTableKeyPressed(evt);
      }
    });

    proteinSimplePane.setViewportView(proteinSimpleTable);

        /*
         int vColIndex = 1;
         TableColumn col = table.getColumnModel().getColumn(vColIndex);
         int width=399;
         col.setPreferredWidth(width);
         */
    TableColumn column = null;

    column = peptideListTable.getColumnModel().getColumn(0);
    column.setPreferredWidth(10);

    column = this.proteinTable.getColumnModel().getColumn(0);
    column.setPreferredWidth(10);

        /*
         for (int i = 0; i < peptideListTable.getColumnCount(); i++) {
         column = peptideListTable.getColumnModel().getColumn(i);

         column.setPreferredWidth(1);


         switch(i)
         {
         case 0:
         column.setPreferredWidth(1);
         break;

         case 1:
         column.setPreferredWidth(1);
         break;

         case 2:
         break;
         case 3:
         break;
         case 4:
         break;
         case 5:
         break;
         case 6:
         break;

         default :
         //column.setPreferredWidth(100);
         break;

         }

         }
         */
    corrPlot = new CorrelationPlot();
    corrPlot.setBackground(new Color(255, 255, 255));
    corrPlot.setSize(250, 250);

    PostOptions options = PostOptions.getInstance();
    options.setFilterFragmentIons(true);
    options.setDisplayFragmentIons(false);

    this.tabbedPanel.setVisible(false);

    proteinRatioDistPanel.setLayout(new java.awt.BorderLayout());
    pepTabbedPanel.addTab("Peptide Ratio Dist", proteinRatioDistPanel);

    tabbedPanel.setEnabledAt(1, false);
  }

  private void initTest() {
    tabbedPanel.setEnabledAt(1, false);

    ////remove this later.  this is for testing..
//        File file = new File("/home/rpark/rpark_on_data/relax_run/BDCHN15");
    File tempFile = new File("/home/rpark/rpark_on_data/relax_run/BDCHN15/N15isotope.ini");

    //File tempFile = new File("E:\\relex\\newFinal_data\\data-independent\\N15isotope.ini");
    isotopeFileField.setText(tempFile.getAbsolutePath());
    currentDirectory = "/home/rpark/rpark_on_data/relax_run/BDCHN15";
    //currentDirectory = "E:\\relex\\newFinal_data\\data-independent";
    this.reportItem.setEnabled(true);

    //String dtaselectFile="e:\\DTASelect.txt";
  }

  /**
   * This method is called from within the constructor to initialize the form.
   * WARNING: Do NOT modify this code. The content of this method is always
   * regenerated by the Form Editor.
   */
  // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
  private void initComponents() {

    extractRadioGrp = new javax.swing.ButtonGroup();
    quantModeRadioGrp = new javax.swing.ButtonGroup();
    jPanel1 = new javax.swing.JPanel();
    chroPanel = new javax.swing.JPanel();
    isotopeLabel = new javax.swing.JLabel();
    isotopeFileField = new javax.swing.JTextField();
    fileSelectBtn = new javax.swing.JButton();
    extractPanel = new javax.swing.JPanel();
    jRadioButton1 = new javax.swing.JRadioButton();
    jRadioButton2 = new javax.swing.JRadioButton();
    jPanel2 = new javax.swing.JPanel();
    fullMassScan = new javax.swing.JRadioButton();
    msmsScan = new javax.swing.JRadioButton();
    jPanel3 = new javax.swing.JPanel();
    jCheckBox1 = new javax.swing.JCheckBox();
    jLabel2 = new javax.swing.JLabel();
    scanBefore = new javax.swing.JTextField();
    scanAfter = new javax.swing.JTextField();
    jLabel3 = new javax.swing.JLabel();
    jPanel4 = new javax.swing.JPanel();
    jCheckBox2 = new javax.swing.JCheckBox();
    jTextField5 = new javax.swing.JTextField();
    jLabel4 = new javax.swing.JLabel();
    jCheckBox3 = new javax.swing.JCheckBox();
    jTextField6 = new javax.swing.JTextField();
    jLabel5 = new javax.swing.JLabel();
    extractBtn = new javax.swing.JButton();
    resetBtn = new javax.swing.JButton();
    chorNoteLabel1 = new javax.swing.JLabel();
    chroNoteLabel2 = new javax.swing.JLabel();
    peptideListBox = new javax.swing.JPanel();
    quanPanel = new javax.swing.JPanel();
    fragIonScrollPanel = new javax.swing.JScrollPane();
    fragIonPanel = new javax.swing.JPanel();
    jToolBar1 = new javax.swing.JToolBar();
    open = new javax.swing.JButton();
    save = new javax.swing.JButton();
    export = new javax.swing.JButton();
    filterBtn = new javax.swing.JButton();
    filePathLabel = new javax.swing.JLabel();
    tabbedPanel = new javax.swing.JTabbedPane();
    proteinPanel = new javax.swing.JPanel();
    proteinListPanel = new javax.swing.JScrollPane();
    proteinTable = new javax.swing.JTable(proteinTableModel) {
      public Component prepareRenderer(TableCellRenderer renderer, int rowIndex, int columnIndex) {
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
    };
    ;
    searchField = new javax.swing.JTextField();
    searchBtn = new javax.swing.JButton();
    peptidePanel = new javax.swing.JPanel();
    peptideList = new javax.swing.JScrollPane();
    proteinSimplePane = new javax.swing.JScrollPane();
    proteinInfoPanel = new javax.swing.JPanel();
    proteinLabel = new javax.swing.JLabel();
    jSplitPane = new javax.swing.JSplitPane();
    paramPanel = new javax.swing.JPanel();
    rrLabel = new javax.swing.JLabel();
    shiftField = new javax.swing.JTextField();
    shiftLabel = new javax.swing.JLabel();
    rrField = new javax.swing.JTextField();
    measuredRatioLabel = new javax.swing.JLabel();
    regressionRatioField = new javax.swing.JTextField();
    areaRatioField = new javax.swing.JTextField();
    measuredRatioLabel1 = new javax.swing.JLabel();
    regLnLabel = new javax.swing.JLabel();
    areaRatioLogField = new javax.swing.JTextField();
    pvalueField = new javax.swing.JTextField();
    pvalueLabel = new javax.swing.JLabel();
    jTextField1 = new javax.swing.JTextField();
    jLabel1 = new javax.swing.JLabel();
    regScoreField = new javax.swing.JTextField();
    regScoreLabel = new javax.swing.JLabel();
    pepTabbedPanel = new javax.swing.JTabbedPane();
    chromatogramPanel = new javax.swing.JPanel();
    correlationPanel = new javax.swing.JPanel();
    mainPanel = new javax.swing.JPanel();
    jLabel7 = new javax.swing.JLabel();
    welcomeLabel = new javax.swing.JLabel();
    jMenuBar1 = new javax.swing.JMenuBar();
    fileMenu = new javax.swing.JMenu();
    openItem = new javax.swing.JMenuItem();
    reportItem = new javax.swing.JMenuItem();
    jSeparator1 = new javax.swing.JSeparator();
    exitItem = new javax.swing.JMenuItem();
    runMenu = new javax.swing.JMenu();
    runItem = new javax.swing.JMenuItem();
    runNonLabel = new javax.swing.JMenuItem();
    confItem = new javax.swing.JMenuItem();
    toolMenu = new javax.swing.JMenu();
    optionItem = new javax.swing.JMenuItem();
    mergeItem = new javax.swing.JMenuItem();
    alignSpectra = new javax.swing.JMenuItem();
    openSpectra = new javax.swing.JMenuItem();
    mrmCsv = new javax.swing.JMenuItem();
    helpMenu = new javax.swing.JMenu();
    versionItem = new javax.swing.JMenuItem();

    chroPanel.setEnabled(false);
    chroPanel.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    isotopeLabel.setText("Select isotope file");
    isotopeLabel.setPreferredSize(new java.awt.Dimension(10, 14));
    chroPanel.add(isotopeLabel, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 10, 150, -1));
    chroPanel.add(isotopeFileField, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 30, 270, -1));

    fileSelectBtn.setLabel("...");
    fileSelectBtn.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        fileSelectBtnActionPerformed(evt);
      }
    });
    chroPanel.add(fileSelectBtn, new org.netbeans.lib.awtextra.AbsoluteConstraints(290, 30, -1, -1));

    extractPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Extraction Options"));
    extractPanel.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    extractRadioGrp.add(jRadioButton1);
    jRadioButton1.setSelected(true);
    jRadioButton1.setText("Single Experiment");
    extractPanel.add(jRadioButton1, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 30, -1, -1));

    extractRadioGrp.add(jRadioButton2);
    jRadioButton2.setText("Multiple Experiment");
    jRadioButton2.setEnabled(false);
    extractPanel.add(jRadioButton2, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 60, -1, -1));

    chroPanel.add(extractPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 80, 180, 110));

    jPanel2.setBorder(javax.swing.BorderFactory.createTitledBorder("Quantitation Mode Options"));
    jPanel2.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    quantModeRadioGrp.add(fullMassScan);
    fullMassScan.setSelected(true);
    fullMassScan.setText("Full Mass Scans");
    fullMassScan.setActionCommand("f");
    jPanel2.add(fullMassScan, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 30, -1, -1));

    quantModeRadioGrp.add(msmsScan);
    msmsScan.setText("MS/MS spectra");
    msmsScan.setActionCommand("m");
    jPanel2.add(msmsScan, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 60, -1, -1));

    chroPanel.add(jPanel2, new org.netbeans.lib.awtextra.AbsoluteConstraints(10, 210, 180, 110));

    jPanel3.setBorder(javax.swing.BorderFactory.createTitledBorder("Time Width Options"));
    jPanel3.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    jCheckBox1.setSelected(true);
    jCheckBox1.setText("Extract by number of MS scans");
    jPanel3.add(jCheckBox1, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 30, -1, -1));

    jLabel2.setText("scans before");
    jPanel3.add(jLabel2, new org.netbeans.lib.awtextra.AbsoluteConstraints(280, 30, -1, -1));

    scanBefore.setText("50");
    jPanel3.add(scanBefore, new org.netbeans.lib.awtextra.AbsoluteConstraints(240, 30, -1, -1));

    scanAfter.setText("50");
    jPanel3.add(scanAfter, new org.netbeans.lib.awtextra.AbsoluteConstraints(240, 60, -1, -1));

    jLabel3.setText("scans after");
    jPanel3.add(jLabel3, new org.netbeans.lib.awtextra.AbsoluteConstraints(280, 60, -1, -1));

    chroPanel.add(jPanel3, new org.netbeans.lib.awtextra.AbsoluteConstraints(210, 80, 400, 110));

    jPanel4.setBorder(javax.swing.BorderFactory.createTitledBorder("Ion Chromatogram m/z width options"));
    jPanel4.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    jCheckBox2.setSelected(true);
    jCheckBox2.setText("Automatic isotope window");
    jPanel4.add(jCheckBox2, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 30, -1, -1));

    jTextField5.setText("98.0");
    jPanel4.add(jTextField5, new org.netbeans.lib.awtextra.AbsoluteConstraints(180, 30, -1, -1));

    jLabel4.setText("ape");
    jPanel4.add(jLabel4, new org.netbeans.lib.awtextra.AbsoluteConstraints(220, 30, -1, -1));

    jCheckBox3.setText("Fixed m/z window");
    jCheckBox3.setEnabled(false);
    jPanel4.add(jCheckBox3, new org.netbeans.lib.awtextra.AbsoluteConstraints(20, 60, -1, -1));

    jTextField6.setEditable(false);
    jTextField6.setText("2.5");
    jPanel4.add(jTextField6, new org.netbeans.lib.awtextra.AbsoluteConstraints(180, 60, -1, -1));

    jLabel5.setText("+/- m/z");
    jPanel4.add(jLabel5, new org.netbeans.lib.awtextra.AbsoluteConstraints(220, 60, -1, -1));

    chroPanel.add(jPanel4, new org.netbeans.lib.awtextra.AbsoluteConstraints(210, 210, 310, 130));

    extractBtn.setText("Extract");
    extractBtn.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        extractBtnActionPerformed(evt);
      }
    });
    chroPanel.add(extractBtn, new org.netbeans.lib.awtextra.AbsoluteConstraints(310, 370, -1, -1));

    resetBtn.setText("Reset");
    resetBtn.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        resetBtnActionPerformed(evt);
      }
    });
    chroPanel.add(resetBtn, new org.netbeans.lib.awtextra.AbsoluteConstraints(390, 370, -1, -1));

    chorNoteLabel1.setText("Note. DTASelect-filter.txt and raw files are");
    chroPanel.add(chorNoteLabel1, new org.netbeans.lib.awtextra.AbsoluteConstraints(340, 30, -1, 20));

    chroNoteLabel2.setText("assumed to be in the same folder");
    chroPanel.add(chroNoteLabel2, new org.netbeans.lib.awtextra.AbsoluteConstraints(370, 50, -1, -1));

    peptideListBox.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptide List"));
    peptideListBox.setLayout(new java.awt.BorderLayout());

    quanPanel.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());

    fragIonScrollPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("FragIon"));

    fragIonPanel.setBackground(new java.awt.Color(255, 255, 255));
    fragIonPanel.setLayout(new org.netbeans.lib.awtextra.AbsoluteLayout());
    fragIonScrollPanel.setViewportView(fragIonPanel);

    quanPanel.add(fragIonScrollPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(1020, 10, 310, 510));

    setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
    setTitle("Census Main");

    jToolBar1.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
    jToolBar1.setToolTipText("open");

    open.setIcon(new javax.swing.ImageIcon(getClass().getResource("/open.gif"))); // NOI18N
    open.setToolTipText("new");
    open.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
    open.setMargin(new java.awt.Insets(2, 2, 2, 2));
    open.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        openActionPerformed(evt);
      }
    });
    jToolBar1.add(open);

    save.setIcon(new javax.swing.ImageIcon(getClass().getResource("/save.gif"))); // NOI18N
    save.setToolTipText("save");
    save.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
    save.setMargin(new java.awt.Insets(2, 2, 2, 2));
    save.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        saveActionPerformed(evt);
      }
    });
    jToolBar1.add(save);

    export.setIcon(new javax.swing.ImageIcon(getClass().getResource("/export.gif"))); // NOI18N
    export.setToolTipText("export report");
    export.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
    export.setMargin(new java.awt.Insets(2, 2, 2, 2));
    export.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        exportActionPerformed(evt);
      }
    });
    jToolBar1.add(export);

    filterBtn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/filter.gif"))); // NOI18N
    filterBtn.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
    filterBtn.setMargin(new java.awt.Insets(2, 2, 2, 2));
    filterBtn.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        filterBtnActionPerformed(evt);
      }
    });
    jToolBar1.add(filterBtn);

    filePathLabel.setToolTipText("File Path");
    jToolBar1.add(filePathLabel);

    proteinTable.addMouseListener(new java.awt.event.MouseAdapter() {
      public void mouseClicked(java.awt.event.MouseEvent evt) {
        proteinTableMouseClicked(evt);
      }
    });
    proteinListPanel.setViewportView(proteinTable);

    searchField.setMaximumSize(null);
    searchField.setPreferredSize(new java.awt.Dimension(39, 19));
    searchField.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        searchFieldActionPerformed(evt);
      }
    });

    searchBtn.setText("Search");
    searchBtn.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        searchBtnActionPerformed(evt);
      }
    });

    org.jdesktop.layout.GroupLayout proteinPanelLayout = new org.jdesktop.layout.GroupLayout(proteinPanel);
    proteinPanel.setLayout(proteinPanelLayout);
    proteinPanelLayout.setHorizontalGroup(
      proteinPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(proteinPanelLayout.createSequentialGroup()
          .add(searchField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 168, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(searchBtn)
          .addContainerGap(1174, Short.MAX_VALUE))
        .add(org.jdesktop.layout.GroupLayout.TRAILING, proteinListPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1399, Short.MAX_VALUE)
    );
    proteinPanelLayout.setVerticalGroup(
      proteinPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(proteinPanelLayout.createSequentialGroup()
          .add(proteinPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(searchBtn)
            .add(searchField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 23, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(proteinListPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 712, Short.MAX_VALUE)
          .add(111, 111, 111))
    );

    tabbedPanel.addTab("Proteins", proteinPanel);

    peptidePanel.setBackground(new java.awt.Color(102, 102, 255));

    peptideList.setBorder(javax.swing.BorderFactory.createTitledBorder("Peptides"));
    peptideList.setMinimumSize(new java.awt.Dimension(31, 200));
    peptideList.setPreferredSize(new java.awt.Dimension(31, 200));

    proteinSimplePane.setBorder(javax.swing.BorderFactory.createTitledBorder("Proteins"));

    proteinInfoPanel.setBorder(javax.swing.BorderFactory.createEmptyBorder(1, 1, 1, 1));

    proteinLabel.setText(" ");

    org.jdesktop.layout.GroupLayout proteinInfoPanelLayout = new org.jdesktop.layout.GroupLayout(proteinInfoPanel);
    proteinInfoPanel.setLayout(proteinInfoPanelLayout);
    proteinInfoPanelLayout.setHorizontalGroup(
      proteinInfoPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(proteinInfoPanelLayout.createSequentialGroup()
          .add(proteinLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1194, Short.MAX_VALUE)
          .addContainerGap())
    );
    proteinInfoPanelLayout.setVerticalGroup(
      proteinInfoPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(proteinLabel)
    );

    jSplitPane.setDividerLocation(700);
    jSplitPane.setMinimumSize(new java.awt.Dimension(100, 81));

    paramPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Scores"));
    paramPanel.setMinimumSize(new java.awt.Dimension(0, 0));
    paramPanel.setPreferredSize(new java.awt.Dimension(4, 19));

    rrLabel.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
    rrLabel.setText("R^2");

    shiftField.setColumns(3);

    shiftLabel.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
    shiftLabel.setText("shift");

    rrField.setColumns(10);
    rrField.setMaximumSize(new java.awt.Dimension(4, 19));

    measuredRatioLabel.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
    measuredRatioLabel.setText("Regression Ratio");

    regressionRatioField.setColumns(10);

    areaRatioField.setColumns(10);
    areaRatioField.setMaximumSize(new java.awt.Dimension(4, 19));

    measuredRatioLabel1.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
    measuredRatioLabel1.setText("Area Ratio");

    regLnLabel.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
    regLnLabel.setText("Ln (Regression Ratio)");

    areaRatioLogField.setColumns(10);

    pvalueLabel.setText("p value");

    jLabel1.setText("Data points");

    regScoreLabel.setText("Regression score (R)");

    org.jdesktop.layout.GroupLayout paramPanelLayout = new org.jdesktop.layout.GroupLayout(paramPanel);
    paramPanel.setLayout(paramPanelLayout);
    paramPanelLayout.setHorizontalGroup(
      paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(org.jdesktop.layout.GroupLayout.TRAILING, paramPanelLayout.createSequentialGroup()
          .addContainerGap()
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(rrLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 132, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(regLnLabel)
            .add(measuredRatioLabel1)
            .add(pvalueLabel)
            .add(jLabel1)
            .add(shiftLabel)
            .add(measuredRatioLabel)
            .add(regScoreLabel))
          .add(23, 23, 23)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
            .add(org.jdesktop.layout.GroupLayout.LEADING, regScoreField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 118, Short.MAX_VALUE)
            .add(org.jdesktop.layout.GroupLayout.LEADING, shiftField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 118, Short.MAX_VALUE)
            .add(pvalueField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 118, Short.MAX_VALUE)
            .add(org.jdesktop.layout.GroupLayout.LEADING, areaRatioField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .add(org.jdesktop.layout.GroupLayout.LEADING, paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
              .add(regressionRatioField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 19, Short.MAX_VALUE)
              .add(areaRatioLogField)
              .add(rrField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .add(org.jdesktop.layout.GroupLayout.LEADING, jTextField1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 118, Short.MAX_VALUE))
          .add(806, 806, 806))
    );
    paramPanelLayout.setVerticalGroup(
      paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(paramPanelLayout.createSequentialGroup()
          .add(2, 2, 2)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(regScoreField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(regScoreLabel))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(rrLabel)
            .add(rrField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(regLnLabel)
            .add(areaRatioLogField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(regressionRatioField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(measuredRatioLabel))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(pvalueField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(pvalueLabel))
          .add(10, 10, 10)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(measuredRatioLabel1)
            .add(areaRatioField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jTextField1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(jLabel1))
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(paramPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
            .add(shiftLabel)
            .add(shiftField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
          .addContainerGap(63, Short.MAX_VALUE))
    );

    jSplitPane.setRightComponent(paramPanel);

    pepTabbedPanel.setMaximumSize(new java.awt.Dimension(512, 302));
    pepTabbedPanel.setMinimumSize(new java.awt.Dimension(0, 0));

    chromatogramPanel.setLayout(new java.awt.BorderLayout());
    pepTabbedPanel.addTab("chromatogram", chromatogramPanel);

    correlationPanel.setBackground(new java.awt.Color(255, 255, 255));
    correlationPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Correlation"));
    pepTabbedPanel.addTab("correlation", correlationPanel);

    jSplitPane.setLeftComponent(pepTabbedPanel);

    org.jdesktop.layout.GroupLayout peptidePanelLayout = new org.jdesktop.layout.GroupLayout(peptidePanel);
    peptidePanel.setLayout(peptidePanelLayout);
    peptidePanelLayout.setHorizontalGroup(
      peptidePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(peptidePanelLayout.createSequentialGroup()
          .add(proteinSimplePane, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 173, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(peptidePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(peptideList, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1208, Short.MAX_VALUE)
            .add(proteinInfoPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, jSplitPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1208, Short.MAX_VALUE))
          .addContainerGap())
    );
    peptidePanelLayout.setVerticalGroup(
      peptidePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(peptidePanelLayout.createSequentialGroup()
          .add(proteinInfoPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(peptideList, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 436, Short.MAX_VALUE)
          .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
          .add(jSplitPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 318, Short.MAX_VALUE)
          .add(71, 71, 71))
        .add(proteinSimplePane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 853, Short.MAX_VALUE)
    );

    tabbedPanel.addTab("Peptides", peptidePanel);

    mainPanel.setLayout(new java.awt.BorderLayout());

    jLabel7.setText("\"\"");
    jLabel7.setVerticalAlignment(javax.swing.SwingConstants.TOP);
    jLabel7.setMaximumSize(new java.awt.Dimension(10, 30));
    jLabel7.setMinimumSize(new java.awt.Dimension(10, 30));
    jLabel7.setPreferredSize(new java.awt.Dimension(10, 50));
    mainPanel.add(jLabel7, java.awt.BorderLayout.NORTH);

    welcomeLabel.setFont(new java.awt.Font("Dialog", 1, 24));
    welcomeLabel.setForeground(new java.awt.Color(0, 0, 102));
    welcomeLabel.setText("Welcome to Census");
    welcomeLabel.setVerticalAlignment(javax.swing.SwingConstants.TOP);
    welcomeLabel.setAlignmentY(0.0F);
    mainPanel.add(welcomeLabel, java.awt.BorderLayout.CENTER);

    fileMenu.setLabel("File");
    fileMenu.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        fileMenuActionPerformed(evt);
      }
    });

    openItem.setText("Open chro file...");
    openItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        OpenChroFileActionPerformed(evt);
      }
    });
    fileMenu.add(openItem);

    reportItem.setText("Export Report...");
    reportItem.setEnabled(false);
    reportItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        reportItemActionPerformed(evt);
      }
    });
    fileMenu.add(reportItem);
    fileMenu.add(jSeparator1);

    exitItem.setLabel("Exit");
    exitItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        exitItemActionPerformed(evt);
      }
    });
    fileMenu.add(exitItem);

    jMenuBar1.add(fileMenu);

    runMenu.setText("Run");

    runItem.setText("run labeled data..");
    runItem.setActionCommand("run..");
    runItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        runItemActionPerformed(evt);
      }
    });
    runMenu.add(runItem);

    runNonLabel.setText("run label free data..");
    runNonLabel.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        runNonLabelActionPerformed(evt);
      }
    });
    runMenu.add(runNonLabel);

    confItem.setText("create config file..");
    confItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        confItemActionPerformed(evt);
      }
    });
    runMenu.add(confItem);

    jMenuBar1.add(runMenu);

    toolMenu.setText("Tools");
    toolMenu.setActionCommand("Menu");

    optionItem.setText("Options...");
    optionItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        optionItemActionPerformed(evt);
      }
    });
    toolMenu.add(optionItem);

    mergeItem.setText("Merge...");
    mergeItem.setToolTipText("merge CenSus output files");
    mergeItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        mergeItemActionPerformed(evt);
      }
    });
    toolMenu.add(mergeItem);

    alignSpectra.setText("align spectra...");
    alignSpectra.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        alignSpectraActionPerformed(evt);
      }
    });
    toolMenu.add(alignSpectra);

    openSpectra.setText("Open spectra file...");
    openSpectra.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        openSpectraActionPerformed(evt);
      }
    });
    toolMenu.add(openSpectra);

    mrmCsv.setText("MRM CSV generator");
    mrmCsv.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        mrmCsvActionPerformed(evt);
      }
    });
    toolMenu.add(mrmCsv);

    jMenuBar1.add(toolMenu);

    helpMenu.setText("Help");

    versionItem.setText("version");
    versionItem.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(java.awt.event.ActionEvent evt) {
        versionItemActionPerformed(evt);
      }
    });
    helpMenu.add(versionItem);

    jMenuBar1.add(helpMenu);

    setJMenuBar(jMenuBar1);

    org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
    getContentPane().setLayout(layout);
    layout.setHorizontalGroup(
      layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(jToolBar1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1417, Short.MAX_VALUE)
        .add(layout.createSequentialGroup()
          .add(tabbedPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1403, Short.MAX_VALUE)
          .addContainerGap(14, Short.MAX_VALUE))
        .add(mainPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1417, Short.MAX_VALUE)
    );
    layout.setVerticalGroup(
      layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
        .add(layout.createSequentialGroup()
          .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(jToolBar1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
            .add(layout.createSequentialGroup()
              .add(25, 25, 25)
              .add(tabbedPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 881, Short.MAX_VALUE)))
          .addContainerGap())
        .add(org.jdesktop.layout.GroupLayout.TRAILING, mainPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 918, Short.MAX_VALUE)
    );

    pack();
  }// </editor-fold>//GEN-END:initComponents

  private void versionItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_versionItemActionPerformed
// TODO add your handling code here:
    JOptionPane.showMessageDialog(this, "Census version " + conf.getVersion() + "\n\nCreated by \nSung Kyu (Robin) Park\nJohn Venable\n\nContact to rpark@scripps.edu for bugs or any questions", "Version", JOptionPane.INFORMATION_MESSAGE);
  }//GEN-LAST:event_versionItemActionPerformed

  private void alignSpectraActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_alignSpectraActionPerformed
// TODO add your handling code here:

    JFileChooser choose = new JFileChooser();
    choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Select config file");
    choose.setApproveButtonText("Select config file");

    if (this.currentDirectory != null) {
      choose.setCurrentDirectory(new File(this.getCurrentDirectory()));
    }

    int returnVal = choose.showOpenDialog(this);

    if (returnVal == choose.CANCEL_OPTION) {
      return;
    }

    final File file = choose.getSelectedFile();
    final String workFolder = file.getParent();

    final ChroProgressDialog chroProgress = new ChroProgressDialog(this, false);
    chroProgress.setLocationRelativeTo(this);
    chroProgress.setResizable(false);
    chroProgress.setVisible(true);

    final RelExMainFrame mFrame = this;
    try {

      SAXBuilder sb = new SAXBuilder();
      Document doc = sb.build(file);
      final Element root = doc.getRootElement();

      Thread t = new Thread() {

        private boolean isSuccessful = false;

        public void run() {

          try {

            Configuration conf = Configuration.getInstance();
            String parentPath = file.getParent();

            conf.readXMLParam(parentPath + File.separator + file.getName());
            conf.setFilePath(parentPath);

            List sampleList = root.getChildren("sample");

            //Configuration.
            String refFileName = "";
            Element refEle = root.getChild("ref");
            if (null != refEle) {
              //String refSamName = refEle.getChildText("sample_name");
              refFileName = refEle.getChildText("file_name");

              conf.setRefFileName(refFileName);
            }

            int index = 0;
            int refIndex = 0;

            Vector<String> fileNameList = new Vector<String>();
            Vector<String> sampleNameList = new Vector<String>();
            Vector<String> pathFileNameList = new Vector<String>();

            //confSam.setRefFileName(refFileName);
            //Hashtable<String, Hashtable> masterHt = new Hashtable<String, Hashtable>();
            HashSet set = new HashSet();

            for (Iterator<Element> itr = sampleList.iterator(); itr.hasNext(); ) {
              Element sam = itr.next();
              String samName = sam.getChildText("name");

              //populate configuration class
              Configuration.Sample confSam = new Configuration.Sample();
              confSam.setName(samName);

              List fileList = sam.getChild("ms_files").getChildren("file");

              for (Iterator<Element> itr1 = fileList.iterator(); itr1.hasNext(); ) {

                Element eachFile = itr1.next();
                String fileName = eachFile.getText();
                //      masterHt.put(fileName, new Hashtable());

                confSam.addFile(fileName);
                pathFileNameList.add(fileName);

                set.add(fileName.substring(0, fileName.lastIndexOf("/")));

                if (fileName.endsWith("ms2")) {
                  fileName = fileName.substring(0, fileName.length() - 3);
                  fileName += "ms1";

                }

                fileNameList.add(fileName);
                sampleNameList.add(samName);

                if (null != refEle && refFileName.equals(fileName)) {
                  refIndex = index;
                }

                index++;
              }

              conf.addExp(confSam);
            }

            long start = System.currentTimeMillis();
            chroalign align = new chroalign();

            int[][][] pathArray = null;

            boolean isAlign = conf.isAlign();

            chroProgress.setProgress(1);

            File alignFile = new File(workFolder + File.separator + "aligned_out.xml");
            File pathoutFile = new File(workFolder + File.separator + "path_out.xml");
            File chrooutFile = new File(workFolder + File.separator + "chro_out.xml");

            int selection = -1;
            if (alignFile.exists() && pathoutFile.exists() && chrooutFile.exists()) {
              Object[] options = {"Yes", "No"};

              selection = JOptionPane.showOptionDialog(mFrame,
                "Alignment out files are found.  Do you want to read them?",
                "Alignment out files found",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                options,
                options[0]
              );
            }

            if (selection != 0) //align
            {
              if (isAlign) //align based on chromatogram profile
              {
                pathArray = align.alignChro(chroProgress, fileNameList.toArray(), sampleNameList.toArray(), refIndex, false, 500, file.getParent());

                NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
                conf.setMapModel(mapModel);
              } else //align based on retention time
              {
                pathArray = align.noAlignChro(chroProgress, fileNameList.toArray(), sampleNameList.toArray(), file.getParent());
                NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
                mapModel.reinitializeMaxIndexByRet(); //this method is called only for the aligning by retention time.  Not nice idea.
                conf.setMapModel(mapModel);
              }

              System.out.println(System.currentTimeMillis());
              String[] targetMS1Files = align.getTargetMS1Files();
              String referenceMS1File = align.getReferenceMS1File();

              conf.setNonlabelFilePaths(set);

              chroProgress.addMessage("Reading Spectra and Running Quantification...");

            }

            this.isSuccessful = true;

          } catch (Exception e) {
            System.out.println("error " + e);
            e.printStackTrace();
            isSuccessful = false;
          }

          SwingUtilities.invokeLater(new Runnable() {
            public void run() {
              //progress.setVisible(false);
              chroProgress.setVisible(false);
              chroProgress.hide();

              if (isSuccessful) {
                SimplePlot sPlot = new SimplePlot();
                //        AlignPlot sPlot = new AlignPlot();
                //      mFrame.getProteinPanel().add(sPlot);
                sPlot.setVisible(true);
                sPlot.Go(workFolder);

              } else {
                JOptionPane.showMessageDialog(mFrame, "Failed to align spectra", "Failed to align spectra", JOptionPane.ERROR_MESSAGE);
              }

            }
          });
        }
      };

      t.start();

    } catch (JDOMException ex) {

      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + ex, JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + ex);
      System.out.println("Error " + ex);
    } catch (IOException ex) {

      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + ex.toString(), JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + ex);
      System.out.println("Error " + ex);
    } catch (Exception e) {
      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + e.toString(), JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + e);
      System.out.println("Error " + e);
    }


  }//GEN-LAST:event_alignSpectraActionPerformed

  public static void mrmCsvTest(java.awt.event.ActionEvent evt) {

    File file = new File("/home/rpark/001_emily_data");
    if (null == file) {
      return;
    }

    //spectrumField.setText(choose.getSelectedFile().getAbsolutePath()); //.getCurrentDirectory().toString());
    File[] fList = file.listFiles();

    boolean isDtaSelect = false;
    boolean isMs2File = false;
    boolean isMs1File = false;

    for (int i = 0; i < fList.length; i++) {
      if ("DTASelect-filter.txt".equals(fList[i].getName())) {
        isDtaSelect = true;
      } else if (fList[i].getName().endsWith("ms1")) {
        isMs1File = true;
      } else if (fList[i].getName().endsWith("ms2")) {
        isMs2File = true;
      }
    }

    StringBuffer errorSb = new StringBuffer();
    if (!isDtaSelect) {
      errorSb.append("DTASlect-filter.txt is not found.\n");
    }
    if (!isMs1File) {
      errorSb.append("ms1 files are not found.\n");
    }
    if (!isMs2File) {
      errorSb.append("ms2 files are not found.\n");
    }

    if (errorSb.length() > 0) {
      return;
    }

//        final ChroProgressDialog chroProgress = new ChroProgressDialog(this, false);
    //      chroProgress.setLocationRelativeTo(this);
    //    chroProgress.setResizable(false);
    //  chroProgress.setVisible(true);
    Configuration conf = Configuration.getInstance();
    conf.setFilePath(file.getAbsolutePath());

    Thread t = new Thread() {
      boolean isSuccessful = false;

      public void run() {
        try {
          ChroGenerator chro = new ChroGenerator(null);

          chro.createMRMFragmentIons(null);
          isSuccessful = true;
        } catch (IOException e) {
          e.printStackTrace();
          isSuccessful = false;

        } catch (Exception e) {
          e.printStackTrace();
          isSuccessful = false;
        }

        SwingUtilities.invokeLater(new Runnable() {
          public void run() {
            if (isSuccessful) {
            }
          }
        });
      }
    };

    try {
      t.start();

    } catch (Exception e) {
      t = null;
    }
  }


  private void mrmCsvActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mrmCsvActionPerformed
// TODO add your handling code here:
    JFileChooser choose = new JFileChooser();
    choose.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    //choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Select Spectrum Directory");

    if (this.currentDirectory != null) {
      choose.setCurrentDirectory(new File(this.currentDirectory));
    }

    int returnVal = choose.showOpenDialog(this);

    File file = choose.getSelectedFile();
    if (null == file) {
      return;
    }

    //spectrumField.setText(choose.getSelectedFile().getAbsolutePath()); //.getCurrentDirectory().toString());
    File[] fList = file.listFiles();

    boolean isDtaSelect = false;
    boolean isMs2File = false;
    boolean isMs1File = false;

    for (int i = 0; i < fList.length; i++) {
      if ("DTASelect-filter.txt".equals(fList[i].getName())) {
        isDtaSelect = true;
      } //else if(fList[i].getName().endsWith("ms1"))
      //    isMs1File = true;
      else if (fList[i].getName().endsWith("ms2")) {
        isMs2File = true;
      }
    }

    StringBuffer errorSb = new StringBuffer();
    if (!isDtaSelect) {
      errorSb.append("DTASlect-filter.txt is not found.\n");
    }
    //if(!isMs1File)
    //    errorSb.append("ms1 files are not found.\n");
    if (!isMs2File) {
      errorSb.append("ms2 files are not found.\n");
    }

    if (errorSb.length() > 0) {
      JOptionPane.showMessageDialog(this, errorSb.toString(), "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);
      return;
    }

    final ChroProgressDialog chroProgress = new ChroProgressDialog(this, false);
    chroProgress.setLocationRelativeTo(this);
    chroProgress.setResizable(false);
    chroProgress.setVisible(true);

    conf.setFilePath(file.getAbsolutePath());

    //final String spectrumFolder = this.spectrumField.getText();
    //final String configFile = this.configField.getText();
    //final String dtaFile = this.idFileField.getText();
    Thread t = new Thread() {
      boolean isSuccessful = false;

      public void run() {
        try {

          ChroGenerator chro = new ChroGenerator(chroProgress);

          chroProgress.addMessageWithLine("Census starts to analyze data...");
          chro.createMRMFragmentIons(chroProgress);
                    /*
                     switch(conf.getQuantLevel())
                     {
                     case 1: //Full
                     chro.createFullscanXmlChro();
                     break;

                     case 2: //MSMS
                     chro.createMsmsXmlChro(chroProgress);
                     break;

                     default :
                     break;
                     }
                     */
          //if( isFull.equals("f") ) //Full Scan
          //    chro.createFullscanXmlChro();
          //else //msms scan
          //    chro.createMsmsXmlChro();

          isSuccessful = true;
        } catch (IOException e) {
          e.printStackTrace();
          isSuccessful = false;

        } catch (Exception e) {
          e.printStackTrace();
          isSuccessful = false;
          //JOptionPane.showMessageDialog(this, "Failed to generate a chro file: " + e, "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);
        }

        SwingUtilities.invokeLater(new Runnable() {
          public void run() {
            chroProgress.setVisible(false);
            chroProgress.hide();

            if (isSuccessful) {

              openMRMCrvFile(conf.getFilePath() + "mrm_frags.xml");

            }
          }
        });
      }
    };

    try {
      t.start();

      if (null != conf.getErrorMessage()) {
        throw new Exception();
      }

    } catch (Exception e) {
      t = null;
      JOptionPane.showMessageDialog(this, "Failed to generate a chro file: " + e, "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);
      //chroProgress.setVisible(false);
    }

  }//GEN-LAST:event_mrmCsvActionPerformed

  private void nonlabelTableKeyPressed(java.awt.event.KeyEvent evt) {

    //JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file", JOptionPane.ERROR_MESSAGE);
  }

  private void nonlabelTableMouseClicked(java.awt.event.MouseEvent evt) {

    //System.out.println("single clicked");
    JTable table = (JTable) evt.getSource();

    generateNonLabelData(currentPeptide, table.getSelectedRow());

    //JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file", JOptionPane.ERROR_MESSAGE);
  }

  private void peptideListTableKeyPressed(java.awt.event.KeyEvent evt) {
// TODO add your handling code here:

    if (evt.getKeyCode() != KeyEvent.VK_DOWN && evt.getKeyCode() != KeyEvent.VK_UP) {
      return;
    }

    int pepSize = currentProtein.getPeptideList().size();

    JTable table = (JTable) evt.getSource();

    //ChroPeptide peptide=null;
    if (evt.getKeyCode() == KeyEvent.VK_DOWN) {
      if (table.getSelectedRow() + 1 >= pepSize) {
        return;
      }

      this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(table.getSelectedRow() + 1);
    } else if (evt.getKeyCode() == KeyEvent.VK_UP) {
      if (table.getSelectedRow() <= 0) {
        return;
      }

      this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(table.getSelectedRow() - 1);
    }

    selectQuantType();
  }

  private void selectQuantType() {
    PostOptions options = PostOptions.getInstance();

    //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.
    if (cr.getExpType() > 0) {
      switch (cr.getExpType()) {
        case CensusConstants.MSMS_DATA_INDEPENDENT:

          if (options.isFilterFragmentIons()) {
            generateInDepFragData(this.currentPeptide);
          } else {
            generateInDepData(this.currentPeptide);
          }

          break;
        case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
        case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
          //generate itraq data
          generateITRAQData(this.currentPeptide);
          break;

        default:
          break;
      }

    } else if (!this.isLabeled()) {
      //if( this.quantLevel==1 )
      generateNonLabelData(currentPeptide);
      //else if(this.quantLevel==2 )
      //  generateNonLabelData(currentPeptide);
    } //else if(this.isDataDependent)
    else if (this.quantLevel == 1) {
      generateDepData(this.currentPeptide);
    } else if (options.isFilterFragmentIons()) {
      generateInDepFragData(this.currentPeptide);
    } else {
      generateInDepData(this.currentPeptide);
    }
  }

  private void peptideListTableMouseClicked(java.awt.event.MouseEvent evt) {
// TODO add your handling code here:

    if (this.isChroFile) {
      JTable table = (JTable) evt.getSource();

      ChroPeptide peptide = (ChroPeptide) currentProtein.getPeptideList().get(table.getSelectedRow());
      this.currentPeptide = peptide;

      selectQuantType();
    }
  }

  private void searchFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_searchFieldActionPerformed
// TODO add your handling code here:
  }//GEN-LAST:event_searchFieldActionPerformed

  public void runtest() {
    this.runNonLabelActionPerformed(null);
  }

  //private void runNonLabelActionPerformed(java.awt.event.ActionEvent evt) {
  //label free
  public void runNonLabelActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_runNonLabelActionPerformed
// TODO add your handling code here:
    //comment out for testing.

    JFileChooser choose = new JFileChooser();
    choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Select config file");
    choose.setApproveButtonText("Select");

    if (this.currentDirectory != null) {
      choose.setCurrentDirectory(new File(this.getCurrentDirectory()));
    }

    int returnVal = choose.showOpenDialog(this);

    if (returnVal == choose.CANCEL_OPTION) {
      return;
    }

    final File file = choose.getSelectedFile();
    final String workFolder = file.getParent();

    final ChroProgressDialog chroProgress = new ChroProgressDialog(this, false);
    chroProgress.setLocationRelativeTo(this);
    chroProgress.setResizable(false);
    chroProgress.setVisible(true);

    // final RelExMainFrame mFrame = this;
//        for testing
//        final File file = new File("/home/rpark/001/nonlabel/IsotopeFree/census_config.xml");
    //      end of testing
    try {

      SAXBuilder sb = new SAXBuilder();
      Document doc = sb.build(file);
      final Element root = doc.getRootElement();

//commented out for testing
      final JFrame mFrame = this;
      Thread t = new Thread() {
        private boolean isSuccessful = true;

        public void run() {

          try {
            //reads conf. file

            Configuration conf = Configuration.getInstance();
            //conf.setTextArea(chroProgress.getProgressText());
            String parentPath = file.getParent();

            conf.readXMLParam(parentPath + File.separator + file.getName());
            conf.setFilePath(parentPath);

//                        System.out.println(conf.isHighRes() + " " + conf.getMassTolerance());
            List sampleList = root.getChildren("sample");

            //Configuration.
            String refFileName = "";
            Element refEle = root.getChild("ref");
            if (null != refEle) {
              //String refSamName = refEle.getChildText("sample_name");
              refFileName = refEle.getChildText("file_name");

              conf.setRefFileName(refFileName);
            }

            int index = 0;
            int refIndex = 0;

            Vector<String> fileNameList = new Vector<String>();
            Vector<String> sampleNameList = new Vector<String>();
            Vector<String> pathFileNameList = new Vector<String>();

            //confSam.setRefFileName(refFileName);
            //Hashtable<String, Hashtable> masterHt = new Hashtable<String, Hashtable>();
            HashSet set = new HashSet();

            for (Iterator<Element> itr = sampleList.iterator(); itr.hasNext(); ) {
              Element sam = itr.next();
              String samName = sam.getChildText("name");

              //populate configuration class
              Configuration.Sample confSam = new Configuration.Sample();
              confSam.setName(samName);

              List fileList = sam.getChild("ms_files").getChildren("file");

              for (Iterator<Element> itr1 = fileList.iterator(); itr1.hasNext(); ) {

                Element eachFile = itr1.next();
                String fileName = eachFile.getText();
                //      masterHt.put(fileName, new Hashtable());

                confSam.addFile(fileName);
                pathFileNameList.add(fileName);

                set.add(fileName.substring(0, fileName.lastIndexOf(File.separator)));

                if (fileName.endsWith("ms2")) {
                  fileName = fileName.substring(0, fileName.length() - 3);
                  fileName += "ms1";

                }

                File tmpFile = new File(fileName);
                if (!tmpFile.exists()) {
                  isSuccessful = false;
                  throw new FileNotFoundException(fileName + " file not found ");
                }

                fileNameList.add(fileName);
                sampleNameList.add(samName);

                if (null != refEle && refFileName.equals(fileName)) {
                  refIndex = index;
                }

                index++;
              }

              conf.addExp(confSam);
            }

            conf.setNonlabelFilenameList(pathFileNameList);

                        /*
                         // move later to the chrogenerator.java
                         List tmpFileNameList = conf.getNonlabelFilenameList();
                         Hashtable<String, int[]> spHt = new Hashtable<String, int[]>();

                         int fileCount=0;
                         for(Iterator<String> itr=tmpFileNameList.iterator(); itr.hasNext(); )
                         {
                         String eachPath = itr.next();
                         eachPath = eachPath.substring(0, eachPath.lastIndexOf(File.separator) + 1);
                         //eachPath += "DTASelect-filter.txt";
                         eachPath += CensusConstants.SEARCH_OUTPUT;

                         DTASelectFilterReader dtaReader = new DTASelectFilterReader(eachPath);

                         for (Iterator<Protein> itr1 = dtaReader.getProteins(); itr1.hasNext(); )
                         {
                         Protein protein = itr1.next();
                         String accession = protein.getLocus();
                         int[] tmpArr = spHt.get(accession);

                         if(null == tmpArr)
                         {
                         tmpArr = new int[tmpFileNameList.size()];
                         tmpArr[fileCount] = Integer.parseInt(protein.getSpectrumCount());

                         spHt.put(accession, tmpArr);
                         }
                         else
                         {
                         tmpArr[fileCount] = Integer.parseInt(protein.getSpectrumCount());
                         }

                         }

                         fileCount++;
                         }

                         conf.setSpHt(spHt);
                         */
            long start = System.currentTimeMillis();
            chroalign align = new chroalign();

            int[][][] pathArray = null;

            boolean isAlign = conf.isAlign();

            chroProgress.setProgress(1);

            File alignFile = new File(workFolder + File.separator + "aligned_out.xml");
            File pathoutFile = new File(workFolder + File.separator + "path_out.xml");
            File chrooutFile = new File(workFolder + File.separator + "chro_out.xml");

            int selection = -1;
            if (alignFile.exists() && pathoutFile.exists() && chrooutFile.exists()) {
              Object[] options = {"Yes", "No"};

              selection = JOptionPane.showOptionDialog(mFrame,
                "Alignment out files are found.  Do you want to read them?",
                "Alignment out files found",
                JOptionPane.YES_NO_OPTION,
                JOptionPane.QUESTION_MESSAGE,
                null,
                options,
                options[0]
              );
            }

            if (selection != 0) //align
            {
              chroProgress.addMessageWithLine("Start aligning spectra...");
              if (isAlign) //align based on chromatogram profile
              {
                //pathArray = align.alignChro(chroProgress.getProgressBar(), fileNameList.toArray(), sampleNameList.toArray(), refIndex, false, 500, file.getParent());
                pathArray = align.alignChro(chroProgress, fileNameList.toArray(), sampleNameList.toArray(), refIndex, false, 500, file.getParent());

                chroProgress.addMessageWithLine("Populating mapping model...");
                NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
                conf.setMapModel(mapModel);
              } else //align based on retention time
              {
                pathArray = align.noAlignChro(chroProgress, fileNameList.toArray(), sampleNameList.toArray(), file.getParent());
                chroProgress.addMessageWithLine("Populating mapping model...");
                NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
                mapModel.reinitializeMaxIndexByRet(); //this method is called only for the aligning by retention time.  Not nice idea.
                conf.setMapModel(mapModel);
              }
              chroProgress.addMessageWithLine("Spectra alignment done.");
            } else //already aligned
            {
              chroProgress.addMessageWithLine("Reading pre-aligned xml file...");

              SAXBuilder builder = new SAXBuilder();
              Document doc = builder.build(pathoutFile);
              Element rootEle = doc.getRootElement();

              int numberOfAlignment = 1;  //number of alignment is one bigger than actual alignment number :-(
              int biggestNum = 0;

              for (Iterator<Element> itr = rootEle.getChildren("dataset").iterator(); itr.hasNext(); ) {
                Element ele = itr.next();

                if (biggestNum < ele.getChildren().size()) {
                  biggestNum = ele.getChildren().size();
                }
                numberOfAlignment++;
              }

              pathArray = new int[numberOfAlignment][2][biggestNum + 1];

              int firstIndex = 0;

              for (Iterator<Element> itr = rootEle.getChildren("dataset").iterator(); itr.hasNext(); ) {
                Element ele = itr.next();

                int thirdIndex = 0;
                for (Iterator<Element> pitr = ele.getChildren("p").iterator(); pitr.hasNext(); ) {
                  Element pele = pitr.next();

                  pathArray[firstIndex][0][thirdIndex] = Integer.parseInt(pele.getAttributeValue("x"));
                  pathArray[firstIndex][1][thirdIndex] = Integer.parseInt(pele.getAttributeValue("y"));
//System.out.println( pathArray[firstIndex][0][thirdIndex]  + "\t" + pathArray[firstIndex][1][thirdIndex] );
                  thirdIndex++;
                }

                firstIndex++;
              }

              chroProgress.addMessageWithLine("Populating mapping model...");

              NonLabelMappingModel mapModel = new NonLabelMappingModel(pathArray, pathFileNameList, refIndex);
              conf.setMapModel(mapModel);

              chroProgress.addMessageWithLine("Pre-alignment file loading done.");

            }

            String[] targetMS1Files = align.getTargetMS1Files();
            String referenceMS1File = align.getReferenceMS1File();

            conf.setNonlabelFilePaths(set);

            chroProgress.addMessageWithLine("Reading Spectra and Running Quantification...");

            ChroGenerator chro = new ChroGenerator(
              //aJProgressBar, //commented out for testing
              chroProgress, //.getProgressBar(),
              //null,
              null,
              referenceMS1File,
              targetMS1Files,
              root,
              pathArray
              //isotopeFileField.getText().trim(),
            );

            //if( 1 == conf.getQuantLevel() )
            chro.createNonlabelXmlChro();
            //else if( 2 == conf.getQuantLevel() )
            //    chro.createNonlabelMsmsXmlChro();

          } catch (FileNotFoundException fe) {
            JOptionPane.showMessageDialog(chroProgress, "Error : " + fe.getMessage(), "Error ", JOptionPane.ERROR_MESSAGE);
            fe.printStackTrace();
            isSuccessful = false;

          } catch (Exception e) {
            JOptionPane.showMessageDialog(chroProgress, "Error : " + e.getMessage(), "Error ", JOptionPane.ERROR_MESSAGE);
            e.printStackTrace();

            isSuccessful = false;
          }

          SwingUtilities.invokeLater(new Runnable() {
            public void run() {
              chroProgress.invalidate();
              chroProgress.repaint();

              chroProgress.setVisible(false);
              //chroProgress.hide();

              if (!isSuccessful) {
                return;
              }

              openChroFile(file.getParent() + File.separator + "census_chro.xml");
            }
          });
        }
      };

      t.start();

    } catch (JDOMException ex) {

      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + ex, JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + ex);
      System.out.println("Error " + ex);
    } catch (IOException ex) {

      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + ex.toString(), JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + ex);
      System.out.println("Error " + ex);
    } catch (Exception e) {
      JOptionPane.showMessageDialog(this, "Error reading the config file.", "Error reading the config file" + e.toString(), JOptionPane.ERROR_MESSAGE);
      chroProgress.addMessageWithLine("Error reading the config file" + e);
      System.out.println("Error " + e);
    }

  }//GEN-LAST:event_runNonLabelActionPerformed

  private void openSpectraActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openSpectraActionPerformed
// TODO add your handling code here:

    //Configuration conf = Configuration.getInstance();
    //conf.setSimpleIndexGenerator(true);
    JFileChooser choose = new JFileChooser();
    choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Select Spectra File");
    //choose.addChoosableFileFilter( new SimpleFileNameFilter("xml", "Chro File") );

    if (currentDirectory != null && !"".equals(currentDirectory)) {
      choose.setCurrentDirectory(new File(currentDirectory));
    }

    int returnVal = choose.showOpenDialog(chroPanel);
    this.specFile = choose.getSelectedFile();

    if (null == specFile || returnVal == choose.CANCEL_OPTION) {
      return;
    }

    currentDirectory = specFile.getAbsolutePath();
    currentDirectory = currentDirectory.substring(0, currentDirectory.lastIndexOf(File.separator));

    filePathLabel.setText(specFile.getAbsolutePath());

    cleanupProteinTableModel();

    this.openSpecFile(specFile.getAbsolutePath());

  }//GEN-LAST:event_openSpectraActionPerformed

  private void searchBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_searchBtnActionPerformed
// TODO add your handling code here:

    int rows = this.proteinTable.getRowCount();
    int columns = this.proteinTable.getColumnCount();

    Object value = null;
    String searchText = this.searchField.getText().trim();

    this.searchIndex = (this.searchIndex >= rows) ? 0 : this.searchIndex;

    for (int i = this.searchIndex; i < rows; i++) {
      value = this.proteinTableModel.getValueAt(i, 0);

      if (null != value && value.toString().contains(searchText)) {
        this.proteinTable.setRowSelectionInterval(i, i);
        this.searchIndex = (this.searchIndex >= rows) ? 0 : i + 1;
        this.proteinTable.scrollRectToVisible(this.proteinTable.getCellRect(i, 0, true));

        return;
      }

      value = this.proteinTableModel.getValueAt(i, columns - 1);

      if (null != value && value.toString().contains(searchText)) {
        this.proteinTable.setRowSelectionInterval(i, i);
        this.searchIndex = (this.searchIndex >= rows) ? 0 : i + 1;
        this.proteinTable.scrollRectToVisible(this.proteinTable.getCellRect(i, 0, true));
        //    this.proteinTable.setColumnSelectionInterval(0,0);
        return;
      }
    }

    if (this.searchIndex != 0) {
      for (int i = 0; i < this.searchIndex; i++) {
        value = this.proteinTableModel.getValueAt(i, 0);

        if (null != value && value.toString().contains(searchText)) {
          this.proteinTable.setRowSelectionInterval(i, i);
          this.searchIndex = (this.searchIndex >= rows) ? 0 : i + 1;
          this.proteinTable.scrollRectToVisible(this.proteinTable.getCellRect(i, 0, true));
          //    this.proteinTable.setColumnSelectionInterval(0,0);
          return;
        }

        value = this.proteinTableModel.getValueAt(i, columns - 1);

        if (null != value && value.toString().contains(searchText)) {
          this.proteinTable.setRowSelectionInterval(i, i);
          this.searchIndex = (this.searchIndex >= rows) ? 0 : i + 1;
          this.proteinTable.scrollRectToVisible(this.proteinTable.getCellRect(i, 0, true));
          //    this.proteinTable.setColumnSelectionInterval(0,0);
          return;
        }
      }
    }

    JOptionPane.showMessageDialog(this, "No protein was found.", "No protein was found.", JOptionPane.ERROR_MESSAGE);
  }//GEN-LAST:event_searchBtnActionPerformed

  public void openSpecFile(final String specFile) {
    this.specFile = new File(specFile);
    this.isChroFile = false;

    try {
      this.welcomeLabel.setVisible(false);
      this.tabbedPanel.setVisible(true);

      final JDialog progress = new JDialog(this);
      JProgressBar aJProgressBar = new JProgressBar(0, 100);
      aJProgressBar.setIndeterminate(true);
      Container cp = progress.getContentPane();
      JLabel jb = new JLabel("Opening DTASelect-filter.out file...",
        SwingConstants.CENTER);
      cp.add(jb, BorderLayout.SOUTH);
      cp.add(aJProgressBar, BorderLayout.NORTH);
      progress.setSize(500, 100);
      progress.setLocationRelativeTo(this);
      progress.pack();

      progress.setResizable(false);
      progress.setVisible(true);

      this.tabbedPanel.remove(this.peptidePanel);

      tabbedPanel.addTab("Peptides", qualPanel);

      final String tempDir = this.currentDirectory;

      final JFrame mFrame = this;
      Thread t = new Thread() {
        public void run() {

          try {

            dtaReader = new DTASelectFilterReader(specFile);
            File f = new File(specFile);
            Configuration conf = Configuration.getInstance();
            conf.setSimpleIndexGenerator(true);

            ht = ChroGenerator.createIndexedFiles(f.getParent() + File.separator, "ms2");

            //qualPanel = new QualificationPanel(mFrame, tempDir, ht);
            qualPanel.setMFrame(mFrame);
            qualPanel.setCurrentDirectory(tempDir);
            qualPanel.setHt(ht);

            proteinList = dtaReader.getChroProteinList();

            Object[] simpleProteinArr = new Object[1];

            for (Iterator<ChroProtein> itr = proteinList.iterator(); itr.hasNext(); ) {
              ChroProtein protein = itr.next();
              proteinTableModel.addRow(protein.getProteinData());

              simpleProteinArr[0] = protein.getLocus();

              proteinSimpleTableModel.addRow(simpleProteinArr);
            }

            tabbedPanel.setSelectedIndex(0);
            tabbedPanel.setEnabledAt(0, true);

            //irisPanel = new IrisPanel(ht, this.currentPeptide);
            //this.peptidePanel.invalidate();
            //this.peptidePanel.validate();
            //this.peptidePanel.repaint();
            //} catch (IOException e)
          } catch (Exception e) {
            System.out.println("error " + e);
            e.printStackTrace();
          }

          SwingUtilities.invokeLater(new Runnable() {
            public void run() {
              progress.setVisible(false);
            }
          });
        }
      };

      t.start();

    } /*
         catch(IOException ioe)
         {
         JOptionPane.showMessageDialog(proteinListPanel, "Failed to open a chro file: " + ioe, "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
         ioe.printStackTrace();
         }
         catch(JDOMException je)
         {
         JOptionPane.showMessageDialog(proteinListPanel, "Failed to open a chro file: " + je, "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
         je.printStackTrace();
         } */ catch (Exception je) {
      JOptionPane.showMessageDialog(proteinListPanel, "Failed to open a chro file: " + je, "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
      je.printStackTrace();
    }

  }

  public void openMRMCrvFile(String chroFile) {

    String filePath = conf.getFilePath();
    File mrmChroFile = new File(filePath);
    if (!mrmChroFile.exists()) {
      JOptionPane.showMessageDialog(this, "MRM Chro file was not generated", "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);
      return;
    }

    try {
      this.welcomeLabel.setVisible(false);

      MRMFragPanel mrmPanel = new MRMFragPanel(chroFile);

      //pack();
      //this.setVisible((true);
      //this.mainPanel.add(mrmPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 50, 1417, 868));
      //mainPanel.add(welcomeLabel, java.awt.BorderLayout.CENTER);
      mainPanel.add(mrmPanel, BorderLayout.CENTER);

      //mainPanel.add(welcomeLabel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 50, 1417, 868));
      //this.getContentPane().add(mrmPanel);
            /*
             org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
             getContentPane().setLayout(layout);
             layout.setHorizontalGroup(
             layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
             .add(jToolBar1, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1417, Short.MAX_VALUE)
             .add(layout.createSequentialGroup()
             .addContainerGap()
             .add(welcomeLabel)
             .addContainerGap(905, Short.MAX_VALUE))
             .add(layout.createSequentialGroup()
             .add(mrmPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 1403, Short.MAX_VALUE)
             .addContainerGap(14, Short.MAX_VALUE))
             );
             layout.setVerticalGroup(
             layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
             .add(layout.createSequentialGroup()
             .add(layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
             .add(layout.createSequentialGroup()
             .add(jToolBar1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
             .add(30, 30, 30)
             .add(welcomeLabel))
             .add(layout.createSequentialGroup()
             .add(25, 25, 25)
             .add(mrmPanel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 831, Short.MAX_VALUE)))
             .addContainerGap())
             );
             */
    } catch (Exception e) {
      e.printStackTrace();
      JOptionPane.showMessageDialog(this, "Failed to open MRM Chro file" + e, "Failed to open MRM Chro file", JOptionPane.ERROR_MESSAGE);

    }
  }

  //    private JTextArea taskOutput;
  private JProgressBar progressBar;

  /**
   * Invoked when task's progress property changes.
   */
  public void propertyChange(java.beans.PropertyChangeEvent evt) {

    if ("progress" == evt.getPropertyName()) {
      int progress = (Integer) evt.getNewValue();
      progressBar.setIndeterminate(false);

      progressBar.setValue(progress);
      //taskOutput.append(String.format(
      //		"Completed %d%% of task.\n", progress));
    }
  }

  public void postProcessOpenChroFile(ChroXmlReader cr, ArrayList<ChroProtein> proteinList, JDialog progressDialog) {
    this.cr = cr;
    this.proteinList = proteinList;

    isDataDependent = cr.isDataDependent();
    quantLevel = cr.getQuantLevel();
    labeled = cr.isLabeled();

    chroFileOpen = true;
    reportItem.setEnabled(true);

    //System.out.println(this.isLabeled());
    if (!isLabeled()) {
      if (cr.getQuantLevel() == 1) {
        fragIonScrollPanel.setVisible(false);
      } else {
        fragIonScrollPanel.setVisible(true);
      }

      correlationPanel.setVisible(false);

      int tabCount = pepTabbedPanel.getTabCount();
      for (int i = 1; i < tabCount; i++) //remove tab except chromatograms
      {
        pepTabbedPanel.removeTabAt(1);
      }

      //this.paramPanel.setVisible(false);
      paramPanel.removeAll();

      nonLabelPanel.setBackground(new java.awt.Color(255, 255, 255));
      nonLabelPanel.setLayout(new BorderLayout());
      nonLabelPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "All", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Dialog", 0, 11), java.awt.Color.black));
      //quanPanel.add(nonLabelPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(710, 190, 580, 300));
      nonlabelTable.setModel(nonlabelTableModel);

      //add action to the table
      nonlabelTable.addKeyListener(new java.awt.event.KeyAdapter() {
        public void keyPressed(java.awt.event.KeyEvent evt) {
          nonlabelTableKeyPressed(evt);
        }
      });

      nonlabelTable.addMouseListener(new java.awt.event.MouseAdapter() {
        public void mouseClicked(java.awt.event.MouseEvent evt) {
          nonlabelTableMouseClicked(evt);
        }
      });

      nonlabelScrollPane.setViewportView(nonlabelTable);

      nonLabelPanel.add(nonlabelScrollPane);
      nonLabelPanel.setMinimumSize(new Dimension(10, 10));
      nonLabelPanel.setMaximumSize(new Dimension(30, 30));

      nonLabelSummaryPanel.setBackground(new java.awt.Color(255, 255, 255));
      nonLabelSummaryPanel.setLayout(new BorderLayout());
      nonLabelSummaryPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Summary", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Dialog", 0, 11), java.awt.Color.black));
      //quanPanel.add(nonLabelSummaryPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(710, 490, 580, 300));

      nonlabelSummaryTable.setModel(nonlabelSummaryTableModel);
      nonlabelSummaryScrollPane.setViewportView(nonlabelSummaryTable);
      nonLabelSummaryPanel.add(nonlabelSummaryScrollPane);
      nonLabelSummaryPanel.setMinimumSize(new Dimension(10, 10));
      nonLabelSummaryPanel.setMaximumSize(new Dimension(30, 30));

      paramPanel.setLayout(new BorderLayout());

      paramPanel.add(nonLabelPanel, BorderLayout.CENTER);

    }

    //iTRAQ
    if (cr.getExpType() == CensusConstants.MSMS_SPECIFIC_SINGLE_MASS || cr.getExpType() == CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS) {
      paramPanel.setVisible(false);
    }

    welcomeLabel.setVisible(false);
    tabbedPanel.setVisible(true);

    plot = new ChromatographPlot(this);
    pepDistPlot = new PeptideDistPlot(this);

    plot.setBackground(new Color(255, 255, 255));
    pepDistPlot.setBackground(new Color(255, 255, 255));

    chromatogramPanel.add(plot);
    proteinRatioDistPanel.add(pepDistPlot);

    chroXmlParser = new ChroXMLParser((Plot) plot);
    peptideDistParser = new PlotMLParser((Plot) pepDistPlot);

    Object[] simpleProteinArr = new Object[2];
    for (Iterator<ChroProtein> itr = proteinList.iterator(); itr.hasNext(); ) {
      ChroProtein protein = itr.next();
      proteinTableModel.addRow(protein.getProteinData());

      simpleProteinArr[0] = protein.getLocus();
      simpleProteinArr[1] = protein.getDescription();

      proteinSimpleTableModel.addRow(simpleProteinArr);

    }

    progressDialog.setVisible(false);
  }

  public void preProcessChroOpen() {

  }

  public void openChroFile(java.net.URI chroFileName) {
    fileUri = chroFileName;
    this.chroFile = null;
    openChroFileImpl();
  }

  public void openChroFile(String chroFileName) {
    this.chroFile = new File(chroFileName);
    openChroFileImpl();
  }

  public void openChroFileImpl() {
    this.isChroFile = true;

    progressBar = new JProgressBar(0, 100);
    progressBar.setValue(0);

    //Call setStringPainted now so that the progress bar height
    //stays the same whether or not the string is shown.
    progressBar.setStringPainted(true);

    final JDialog progressDialog = new JDialog(this);
    Container cp1 = progressDialog.getContentPane();
    cp1.add(progressBar, BorderLayout.NORTH);
    progressDialog.setSize(500, 100);
    progressDialog.setLocationRelativeTo(this);
    progressDialog.pack();

    progressDialog.setResizable(false);
    progressDialog.setVisible(true);

    ProgressTask pTask = new ProgressTask(this, chroFile, progressDialog);
    pTask.addPropertyChangeListener(this);
    pTask.execute();

    this.reportItem.setEnabled(true);

    try {
      SAXBuilder sb = new SAXBuilder();

    } catch (Exception ex) {
      JOptionPane.showMessageDialog(this, "Error reading the census chro file.", "Error reading the chro file" + ex, JOptionPane.ERROR_MESSAGE);
      System.out.println("Error " + ex);
    }

  }

  private void runItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_runItemActionPerformed
// TODO add your handling code here:
    OpenConfigDialog dialog = new OpenConfigDialog(this, true);
    dialog.pack();
    dialog.setLocationRelativeTo(this);
    dialog.setVisible(true);
    dialog.setResizable(false);

  }//GEN-LAST:event_runItemActionPerformed

  private void confItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_confItemActionPerformed
    // TODO add your handling code here:
    SelectConfDialog dialog = new SelectConfDialog(this, true);
    dialog.pack();
    dialog.setLocationRelativeTo(this);
    dialog.setVisible(true);
    dialog.setResizable(false);

  }//GEN-LAST:event_confItemActionPerformed

  private void mergeItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mergeItemActionPerformed
    // TODO add your handling code here:

    MergeDialog dialog = new MergeDialog(this, true);
    dialog.pack();
    dialog.setLocationRelativeTo(this);
    dialog.setVisible(true);
    dialog.setResizable(false);

  }//GEN-LAST:event_mergeItemActionPerformed

  public void mergeFiles(
    final Object[] arr,
    final File outputFile,
    final String dbFilePath,
    final FilterModel fModel,
    final double correctFactorValue
  ) {
    try {
      //final String[] arr1 = {"/home/rpark/a1.txt", "/home/rpark/a2.txt", "/home/rpark/a3.txt",  };

      //display progress bar
      final JDialog progress = new JDialog(this);
      //final ChroPeptide tempPeptide = this.currentPeptide;

      final JProgressBar aJProgressBar = new JProgressBar(0, 100);
      aJProgressBar.setStringPainted(true);
      Container cp = progress.getContentPane();
      JLabel jb = new JLabel("Merging files...", SwingConstants.CENTER);
      cp.add(jb, BorderLayout.SOUTH);
      cp.add(aJProgressBar, BorderLayout.NORTH);
      progress.setSize(500, 100);
      progress.setLocationRelativeTo(this);
      progress.pack();
      progress.setResizable(false);
      progress.setVisible(true);

      Thread t = new Thread() {
        private boolean isSuccessful = true;
        private String errorMessage = "";
        private PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

        public void run() {

          printHeader(p);

          p.println("H\tMerged Data");

          for (int i = 0; i < arr.length; i++) {
            p.print("H\tMerged File : ");
            p.println(arr[i]);
          }

          p.print("H\tCorrection Factor Value : ");
          p.println(correctFactorValue);

          if (fModel.isDetSelect()) {
            p.print("H\tDeterminant Factor : ");
            p.println(fModel.getDetValue());
          } else {
            p.println("H\tNo Determinant Factor");
          }

          if (fModel.isPValueSelect()) {
            p.print("H\tOutlier pValue : ");
            p.println(fModel.getPValue());
          } else {
            p.println("H\tNo Outlier pValue");
          }

          if (fModel.isFilterFragmentIons()) {
            p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
          }

          p.println("H\tPLINE\tLOCUS\tWEIGHTED_AVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
          p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME");

          StringBuffer result = new StringBuffer();

          int totalCount = 0;
          int quantifiedCount = 0;
          int redunProteinCount = 0;
          int uniqueProteinCount = 0;

          CenSusReportReader cReader = new CenSusReportReader();

          try {

            List<MergeProteinModel> proteinList = cReader.merge(arr, aJProgressBar, dbFilePath, fModel, correctFactorValue);

            for (Iterator<MergeProteinModel> itr = proteinList.iterator(); itr.hasNext(); ) {
              MergeProteinModel protein = itr.next();
              if (protein.getPeptides().size() <= 0) {
                continue;
              }

              if (fModel.isPValueSelect()) {
                protein.removeOutlier(fModel.getPValue());
              }

              for (Iterator<ChroProtein> pItr = protein.getProteins().iterator(); pItr.hasNext(); ) {
                ChroProtein cPro = pItr.next();

                result.append("P\t");
                result.append(cPro.getLocus());
                result.append("\t");

                if (Double.compare(protein.getAverageRatio(), Double.NaN) == 0) {
                  result.append("NA");
                } else {
                  result.append(CensusHelper.format.format(protein.getAverageRatio()));
                }
                result.append("\t");

                if (Double.compare(protein.getStdev(), Double.NaN) == 0) {
                  result.append("NA");
                } else {
                  result.append(CensusHelper.format.format(protein.getStdev()));
                }

                result.append("\t");
                result.append(protein.getPeptides().size());

                result.append("\t");
                //result.append( cPro.getSpectrumCount() );
                result.append(protein.getTotalSpecCount());

                result.append("\t");
                result.append(cPro.getDescription());
                result.append("\n");
                redunProteinCount++;

              }

              uniqueProteinCount++;

              Vector pepList = new Vector();
              pepList.addAll(protein.getPeptides());

              for (Iterator<MergeProteinModel.Peptide> pepItr = pepList.iterator(); pepItr.hasNext(); ) {
                MergeProteinModel.Peptide pep = pepItr.next();

                result.append("S\t");
                result.append(pep.isUnique() ? "U" : "");
                result.append("\t");
                result.append(pep.getSequence());
                result.append("\t");
                result.append(CensusHelper.format.format(pep.getRatio()));
                result.append("\t");
                result.append(CensusHelper.format.format(pep.getRegFactor()));
                result.append("\t");
                result.append(CensusHelper.format.format(pep.getRegFactor() * pep.getRegFactor()));
                result.append("\t");
                result.append((pep.getSamIntensity() < 0) ? "" : pep.getSamIntensity());
                result.append("\t");
                result.append((pep.getRefIntensity() < 0) ? "" : pep.getRefIntensity());
                result.append("\t");
                result.append(pep.getAreaRatio());
                result.append("\t");
                result.append(pep.getProfileScore());
                result.append("\t");
                result.append(pep.getFileName());
                result.append("\n");

                quantifiedCount++;
              }
            }
          } catch (IOException e) {
            isSuccessful = false;
            errorMessage = e.getMessage();

            e.printStackTrace();
          } catch (Exception e) {
            isSuccessful = false;
            errorMessage = e.getMessage();

            e.printStackTrace();
          }

          final int finalTotalCount = cReader.getTotalPeptideCount();
          final int finalQuantifiedCount = quantifiedCount; //quantifiedCount;

          SwingUtilities.invokeLater(new Runnable() {
                                       public void run() {
                                         progress.setVisible(false);
                                         progress.hide();

                                         if (isSuccessful) {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             "Merging files completed.",
                                             "Merging files completed.",
                                             JOptionPane.PLAIN_MESSAGE);
                                         } else {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             errorMessage,
                                             "Report file Creation",
                                             JOptionPane.ERROR_MESSAGE);

                                           outputFile.delete();
                                         }

                                       }

                                     }
          );

          p.print("H\t");
          p.print("Total Redundant Proteins\t");
          p.println(redunProteinCount);
          p.print("H\t");
          p.print("Total Unique Proteins\t");
          p.println(uniqueProteinCount);
          p.print("H\t");
          p.print("Total peptides\t");
          p.println(finalTotalCount);
          p.print("H\t");
          p.print("Quantified peptides\t");
          p.print(finalQuantifiedCount);
          p.print("\n");
          p.print("H\t");
          p.print("Quantification efficiency\t");
          p.print(CensusHelper.format.format((double) finalQuantifiedCount / finalTotalCount * 100));
          p.print(" %\n");
          p.print(result.toString());

          if (null != p) {
            p.close();
          }

        }
      };

      t.start();
    } catch (IOException ioe) {
      JOptionPane.showMessageDialog(proteinListPanel, "Failed to open a chro file: " + ioe, "Failed to open a chro file", JOptionPane.ERROR_MESSAGE);
      ioe.printStackTrace();
    }

  }

  private void exportActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exportActionPerformed
// TODO add your handling code here:
    reportItemActionPerformed(evt);
  }//GEN-LAST:event_exportActionPerformed

  private void saveActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveActionPerformed
// TODO add your handling code here:

    try {
      //ChroGenerator.createFullscanXmlChro(proteinList, isDataDependent);

    } catch (Exception e) {
      System.out.println("Error : " + e.toString());
    }


  }//GEN-LAST:event_saveActionPerformed

  public void toggleFilterIcon(boolean isFilter) {

    if (isFilter) {
      filterBtn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/filter.gif")));
    } else {
      filterBtn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/nofilter.gif")));
    }
  }

  private void filterBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_filterBtnActionPerformed
// TODO add your handling code here:

    PostOptions options = PostOptions.getInstance();

    if (options.isFilterFragmentIons()) {
      options.setFilterFragmentIons(false);
      toggleFilterIcon(false);
      //filterBtn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/nofilter.gif")));
    } else {
      options.setFilterFragmentIons(true);
      toggleFilterIcon(true);
      //filterBtn.setIcon(new javax.swing.ImageIcon(getClass().getResource("/filter.gif")));
    }

    if (!this.isLabeled() && this.quantLevel == 1) {
      generateNonLabelData(currentPeptide);
    } else if (this.isDataDependent) {
      generateDepData(currentPeptide);
    } else if (options.isFilterFragmentIons()) {
      generateInDepFragData(currentPeptide);
    } else {
      generateInDepData(currentPeptide);
    }

    if (cr.getQuantLevel() == 1) {
      updatePeptideInfo();
    } else if (cr.getQuantLevel() == 2) {
      if (options.isFilterFragmentIons()) {
        updateMS2PeptideFilterInfo(null, -1);
      } else {
        updatePeptideInfo();
      }
    }

  }//GEN-LAST:event_filterBtnActionPerformed

  private void openActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openActionPerformed
// TODO add your handling code here:
    this.OpenChroFileActionPerformed(evt);
  }//GEN-LAST:event_openActionPerformed

  private void optionItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_optionItemActionPerformed

    //options.setFilterFragmentIons(this.filterFragmentIons);
    //options.setDisplayFragmentIons(this.isDisplayFragmentIons());
    OptionsDialog dialog = new OptionsDialog(this, true);
    dialog.pack();
    dialog.setLocationRelativeTo(this);
    dialog.setVisible(true);
    dialog.setResizable(false);

    // TODO add your handling code here:
  }//GEN-LAST:event_optionItemActionPerformed

  private void exportITRAQSingleReport(final File file, final double intensityThreshold) {
    exportITRAQSingleReport(file, intensityThreshold, true);

  }

  public static String exportITRAQSingleReportCoreRef(
    PrintStream p,
    double intensityThreshold,
    JProgressBar aJProgressBar,
    String errorMessage,
    boolean isSuccessful,
    ArrayList<ChroProtein> proteinList,
    String path) {

    Configuration conf = Configuration.getInstance();
    //List massMonitorList = conf.getMsmsMassArr();
    List<ReportIon> massMonitorList = conf.getReportIonList();

    StringBuffer result = new StringBuffer();

    printHeader(p);
    p.println("H\tCensus msms analysis");
    p.println("H\tIntensity Threshold\t" + intensityThreshold);
    p.print("H\tPLINE\tLOCUS\tSPEC_COUNT\tPEP_NUM\t");

    for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
      ReportIon each = itr.next();
      p.print("average m/z_");
      p.print(each);
      p.print("\t");
      p.print("norm_average m/z_");
      p.print(each);
      p.print("\t");
    }

    String massStr = massMonitorList.get(0).toString();
    //			if(massMonitorList.size()>1)
    //for(int i=0;i<massMonitorList.size();i++)
    for (Iterator<ReportIon> ditr = massMonitorList.iterator(); ditr.hasNext(); ) {
      p.print("norm_ratio(" + ditr.next() + ")");
      p.print("\t");
    }

    p.println("DESCRIPTION");

    p.print("H\tSLINE\tUNIQUE\tSEQUENCE\t"); //RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tSIGNAL_TO_NOISE_RATIO\tFILE_NAME");

    for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
      ReportIon each = itr.next();

      p.print("m/z_");
      p.print(each);
      p.print("_int");
      p.print("\t");
      p.print("norm_m/z_");
      p.print(each);
      p.print("_int");
      p.print("\t");
    }

    for (Iterator<ReportIon> ditr = massMonitorList.iterator(); ditr.hasNext(); ) {
      p.print("norm_ratio(" + ditr.next() + ")");
      p.print("\t");
    }

    p.print("SpC\tScanNum\tCState\tFilename");
    p.println();

    //                        int totalCount=0;
    int quantifiedCount = 0;

    try {
      double eachSeg = (double) 100 / proteinList.size();
      double percent = 0;
      boolean isSameGroup = false;

      TDoubleArrayList[] normIntList = new TDoubleArrayList[massMonitorList.size()];

      for (int ni = 0; ni < normIntList.length; ni++) {
        normIntList[ni] = new TDoubleArrayList();
      }

      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();

        List<ChroPeptide> peptideList = protein.getPeptideList();
        //List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();
        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
          ChroPeptide peptide = pepItr.next();

          ChroData cdata = peptide.getData(0);
          peptide.setTotalIntArr(cdata.getIntensityArr());
        }
      }

      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();

        int proLength = -1;
        if (!"N/A".equals(protein.getLength())) {
          proLength = Integer.parseInt(protein.getLength());
        }

        String[] specArr = protein.getSpectrumCount().split(",");

        isSameGroup = true;
        StringBuffer proteinSb = new StringBuffer();

        List<ChroPeptide> peptideList = protein.getPeptideList();

        proteinSb.append("P\t");
        proteinSb.append(protein.getLocus()).append("\t");
        proteinSb.append(protein.getSpectrumCount()).append("\t");  //
        proteinSb.append(peptideList.size()).append("\t"); //bug??

        double devSum = 0;

        StringBuffer pepSb = new StringBuffer();

        int tmpCount = 0;

        long[] intSumArr = new long[massMonitorList.size()];
        long[] normIntSumArr = new long[massMonitorList.size()];
        double[] ratioSumArr = new double[massMonitorList.size()];

        for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
          ChroPeptide each = tempItr.next();
          long[] larr = each.getTotalIntArr();
          double avgPepInt = 0;

          long tmpSum = 0;
          for (long l : larr) {
            tmpSum += l;
          }

          pepSb.append("S\t");
          pepSb.append(each.isUnique() ? "U" : "");
          pepSb.append("\t");
          pepSb.append(each.getSequence());
          pepSb.append("\t");

          long[] normArr = new long[larr.length];

          pepSb.append(larr[0]).append("\t");
          //System.out.println(larr[0]);

          for (int i = 1; i < larr.length; i++) {
            pepSb.append(larr[i]).append("\t");
            intSumArr[i] += larr[i];
            //avgPepInt += larr[i];
            //double logint = Math.log(larr[i]) + normCorrectionList.get(i);
            double logint = Math.log(larr[i]);
            double revInt = Math.exp(logint);
            normIntSumArr[i] += revInt;
            normArr[i] = (long) revInt;
            avgPepInt += revInt;

//System.out.println("===" + larr[i] + " " + revInt);
            pepSb.append((int) revInt).append("\t");
//			normIntList[i].add( logint );

          }

          avgPepInt = avgPepInt / normArr.length;
//System.out.println("===" + avgPepInt);

          for (int i = 0; i < normArr.length; i++) {
            if (avgPepInt > 0) {
              ratioSumArr[i] += (double) normArr[i] / avgPepInt;
            }
          }

          for (int i = 0; i < normArr.length; i++) {
            if (avgPepInt <= 0) {
              pepSb.append("NA\t");
            } else {
              pepSb.append(CensusHelper.format.format((double) normArr[i] / avgPepInt)).append("\t");

//System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
            }
          }

          pepSb.append(each.getSpecCount()).append("\t");
          pepSb.append(each.getScanNum()).append("\t");
          pepSb.append(each.getChargeState()).append("\t");
          pepSb.append(each.getFileName()).append("\n");

          tmpCount++;
        }

        if (tmpCount <= 0) {
          continue;
        }

        double ratioTotalSum = 0;

//		for(long l : intSumArr) {
        for (int i = 0; i < intSumArr.length; i++) {
          proteinSb.append(CensusHelper.format.format((double) intSumArr[i] / tmpCount)).append("\t")
            .append(CensusHelper.format.format((double) normIntSumArr[i] / tmpCount)).append("\t");
        }

        for (double d : ratioSumArr) {
          ratioTotalSum += d;
        }

        if (ratioSumArr.length > 1) {
          for (int i = 0; i < ratioSumArr.length; i++) {

            if (ratioTotalSum <= 0) {
              proteinSb.append("NA\t");
            } else {
              proteinSb.append(CensusHelper.format.format(ratioSumArr[i] / tmpCount)).append("\t");
            }

          }
        }

        proteinSb.append(protein.getDescription());
        proteinSb.append("\n");

        percent += eachSeg;
        if (null != aJProgressBar) {
          aJProgressBar.setValue((int) percent);
        }

        if (pepSb.length() <= 0) {
          continue;
        }

        result.append(proteinSb.toString());
        result.append(pepSb.toString());

      }

            /*
             while( null != (eachLine = br.readLine()) ) {
             if(!eachLine.startsWith("S"))
             continue;

             String[] arr = eachLine.split("\t");
             list.add( Math.log( Double.parseDouble(arr[3])) );
             }


             */
    } catch (Exception e) {
      isSuccessful = false;
      errorMessage = e.getMessage();
      e.printStackTrace();
    }

    return result.toString();
  }

  public static Map<String,String> extractTMTPhosphoLocalizationScore(String path) throws IOException {
    BufferedReader br = new BufferedReader(new FileReader(path));
    Map<String,String> localScoreMap = new HashMap<>();
    String line;
    boolean start= false;
    int localLoc = -1;
    int sequenceLoc = -1;
    while((line=br.readLine())!=null)
    {
      if(!start)
      {
        if(line.startsWith("Unique"))
        {
          start = true;
          String [] arr = line.split("\t");
          for(int i=0; i<arr.length; i++)
          {
            if(arr[i].startsWith("Localization Score"))
            {
              localLoc =i;
            }
            else if(arr[i].startsWith("Sequence"))
            {
              sequenceLoc = i;
            }
          }
          if(localLoc<0) return null;
        }
      }
      else if(start && (line.startsWith("*")||line.startsWith("\t")))
      {
        String [] arr = line.split("\t");
        String[] prekey = arr[1].split("\\.");
        String filename = prekey[0];
        String cs = prekey[prekey.length-1];
        String sequence = arr[sequenceLoc];
        String key = sequence+filename+cs;
        String values = arr[localLoc];
        localScoreMap.put(key,values);
      }
    }

    return localScoreMap;
  }

  //itraq TMT tmt iTRAQ
  public static void exportITRAQSingleReportCore(
    PrintStream p,
    double intensityThreshold,
    int intThresholdType,
    JProgressBar aJProgressBar,
    String errorMessage,
    boolean isSuccessful,
    ArrayList<ChroProtein> proteinList,
    String path,
    String type, ReportParam rParam) {

    Configuration conf = Configuration.getInstance();
    List<ReportIon> massMonitorList = conf.getReportIonList();

   // StringBuffer result = new StringBuffer();
    int Protein_length = 0;
    double totalIntensity = 0.0;
    printHeader(p);
    p.println("H\tCensus msms analysis");

    String phosphoPath =  path+"/DTASelect-filter.txt.phospho";
    File phospho = new File(phosphoPath);
    Map<String,String> localMap = new HashMap<>();
    if(phospho.exists())
    {
      try {
        localMap = extractTMTPhosphoLocalizationScore(phosphoPath);
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

    if (1 == intThresholdType) {
      //List<Double> intList = new ArrayList<Double>();
      DescriptiveStatistics stats = new DescriptiveStatistics();
      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();
        Protein_length = Integer.parseInt(protein.getLength());
        List<ChroPeptide> peptideList = protein.getPeptideList();
        //List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
          ChroPeptide peptide = pepItr.next();

          ChroData cdata = peptide.getData(0);
          peptide.setTotalIntArr(cdata.getIntensityArr());
        }

        for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
          ChroPeptide each = tempItr.next();
          long[] larr = each.getTotalIntArr();

          long tmpSum = 0;

          for (long l : larr) {
            tmpSum += l;
          }
          //intList.add( Math.log(tmpSum)/Math.log(10) );

          if (tmpSum > 0) {
            stats.addValue(Math.log(tmpSum) / Math.log(10));
          }
        }

      }

      p.println("H\tIntensity ThresholdType\tz-score\t" + intensityThreshold);
      intensityThreshold = stats.getMean() - intensityThreshold * stats.getStandardDeviation();
      if (intensityThreshold < 0) {
        intensityThreshold = 0;
      }

      intensityThreshold = Math.pow(10, intensityThreshold);
      p.println("H\tIntensity Threshold\t" + intensityThreshold);

    } else if (2 == intThresholdType) {
      p.println("H\tIntensity Threshold Type\tfixed_value");
      p.println("H\tIntensity Threshold\t" + intensityThreshold);
      //p.println("H\tAverage intensity in protein line uses only non-zero intensities from peptides\t" + intensityThreshold);
    }

    p.println("H\tPrecursor purity threshold\t" + rParam.getPurityThreshold());
    p.println("H\tPrecursor SN filter\t" + rParam.getSignalToNoiseThreshold());


    p.print("H\tPLINE\tLOCUS\tSPEC_COUNT\tSEQ_COUNT\tSEQ_COVERAGE\tLENGTH\tMOLWT\tpI\tPEP_NUM\t");

    for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
      ReportIon each = itr.next();
      p.print("total m/z_");
      p.print(each.getMass());
      p.print("\t");
      p.print("average m/z_");
      p.print(each.getMass());
      p.print("\t");
      p.print("norm_total m/z_");
      p.print(each.getMass());
      p.print("\t");
      p.print("norm_average m/z_");
      p.print(each.getMass());
      p.print("\t");
      p.print("norm_median m/z_");
      p.print(each.getMass());
      p.print("\t");
    }

//	String massStr = massMonitorList.get(0).toString();
    //			if(massMonitorList.size()>1)
    //for(int i=0;i<massMonitorList.size();i++)
	/*
         for(Iterator<ReportIon> ditr=massMonitorList.iterator(); ditr.hasNext(); ) {
         p.print("norm_ratio(" + ditr.next().getMass() + ")");
         p.print("\t");
         } */
    p.print("TOTAL_PEPTIDE_INTENSITY_DIVIDEDBY_PROTEIN_LENGTH_LOG10");
    p.print("\t");
    p.print("LENGTH");
    p.print("\t");
    p.println("DESCRIPTION");

    p.print("H\tSLINE\tUNIQUE\tSEQUENCE\t"); //RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tSIGNAL_TO_NOISE_RATIO\tFILE_NAME");

    for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
      ReportIon each = itr.next();

      p.print("m/z_");
      p.print(each.getMass());
      p.print("_int");
      p.print("\t");
      p.print("norm_m/z_");
      p.print(each.getMass());
      p.print("_int");
      p.print("\t");
    }


        /*
         for(Iterator<ReportIon> ditr=massMonitorList.iterator(); ditr.hasNext(); ) {
         p.print("norm_ratio(" + ditr.next().getMass() + ")");
         p.print("\t");
         } */
    p.print("Localization_Score\t");
    p.print("TMT_purity\t");
    p.print("Signal-noise\t");
    // p.print("Report_Ion_Signal-noise\t");
    p.print("SpC\tIon Count\tScanNum\tCState\tFilename");
    p.print("\tPTMIndex\tPTMIndexProteins");

    p.println();

    //                        int totalCount=0;
    int quantifiedCount = 0;
    String outlierLevel = Configuration.getInstance().getTmtOutlierLevel();

//System.out.println("outlier level : " +  outlierLevel);
    try {
      double eachSeg = (double) 100 / proteinList.size();
      double percent = 0;
      boolean isSameGroup = false;

      List intensityList = new ArrayList();
      TDoubleArrayList[] normIntList = new TDoubleArrayList[massMonitorList.size()];
      TLongArrayList[] rawIntList = new TLongArrayList[massMonitorList.size()];

      for (int ni = 0; ni < normIntList.length; ni++) {
        normIntList[ni] = new TDoubleArrayList();
        rawIntList[ni] = new TLongArrayList();
      }

      //run for loop to calculate dist
      Map<String, Integer> seqCountMap = new HashMap<>();
      Map<String, Integer> ionCountMap = new HashMap<>();

      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();

        List<ChroPeptide> peptideList = protein.getPeptideList();
        //List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();
        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
          ChroPeptide peptide = pepItr.next();
          String seq = peptide.getSequence();
          String ionSeq = seq.concat(Integer.toString(peptide.getChargeState()));
          Integer count = seqCountMap.get(seq);
          if (count != null) {
            seqCountMap.put(seq, count + 1);
          } else {
            seqCountMap.put(seq, 1);
          }
          count = ionCountMap.get(ionSeq);
          if (count != null) {
            ionCountMap.put(ionSeq, count + 1);
          } else {
            ionCountMap.put(ionSeq, 1);
          }

          ChroData cdata = peptide.getData(0);
          peptide.setTotalIntArr(cdata.getIntensityArr());
        }

        for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
          ChroPeptide each = tempItr.next();
          long[] larr = each.getTotalIntArr();

          long tmpSum = 0;

          for (long l : larr) {
            tmpSum += l;
          }


          if (intensityThreshold < 0) {

            long baseline = getNoiseBaseline(each.getScanNum(), path + File.separator + each.getFileName() + ".ms2");

            //System.out.println("------------" + each.getScanNum() + " " + each.getFileName() + " " + baseline);
            boolean biggerThanBaseline = false;
            for (long l : larr) {
              if (l > baseline) {
                biggerThanBaseline = true;
              }
            }

            //if(biggerThanBaseline)
            //	System.out.println("------------keep");
            //else
            //	System.out.println("------------no");
          }

          if (tmpSum < intensityThreshold) {
            continue;
          }

          for (int i = 0; i < larr.length; i++) {
            if (larr[i] <= 0) {
              continue;
            }

            //double logint = Math.log(larr[i])/Math.log(2);
            double logint = Math.log(larr[i]);
            //System.out.println(i+ "\t" + larr[i] + "\t" + logint + "\t" );

            rawIntList[i].add(larr[i]);
            normIntList[i].add(logint);
          }

        }
      }

      List<Double> correctionList = new ArrayList<Double>();
      double correctionAvg = 0;

      if ("mode".equals(type)) {
        for (TDoubleArrayList norm : normIntList) {
          edu.scripps.pms.util.stats.Histogram hist = new edu.scripps.pms.util.stats.Histogram(120, 0, 20);
          hist.setData(norm.toNativeArray());

          double[] bins = hist.getBins();
          int[] freqArr = hist.getFreqArrSmooth(7);

          double maxBin = 0;
          double maxValue = 0;

          for (int i = 0; i < freqArr.length - 1; i++) {

            if (maxValue < freqArr[i]) {
              maxValue = freqArr[i];
              maxBin = bins[i];
            }

            //		if(count==5)
            //	System.out.println("++==\t" + bins[i] +"\t" + freqArr[i]);
          }

          int underFlow = hist.getUnderFlows() + freqArr[0];
          int overFlow = freqArr[freqArr.length - 1] + hist.getOverFlows();
          //System.out.println(">>>==\t" + underFlow + " " + overFlow);
          correctionAvg += maxBin;
          correctionList.add(maxBin);
        }

        correctionAvg /= correctionList.size();

        List<Double> normCorrectionList = new ArrayList<Double>();
        for (int i = 0; i < correctionList.size(); i++) {

          normCorrectionList.add(correctionAvg - correctionList.get(i));
        }

        //System.out.println("=-==" + correctionAvg);
        //System.out.println("=-==" + normCorrectionList);
        for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
          ChroProtein protein = proItr.next();

          int proLength = -1;
          if (!"N/A".equals(protein.getLength())) {
            proLength = Integer.parseInt(protein.getLength());
          }
          StringBuilder proteinCalcSb = new StringBuilder();
          String[] specArr = protein.getSpectrumCount().split(",");

          isSameGroup = true;
          StringBuffer proteinSb = new StringBuffer();

          List<ChroPeptide> peptideList = protein.getPeptideList();
          int numUniquePeptides =0;
          int numPeptides =0;
          proteinSb.append("P\t");
          proteinSb.append(protein.getLocus()).append("\t");
          proteinSb.append(protein.getSpectrumCount()).append("\t");  //
          proteinSb.append(peptideList.size()).append("\t"); //bug??

          double devSum = 0;

          StringBuffer pepSb = new StringBuffer();

          int tmpCount = 0;

          long[] intSumArr = new long[massMonitorList.size()];
          long[] normIntSumArr = new long[massMonitorList.size()];
          int[] normIntCountArr = new int[massMonitorList.size()];
          int[] intCountArr = new int[massMonitorList.size()];
          double[] ratioSumArr = new double[massMonitorList.size()];
          Set<String>  sequenceSet = new HashSet<>();

          for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
            ChroPeptide each = tempItr.next();

            long[] larr = each.getTotalIntArr();
            double avgPepInt = 0;

            long tmpSum = 0;
            for (long l : larr) {
              tmpSum += l;
            }

            if (tmpSum < intensityThreshold) {
              continue;
            }

            if (each.getTmtPurity() < rParam.getPurityThreshold()) {
              continue;
            }
            if (rParam.getSignalToNoiseThreshold() > 0 && each.getSignalNoise() < rParam.getSignalToNoiseThreshold()) {
              continue;
            }
            if(rParam.isIsUniquePeptide() && !each.isUnique())continue;
            sequenceSet.add(each.getSequence());
            pepSb.append("S\t");
            pepSb.append(each.isUnique() ? "U" : "");
            pepSb.append("\t");
            pepSb.append(each.getSequence());
            pepSb.append("\t");

            long[] normArr = new long[larr.length];
            //if(protein.getLocus().equals("IPI00396680.3"))
            for (int i = 0; i < larr.length; i++) {
              pepSb.append(larr[i]).append("\t");

              //if(protein.getLocus().equals("IPI00396680.3"))
              //	System.out.println(i + "-----\t" + larr[i]);
              //avgPepInt += larr[i];
              double logint = Math.log(larr[i]) + normCorrectionList.get(i);

              double revInt = Math.exp(logint);

              if (larr[i] > 0) {
                //if(protein.getLocus().equals("IPI00396680.3"))
                //	System.out.println("\t" + i + "--\t" + Math.log(larr[i]) + " " + revInt);
                intSumArr[i] += larr[i];
                intCountArr[i]++;
                normIntSumArr[i] += revInt;
                avgPepInt += revInt;
                normArr[i] = (long) revInt;
                pepSb.append((int) revInt).append("\t");

                normIntCountArr[i]++;

              } else {
                pepSb.append("0").append("\t");
              }

              //if(i==1)
              //System.out.println("==22\t" + normIntSumArr[i]);
              //System.out.println("===" + larr[i] + " " + revInt);
              //			normIntList[i].add( logint );
            }

            avgPepInt = avgPepInt / normArr.length;
            //System.out.println("===" + avgPepInt);

            for (int i = 0; i < normArr.length; i++) {
              if (avgPepInt > 0) {
                ratioSumArr[i] += (double) normArr[i] / avgPepInt;
              }
            }

                        /*
                         for(int i=0;i<normArr.length;i++)
                         if(avgPepInt<=0)
                         pepSb.append("0\t");
                         else {
                         pepSb.append( CensusHelper.format.format((double)normArr[i]/avgPepInt) ).append("\t");

                         //System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
                         } */
            if(each.isUnique())numUniquePeptides++;
            numPeptides++;
            String seq = each.getSequence();
            String ionSeq = seq.concat(Integer.toString(each.getChargeState()));
            int specCount = seqCountMap.get(seq);
            int ionCount = ionCountMap.get(ionSeq);
            double purity = each.getTmtPurity();
            double signal = each.getSignalNoise();
            pepSb.append(purity).append("\t");
            pepSb.append(signal).append("\t");
            pepSb.append(specCount).append("\t");
            pepSb.append(ionCount).append("\t");
            pepSb.append(each.getScanNum()).append("\t");
            pepSb.append(each.getChargeState()).append("\t");
            pepSb.append(each.getFileName()).append("\n");

            tmpCount++;
          }

          if (tmpCount <= 0) {
            continue;
          }

          double ratioTotalSum = 0;

          //System.out.println("==========" + (double)normIntSumArr[1] +   " " + normIntSumArr[1]+ " " + tmpCount + " " + normIntCountArr[1]);
          //		for(long l : intSumArr) {
          for (int i = 0; i < intSumArr.length; i++) {
            //proteinSb.append( CensusHelper.format.format((double)intSumArr[i]/tmpCount) ).append("\t")
            //	.append( CensusHelper.format.format((double)normIntSumArr[i]/tmpCount) ).append("==\t");
            //if(protein.getLocus().equals("IPI00396680.3"))

            if (intCountArr[i] > 0) {
              proteinSb.append(CensusHelper.format.format((double) intSumArr[i] / intCountArr[i])).append("\t");
              proteinCalcSb.append(CensusHelper.format.format((double) intSumArr[i] / intCountArr[i])).append("\t");
            } else {
              proteinSb.append("0\t");
              proteinCalcSb.append("0\t");
            }

            if (normIntCountArr[i] > 0) {
              proteinSb.append(CensusHelper.format.format((double) normIntSumArr[i] / normIntCountArr[i])).append("\t");
              proteinCalcSb.append(CensusHelper.format.format((double) normIntSumArr[i] / normIntCountArr[i])).append("\t");
            } else {
              proteinSb.append("0\t");
              proteinCalcSb.append("0\t");

            }

          }

                    /*
                     for(double d : ratioSumArr) {
                     ratioTotalSum += d;
                     }


                     if(ratioSumArr.length>1) {
                     for(int i=0;i<ratioSumArr.length;i++) {

                     if(ratioTotalSum<=0)
                     proteinSb.append("NA\t");
                     else
                     proteinSb.append( CensusHelper.format.format(ratioSumArr[i]/tmpCount) ).append("\t");

                     }
                     }
                     else proteinSb.append("NA\t");
                     *}

                     */


          proteinSb.append(protein.getDescription());
          proteinSb.append("\n");
          List<String> redundantProteinRowList = new ArrayList<>();
          for(ChroProtein redundantProt : protein.getRedunList())
          {
            StringBuilder redundantProteinSB = new StringBuilder();
            redundantProteinSB.append("P\t").append(protein.getLocus()).append("\t");
            redundantProteinSB.append(protein.getSpectrumCount()).append("\t");
            redundantProteinSB.append(peptideList.size()).append("\t");
            redundantProteinSB.append(proteinCalcSb);
            redundantProteinSB.append(redundantProt.getDescription()).append("\n");
            redundantProteinRowList.add(redundantProteinSB.toString());
            redundantProteinRowList.add("\n");

          }


          percent += eachSeg;
          if (null != aJProgressBar) {
            aJProgressBar.setValue((int) percent);
          }

          if (pepSb.length() <= 0) {
            continue;
          }
          if(numUniquePeptides<rParam.getMinimumNumberOfUniquePeptides()) continue;
        //  if(numPeptides<rParam.getMinimumPeptidePerProtein()) continue;
          if(sequenceSet.size()<rParam.getMinimumPeptidePerProtein()) continue;



          //result.append(proteinSb.toString());
          //result.append(pepSb.toString());
          for(String redun: redundantProteinRowList)
          {
            p.print(redun);
          }
          p.print(proteinSb.toString());
          p.print(pepSb.toString());

        }

      } else if ("mean".equals(type)) { //not done yet

        //System.out.println("typedddddd.." + type);
        long[] sumIntArr = new long[rawIntList.length];

        int tcount = 0;
        for (TLongArrayList norm : rawIntList) {

          long[] arr = norm.toNativeArray();

          for (long l : arr) {
            //   System.out.println(l);
            sumIntArr[tcount] += l;
          }

          tcount++;

        }

        double averageInt = 0;
        for (long l : sumIntArr) {
          averageInt += l;
          //System.out.println(l);
        }

        averageInt = averageInt / sumIntArr.length;

        //   System.out.println("--------" + averageInt);
                /*
                 System.out.println("===" + tcount);
                 for(int i=0;i<avgIntArr.length;i++) {
                 System.out.println(avgIntArr[i]);

                 }*/
        for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
          ChroProtein protein = proItr.next();

          int proLength = -1;
          if (!"N/A".equals(protein.getLength())) {
            proLength = Integer.parseInt(protein.getLength());
          }

          String[] specArr = protein.getSpectrumCount().split(",");

          isSameGroup = true;
          StringBuffer proteinSb = new StringBuffer();

          List<ChroPeptide> peptideList = protein.getPeptideList();
          int numUniquePeptides =0;
          int numPeptides = 0;
          proteinSb.append("P\t");
          proteinSb.append(protein.getLocus()).append("\t");
          proteinSb.append(protein.getSpectrumCount()).append("\t");  //
          proteinSb.append(protein.getSeqCount()).append("\t");
          proteinSb.append(protein.getSeqCoverage()).append("\t");
          proteinSb.append(protein.getLength()).append("\t");
          proteinSb.append(protein.getMolWt()).append("\t");
          proteinSb.append(protein.getPI()).append("\t");
          proteinSb.append(peptideList.size()).append("\t"); //bug??



          double devSum = 0;

          StringBuffer pepSb = new StringBuffer();

          int tmpCount = 0;

          long[] intSumArr = new long[massMonitorList.size()];
          long[] normIntSumArr = new long[massMonitorList.size()];
          int[] normIntCountArr = new int[massMonitorList.size()];
          int[] intCountArr = new int[massMonitorList.size()];
          double[] ratioSumArr = new double[massMonitorList.size()];
          List<long[]> normList = new ArrayList<>();
          Set<String>  sequenceSet = new HashSet<>();
          if ("PEPTIDE".equals(outlierLevel)) { //find outliers among same peptides with different charge statss or different salt steps for same protien

            Hashtable<String, TMTOutlierModel> outlierHt = new Hashtable<String, TMTOutlierModel>();

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              String seq = each.getSequence();

              TMTOutlierModel om = outlierHt.get(seq);
              if (null == om) {
                om = new TMTOutlierModel(massMonitorList.size());
                om.addPeptide(each);
                outlierHt.put(seq, om);
              } else {
                om.addPeptide(each);
              }

            }

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              long[] larr = each.getTotalIntArr();
              double avgPepInt = 0;
              if (rParam.getPurityThreshold() > 0 && each.getTmtPurity() < rParam.getPurityThreshold()) {
                continue;
              }
              if (each.getSignalNoise() < rParam.getSignalToNoiseThreshold()) {
                continue;
              }
              //System.out.println("====" + each.getSequence());
              long tmpSum = 0;


              for (long l : larr) {
                tmpSum += l;
              }

              if (tmpSum <= intensityThreshold) {
                continue;
              }
              sequenceSet.add(each.getSequence());

              TMTOutlierModel om = outlierHt.get(each.getSequence());

              double[] referencePattern = om.getReferencePattern(intensityThreshold);
              double[] correctedPattern = new double[referencePattern.length];
              double correctedSum = 0;
              for (int i = 0; i < larr.length; i++) {
                correctedPattern[i] = (int) (larr[i] * averageInt / sumIntArr[i]);
                correctedSum += correctedPattern[i];
              }

              for (int i = 0; i < correctedPattern.length; i++) {
                correctedPattern[i] = correctedPattern[i] / correctedSum;

              }

              SimpleRegression sr = new SimpleRegression();

              for (int i = 0; i < correctedPattern.length; i++) {
                sr.addData(referencePattern[i], correctedPattern[i]);

                //   System.out.println(referencePattern[i] + " " +  correctedPattern[i]);
              }

              //double reg = sr.getR();
              if (sr.getR() < 0.2) {
                continue;
              }
              // else System.out.print("===" + sr.getR() + " ");
              if(each.isUnique())numUniquePeptides++;
              if(rParam.isIsUniquePeptide() && !each.isUnique())continue;
              numPeptides++;
              pepSb.append("S\t");
              pepSb.append(each.isUnique() ? "U" : "");
              pepSb.append("\t");
              pepSb.append(each.getSequence());
              pepSb.append("\t");

              long[] normArr = new long[larr.length];

              for (int i = 0; i < larr.length; i++) {
                pepSb.append(larr[i]).append("\t");

                //System.out.println(avgIntArr[i] + " aaaaaaaa" + sumIntArr[i]);
                int correctedValue = (int) (larr[i] * averageInt / sumIntArr[i]);
                //System.out.println(larr[i] + " " + averageInt + " " + sumIntArr[i] + " " + correctedValue + " " + i);

                avgPepInt += correctedValue;

                if (larr[i] > 0) {

                  intSumArr[i] += larr[i];
                  intCountArr[i]++;
                  normIntSumArr[i] += correctedValue;
                  avgPepInt += correctedValue;
                  normArr[i] = (long) correctedValue;
                  pepSb.append(correctedValue).append("\t");

                  normIntCountArr[i]++;

                } else {
                  pepSb.append("0").append("\t");
                }
              }

              avgPepInt = avgPepInt / normArr.length;

              for (int i = 0; i < normArr.length; i++) {
                if (avgPepInt > 0) {
                  ratioSumArr[i] += (double) normArr[i] / avgPepInt;
                }
              }

                            /*
                             for(int i=0;i<normArr.length;i++)
                             if(avgPepInt<=0)
                             pepSb.append("0\t");
                             else {
                             pepSb.append(CensusHelper.format.format((double)normArr[i]/avgPepInt) ).append("\t");

                             //System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
                             }*/
              String seq = each.getSequence();
              String ionSeq = seq.concat(Integer.toString(each.getChargeState()));
              int specCount = seqCountMap.get(seq);
              int ionCount = ionCountMap.get(ionSeq);
              double purity = each.getTmtPurity();
              double signal = each.getSignalNoise();
              String key = seq+each.getFileName()+each.getChargeState();
              String localScore = localMap.get(key);
              if(localScore == null) localScore = "NA";

              pepSb.append(localScore).append("\t");
              pepSb.append(purity).append("\t");
              pepSb.append(signal).append("\t");
              pepSb.append(specCount).append("\t");
              pepSb.append(ionCount).append("\t");
              pepSb.append(each.getScanNum()).append("\t");
              pepSb.append(each.getChargeState()).append("\t");
              pepSb.append(each.getFileName()).append("\t");
              pepSb.append(each.getPtmIndex()).append("\t");
              pepSb.append(each.getPtmIndexProtein()).append("\n");

              tmpCount++;
            }

          } else if ("PROTEIN".equals(outlierLevel)) {    //find outliers amount peptides in same protein

            double[] referencePattern = new double[massMonitorList.size()];

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();

              long[] larr = each.getTotalIntArr();
              //double avgPepInt = 0;

              long tmpSum = 0;
              for (long l : larr) {
                tmpSum += l;
              }

              //double tmpAvg = tmpSum/larr.length;
              if (tmpSum <= 0 || tmpSum <= intensityThreshold) {
                continue;
              }

              for (int i = 0; i < larr.length; i++) {
                referencePattern[i] += larr[i] / (double) tmpSum;

                //System.out.println("----" + i + " " + referencePattern[i] + " " + larr[i] + " " + tmpSum);
              }
              //System.out.println("");
            }

            for (int i = 0; i < referencePattern.length; i++) {
              referencePattern[i] = referencePattern[i] / peptideList.size();

              // System.out.println("====" + referencePattern[i] + " " + peptideList.size());
            }


                        /*
                         for(double d:referencePattern) {
                         System.out.print("==" + d + " " );
                         }
                         */

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              long[] larr = each.getTotalIntArr();
              double avgPepInt = 0;

              if (rParam.getPurityThreshold() > 0 && each.getTmtPurity() < rParam.getPurityThreshold()) {
                continue;
              }
              if (each.getSignalNoise() < rParam.getSignalToNoiseThreshold()) {
                continue;
              }
              //System.out.println("====" + each.getSequence());
              long tmpSum = 0;
              for (long l : larr) {
                tmpSum += l;
              }

              if (tmpSum <= intensityThreshold) ;
//                                continue;
              sequenceSet.add(each.getSequence());
              double[] correctedPattern = new double[referencePattern.length];
              double correctedSum = 0;
              for (int i = 0; i < larr.length; i++) {
                correctedPattern[i] = (int) (larr[i] * averageInt / sumIntArr[i]);
                correctedSum += correctedPattern[i];
              }

              for (int i = 0; i < correctedPattern.length; i++) {
                correctedPattern[i] = correctedPattern[i] / correctedSum;

              }

              SimpleRegression sr = new SimpleRegression();

              for (int i = 0; i < correctedPattern.length; i++) {
                sr.addData(referencePattern[i], correctedPattern[i]);
              }
              if(each.isUnique())numUniquePeptides++;
              if(rParam.isIsUniquePeptide() && !each.isUnique())continue;
              numPeptides++;
              pepSb.append("S\t");
              pepSb.append(each.isUnique() ? "U" : "");
              pepSb.append("\t");
              pepSb.append(each.getSequence());
              pepSb.append("\t");

              long[] normArr = new long[larr.length];

              for (int i = 0; i < larr.length; i++) {
                pepSb.append(larr[i]).append("\t");

                //System.out.println(avgIntArr[i] + " aaaaaaaa" + sumIntArr[i]);
                int correctedValue = (int) (larr[i] * averageInt / sumIntArr[i]);
                //System.out.println(larr[i] + " " + averageInt + " " + sumIntArr[i] + " " + correctedValue + " " + i);

                avgPepInt += correctedValue;

                if (larr[i] > 0) {

                  intSumArr[i] += larr[i];
                  intCountArr[i]++;
                  normIntSumArr[i] += correctedValue;
                  avgPepInt += correctedValue;
                  normArr[i] = (long) correctedValue;
                  pepSb.append(correctedValue).append("\t");

                  normIntCountArr[i]++;

                } else {
                  pepSb.append("0").append("\t");
                }
              }

              avgPepInt = avgPepInt / normArr.length;

              for (int i = 0; i < normArr.length; i++) {
                if (avgPepInt > 0) {
                  ratioSumArr[i] += (double) normArr[i] / avgPepInt;
                }
              }

                            /*
                             for(int i=0;i<normArr.length;i++)
                             if(avgPepInt<=0)
                             pepSb.append("0\t");
                             else {
                             pepSb.append(CensusHelper.format.format((double)normArr[i]/avgPepInt) ).append("\t");

                             //System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
                             }*/
              String seq = each.getSequence();
              String ionSeq = seq.concat(Integer.toString(each.getChargeState()));
              int specCount = seqCountMap.get(seq);
              int ionCount = ionCountMap.get(ionSeq);
              double purity = each.getTmtPurity();
              double signal = each.getSignalNoise();
              String key = seq+each.getFileName()+each.getChargeState();
              String localScore = localMap.get(key);
              if(localScore == null) localScore = "NA";

              pepSb.append(localScore).append("\t");
              pepSb.append(purity).append("\t");
              pepSb.append(signal).append("\t");
              pepSb.append(specCount).append("\t");
              pepSb.append(ionCount).append("\t");
              pepSb.append(each.getScanNum()).append("\t");
              pepSb.append(each.getChargeState()).append("\t");

              pepSb.append(each.getFileName()).append("\t");
              pepSb.append(each.getPtmIndex()).append("\t");
              pepSb.append(each.getPtmIndexProtein()).append("\n");


              tmpCount++;
            }

          } else if ("DMCCLAT".equals(outlierLevel)) { //this is Dan (dmcclat)'s approaceh.  It is still under testing with some hard code

            Hashtable<String, TMTRatioOutlierModel> ratioOutlierHt = new Hashtable<String, TMTRatioOutlierModel>();

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              String seq = each.getSequence();
              long[] larr = each.getTotalIntArr();

              int zeroCount = 0;

              for (long l : larr) {
                if (l <= 0) {
                  zeroCount++;
                }
              }

              if (zeroCount > 0) {
                continue;
              }

              TMTRatioOutlierModel om = ratioOutlierHt.get(seq);
              if (null == om) {
                om = new TMTRatioOutlierModel(averageInt, sumIntArr);
                om.addPeptide(each);
                ratioOutlierHt.put(seq, om);
              } else {
                om.addPeptide(each);
              }
            }
            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              long[] larr = each.getTotalIntArr();
              if (rParam.getPurityThreshold() > 0 && each.getTmtPurity() < rParam.getPurityThreshold()) {
                continue;
              }
              if (each.getSignalNoise() < rParam.getSignalToNoiseThreshold()) {
                continue;
              }
              int zeroCount = 0;
              long tmpSum = 0;
              for (long l : larr) {
                if (l <= 0) {
                  zeroCount++;
                }
                tmpSum += l;
              }

              if (tmpSum <= intensityThreshold || zeroCount > 0) {
                continue;
              }
              sequenceSet.add(each.getSequence());

              TMTRatioOutlierModel om = ratioOutlierHt.get(each.getSequence());

              //System.out.println("====" + each.getSequence());
              boolean removeOutlier = om.isOutlier2(each, 0.05);
              if (removeOutlier) {
                continue;
              }

              // else System.out.print("===" + sr.getR() + " ");
              if(each.isUnique())numUniquePeptides++;
              if(rParam.isIsUniquePeptide() && !each.isUnique())continue;
              numPeptides++;
              pepSb.append("S\t");
              pepSb.append(each.isUnique() ? "U" : "");
              pepSb.append("\t");
              pepSb.append(each.getSequence());
              pepSb.append("\t");

              long[] normArr = new long[larr.length];
              double avgPepInt = 0;

              for (int i = 0; i < larr.length; i++) {
                pepSb.append(larr[i]).append("\t");

                //System.out.println(avgIntArr[i] + " aaaaaaaa" + sumIntArr[i]);
                int correctedValue = (int) (larr[i] * averageInt / sumIntArr[i]);
                //System.out.println("dmcclat==\t" + larr[i] + " " + averageInt + " " + sumIntArr[i] + " " + correctedValue + " " + i);

                avgPepInt += correctedValue;

                if (larr[i] > 0) {

                  intSumArr[i] += larr[i];
                  intCountArr[i]++;
                  normIntSumArr[i] += correctedValue;
                  avgPepInt += correctedValue;
                  normArr[i] = (long) correctedValue;
                  pepSb.append(correctedValue).append("\t");

                  normIntCountArr[i]++;

                } else {
                  pepSb.append("0").append("\t");
                }
              }

              avgPepInt = avgPepInt / normArr.length;

              for (int i = 0; i < normArr.length; i++) {
                if (avgPepInt > 0) {
                  ratioSumArr[i] += (double) normArr[i] / avgPepInt;
                }
              }

                            /*
                             for(int i=0;i<normArr.length;i++)
                             if(avgPepInt<=0)
                             pepSb.append("0\t");
                             else {
                             pepSb.append(CensusHelper.format.format((double)normArr[i]/avgPepInt) ).append("\t");

                             //System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
                             }*/
              String seq = each.getSequence();
              String ionSeq = seq.concat(Integer.toString(each.getChargeState()));
              int specCount = seqCountMap.get(seq);
              int ionCount = ionCountMap.get(ionSeq);
              double purity = each.getTmtPurity();
              double signal = each.getSignalNoise();
              String key = seq+each.getFileName()+each.getChargeState();
              String localScore = localMap.get(key);
              if(localScore == null) localScore = "NA";

              pepSb.append(localScore).append("\t");
              pepSb.append(purity).append("\t");
              pepSb.append(signal).append("\t");
              pepSb.append(specCount).append("\t");
              pepSb.append(ionCount).append("\t");
              pepSb.append(each.getScanNum()).append("\t");
              pepSb.append(each.getChargeState()).append("\t");
              pepSb.append(each.getFileName()).append("\t");
              pepSb.append(each.getPtmIndex()).append("\t");
              pepSb.append(each.getPtmIndexProtein()).append("\n");

              tmpCount++;
            }

          } else {  //no outlier

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();

              long[] larr = each.getTotalIntArr();
              //double avgPepInt = 0;

              long tmpSum = 0;
              for (long l : larr) {
                tmpSum += l;
              }

              //double tmpAvg = tmpSum/larr.length;
              if (tmpSum <= 0 || tmpSum <= intensityThreshold) {
                continue;
              }

            }

            for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
              ChroPeptide each = tempItr.next();
              long[] larr = each.getTotalIntArr();
              double avgPepInt = 0;
              if (rParam.getPurityThreshold() > 0 && each.getTmtPurity() < rParam.getPurityThreshold()) {
                continue;
              }
              if (each.getSignalNoise() < rParam.getSignalToNoiseThreshold()) {
                continue;
              }
              //System.out.println("====" + each.getSequence());
              long tmpSum = 0;
              for (long l : larr) {
                tmpSum += l;
              }

              if (tmpSum <= intensityThreshold) {
                continue;
              }
              sequenceSet.add(each.getSequence());

              if(each.isUnique())numUniquePeptides++;
              if(rParam.isIsUniquePeptide() && !each.isUnique())continue;
              numPeptides++;
              pepSb.append("S\t");
              pepSb.append(each.isUnique() ? "U" : "");
              pepSb.append("\t");
              pepSb.append(each.getSequence());
              pepSb.append("\t");

              long[] normArr = new long[larr.length];

              for (int i = 0; i < larr.length; i++) {
                pepSb.append(larr[i]).append("\t");

                //System.out.println(avgIntArr[i] + " aaaaaaaa" + sumIntArr[i]);
                int correctedValue = (int) (larr[i] * averageInt / sumIntArr[i]);
                //System.out.println("==\t" + larr[i] + " " + averageInt + " " + sumIntArr[i] + " " + correctedValue + " " + i);

                avgPepInt += correctedValue;

                intCountArr[i]++;
                normIntCountArr[i]++;

                if (larr[i] > 0) {

                  intSumArr[i] += larr[i];

                  normIntSumArr[i] += correctedValue;
                  avgPepInt += correctedValue;
                  normArr[i] = (long) correctedValue;
                  pepSb.append(correctedValue).append("\t");


                } else {
                  pepSb.append("0").append("\t");
                }
              }

              normList.add(normArr);

              avgPepInt = avgPepInt / normArr.length;

              for (int i = 0; i < normArr.length; i++) {
                if (avgPepInt > 0) {
                  ratioSumArr[i] += (double) normArr[i] / avgPepInt;
                }
              }

                            /*
                             for(int i=0;i<normArr.length;i++)
                             if(avgPepInt<=0)
                             pepSb.append("0\t");
                             else {
                             pepSb.append(CensusHelper.format.format((double)normArr[i]/avgPepInt) ).append("\t");

                             //System.out.println("===" + CensusHelper.format.format((double)normArr[i]/avgPepInt) );
                             }*/
              String seq = each.getSequence();
              String ionSeq = seq.concat(Integer.toString(each.getChargeState()));
              int specCount = seqCountMap.get(seq);
              int ionCount = ionCountMap.get(ionSeq);
              double purity = each.getTmtPurity();
              double signal = each.getSignalNoise();
              String key = seq+each.getFileName()+each.getChargeState();
              String localScore = localMap.get(key);
              if(localScore == null) localScore = "NA";

              pepSb.append(localScore).append("\t");
              pepSb.append(purity).append("\t");
              pepSb.append(signal).append("\t");
              pepSb.append(specCount).append("\t");
              pepSb.append(ionCount).append("\t");
              pepSb.append(each.getScanNum()).append("\t");
              pepSb.append(each.getChargeState()).append("\t");
              pepSb.append(each.getFileName()).append("\t");
              pepSb.append(each.getPtmIndex()).append("\t");
              pepSb.append(each.getPtmIndexProtein()).append("\n");


              tmpCount++;
            }

          }

          if (tmpCount <= 0) {
            continue;
          }

          double ratioTotalSum = 0;

          //System.out.println("==========" + (double)normIntSumArr[1] +   " " + normIntSumArr[1]+ " " + tmpCount + " " + normIntCountArr[1]);
          //		for(long l : intSumArr) {
          StringBuilder proteinCalcSB = new StringBuilder();
          for (int i = 0; i < intSumArr.length; i++) {
            //proteinSb.append( CensusHelper.format.format((double)intSumArr[i]/tmpCount) ).append("\t")
            //	.append( CensusHelper.format.format((double)normIntSumArr[i]/tmpCount) ).append("==\t");
            //if(protein.getLocus().equals("IPI00396680.3"))
            // rpark.statistics.CommonStat.getMedianValue(ratioSumArr);
            if (intCountArr[i] > 0) {
              proteinSb.append(CensusHelper.format.format((double) intSumArr[i])).append("\t");
              proteinCalcSB.append(CensusHelper.format.format((double) intSumArr[i])).append("\t");
              totalIntensity += intSumArr[i];
            } else {
              proteinCalcSB.append("0\t");
              proteinSb.append("0\t");
            }


            if (intCountArr[i] > 0) {
              proteinSb.append(CensusHelper.format.format((double) intSumArr[i] / intCountArr[i])).append("\t");
              proteinCalcSB.append(CensusHelper.format.format((double) intSumArr[i] / intCountArr[i])).append("\t");
            } else {
              proteinCalcSB.append("0\t");

              proteinSb.append("0\t");
            }

            if (normIntCountArr[i] > 0) {
              proteinSb.append(CensusHelper.format.format((double) normIntSumArr[i])).append("\t");
              proteinCalcSB.append(CensusHelper.format.format((double) normIntSumArr[i])).append("\t");
            } else {
              proteinCalcSB.append("0\t");

              proteinSb.append("0\t");
            }
            if (normIntCountArr[i] > 0) {
              proteinCalcSB.append(CensusHelper.format.format((double) normIntSumArr[i] / normIntCountArr[i])).append("\t");
              proteinSb.append(CensusHelper.format.format((double) normIntSumArr[i] / normIntCountArr[i])).append("\t");
            } else {
              proteinCalcSB.append("0\t");
              proteinSb.append("0\t");
            }

            double[] arr = new double[normList.size()];
            for (int k = 0; k < normList.size(); k++) {
              arr[k] = normList.get(k)[i];
            }

            double median = CommonStat.getMedianValue(arr);
            proteinSb.append(median).append("\t");
            proteinCalcSB.append(median).append("\t");
            //   System.out.println(""+median);
          }

          for (double d : ratioSumArr) {
            ratioTotalSum += d;
          }


                    /*
                     if(ratioSumArr.length>1) {
                     for(int i=0;i<ratioSumArr.length;i++) {

                     if(ratioTotalSum<=0)
                     proteinSb.append("NA\t");
                     else
                     proteinSb.append( CensusHelper.format.format(ratioSumArr[i]/tmpCount) ).append("\t");

                     }
                     }
                     */
          percent += eachSeg;
          if (null != aJProgressBar) {
            aJProgressBar.setValue((int) percent);
          }

          if (pepSb.length() <= 0) {
            continue;
          }
          for(ChroProtein redun: protein.getRedunList())
          {
            p.append("P\t");
            p.append(redun.getLocus()).append("\t");
            p.append(redun.getSpectrumCount()).append("\t");  //
            p.append(redun.getSeqCount()).append("\t");
            p.append(redun.getSeqCoverage()).append("\t");
            p.append(redun.getLength()).append("\t");
            p.append(redun.getMolWt()).append("\t");
            p.append(redun.getPI()).append("\t");
            p.append(Integer.toString(peptideList.size())).append("\t"); //
            p.append(proteinCalcSB.toString());
            p.append(Double.toString(Math.log10(totalIntensity / (Double.parseDouble(protein.getLength()))))).append("\t");
            p.append(Double.toString(Double.parseDouble(protein.getLength()))).append("\t");
            p.append(redun.getDescription());
            p.append("\n");
          }

          proteinSb.append(Math.log10(totalIntensity / (Double.parseDouble(protein.getLength())))).append("\t");
          proteinSb.append(Double.parseDouble(protein.getLength())).append("\t");
          proteinSb.append(protein.getDescription());
          if(numUniquePeptides<rParam.getMinimumNumberOfUniquePeptides()) continue;
          //if(numPeptides<rParam.getMinimumPeptidePerProtein()) continue;
          if(sequenceSet.size()<rParam.getMinimumPeptidePerProtein()) continue;

          //result.append(proteinSb.toString()).append("\n");
          //result.append(pepSb.toString());

          p.println(proteinSb.toString());
          p.print(pepSb.toString());


          totalIntensity = 0;
        }
        //   System.out.println("----------------------------");

      } else {

        System.out.println("normalization type is required");
      }

            /*
             while( null != (eachLine = br.readLine()) ) {
             if(!eachLine.startsWith("S"))
             continue;

             String[] arr = eachLine.split("\t");
             list.add( Math.log( Double.parseDouble(arr[3])) );
             }
             */

    } catch (Exception e) {
      isSuccessful = false;
      errorMessage = e.getMessage();
      e.printStackTrace();
    }

   // return result.toString();
  }

  private static long getNoiseBaseline(int scanNum, String fileName) throws Exception {

    File indexFile = new File(fileName + ".index");
    IndexedFile iFile = new IndexedFile(indexFile, fileName);
    double[][] arr = CalcUtilGeneric.getSpectrumArr(iFile, scanNum);

    //         System.out.println(arr.length);
///                for(int i=0;i<arr.length;i++)
    //                     for(int j=0;j<arr[i].length;j++)
    //                           System.out.println(arr[i][j]);
    DescriptiveStatistics stat = new DescriptiveStatistics();

    for (int i = 0; i < arr[1].length; i++) {
      stat.addValue(arr[1][i]);

    }

    NormalDistribution ndist = new NormalDistribution(stat.getMean(), stat.getStandardDeviation());

    long threshold = 0;
    for (double i = 0; i < 1000; i++) {

      double pvalue = 1 - ndist.cumulativeProbability(i);
      if (pvalue < 0.1) {
        threshold = (long) i;
        break;
      }
    }

    return threshold;
  }

  private void exportITRAQSingleReport(final File file, final double intensityThreshold, boolean gui) {
    try {
      //display progress bar
      final JDialog progress = new JDialog(this);
      final String tempFolder = this.currentDirectory + File.separator;
      final ChroPeptide tempPeptide = this.currentPeptide;

      //final JPanel tempPepPanel = this.peptidePanel;
      final JProgressBar aJProgressBar = new JProgressBar(0, 100);
      //aJProgressBar.setIndeterminate(true);
      aJProgressBar.setStringPainted(true);

      Container cp = progress.getContentPane();

      //cp.setSize(1000, 500);
      JLabel jb = new JLabel("Generating report...",
        SwingConstants.CENTER);
      cp.add(jb, BorderLayout.SOUTH);
      cp.add(aJProgressBar, BorderLayout.NORTH);
      progress.setSize(500, 100);
      progress.setLocationRelativeTo(this);
      progress.pack();

      progress.setResizable(false);
      progress.setVisible(true);

      Configuration conf = Configuration.getInstance();

      Thread t = new Thread() {
        private boolean isSuccessful = true;
        private String errorMessage = "";
        private PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));

        public void run() {

          //String result = exportITRAQSingleReportCore(p, intensityThreshold, 0, aJProgressBar, errorMessage, isSuccessful, proteinList, "", "mean", null);
          exportITRAQSingleReportCore(p, intensityThreshold, 0, aJProgressBar, errorMessage, isSuccessful, proteinList, "", "mean", null);
          SwingUtilities.invokeLater(new Runnable() {
                                       public void run() {
                                         progress.setVisible(false);

                                         if (isSuccessful) {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             "Report file was successfully created.",
                                             "Report file Creation",
                                             JOptionPane.PLAIN_MESSAGE);
                                         } else {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             errorMessage,
                                             "Report file Creation",
                                             JOptionPane.ERROR_MESSAGE);
                                         }

                                       }

                                     }
          );

          //p.print(result);

          if (null != p) {
            p.close();
          }

        }

      };

      t.start();

    } catch (IOException e) {
      System.out.println("Failed to write file" + e);
      e.printStackTrace();

    }

  }

  public void exportITRAQMultipleReport(final File file) {
    try {
      //display progress bar
      final JDialog progress = new JDialog(this);
      final String tempFolder = this.currentDirectory + File.separator;
      //final ChroPeptide tempPeptide = this.currentPeptide;

      //final JPanel tempPepPanel = this.peptidePanel;
      final JProgressBar aJProgressBar = new JProgressBar(0, 100);
      final List<ReportIon> massMonitorList = conf.getReportIonList();
      //aJProgressBar.setIndeterminate(true);
      aJProgressBar.setStringPainted(true);

      Container cp = progress.getContentPane();

      //cp.setSize(1000, 500);
      JLabel jb = new JLabel("Generating report...",
        SwingConstants.CENTER);
      cp.add(jb, BorderLayout.SOUTH);
      cp.add(aJProgressBar, BorderLayout.NORTH);
      progress.setSize(500, 100);
      progress.setLocationRelativeTo(this);
      progress.pack();

      progress.setResizable(false);
      progress.setVisible(true);

      Configuration conf = Configuration.getInstance();

      Thread t = new Thread() {
        private boolean isSuccessful = true;
        private String errorMessage = "";
        private PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));

        public void run() {
          printHeader(p);
          p.println("H\tCensus msms analysis");
          p.print("H\tPLINE\tLOCUS\tSPEC_COUNT\t");

          StringBuffer result = new StringBuffer();

          for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
            ReportIon each = itr.next();

            p.print("m/z_");
            p.print(each);
            p.print("_total_int\t");
          }

          p.println("DESCRIPTION");

          p.print("H\tSLINE\tUNIQUE\tSEQUENCE\t"); //RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tSIGNAL_TO_NOISE_RATIO\tFILE_NAME");

          for (Iterator<ReportIon> itr = massMonitorList.iterator(); itr.hasNext(); ) {
            ReportIon each = itr.next();

            p.print("area_m/z_");
            p.print(each);
            p.print("");
            p.print("\t");
          }

          String massStr = massMonitorList.get(0).toString();
          for (int i = 1; i < massMonitorList.size(); i++) {
            p.print("ratio(" + massStr + "/" + massMonitorList.get(i) + ")");
          }

          p.println();

          int totalCount = 0;
          int quantifiedCount = 0;

          try {
            double eachSeg = (double) 100 / proteinList.size();
            double percent = 0;

            boolean isSameGroup = false;

            for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
              ChroProtein protein = proItr.next();

              //**  for normalization **
              int proLength = -1;
              if (!"N/A".equals(protein.getLength())) {
                proLength = Integer.parseInt(protein.getLength());
              }

              String[] specArr = protein.getSpectrumCount().split(",");

              isSameGroup = true;
              StringBuffer proteinSb = new StringBuffer();

              List<ChroPeptide> peptideList = protein.getPeptideList();
              //List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

              for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
                ChroPeptide peptide = pepItr.next();

                totalCount++;

                List dataList = peptide.getDataList();

                long[] intensitySumArr = new long[massMonitorList.size()];

                int startRange = Integer.parseInt(peptide.getStartRange());
                int endRange = Integer.parseInt(peptide.getEndRange());

                //LinearRegression reg = null;
                for (Iterator<ChroiTRAQLabelData> dataItr = dataList.iterator(); dataItr.hasNext(); ) {
                  ChroiTRAQLabelData eachData = dataItr.next();

                  int scanNum = eachData.getScanNum();
                  long[] intenArr = eachData.getIntensityArr();

                  for (int i = 0; i < intenArr.length; i++) {
                    if (scanNum >= startRange && scanNum <= endRange) {
                      intensitySumArr[i] += intenArr[i];
                    }
                  }

                }

                peptide.setTotalIntArr(intensitySumArr);

                //tempPepList.add(peptide);
              }

              proteinSb.append("P\t");
              proteinSb.append(protein.getLocus()).append("\t");
              proteinSb.append(protein.getSpectrumCount()).append("\t");

              double devSum = 0;

              StringBuffer pepSb = new StringBuffer();

              int tmpCount = 0;

              for (Iterator<ChroPeptide> tempItr = peptideList.iterator(); tempItr.hasNext(); ) {
                ChroPeptide each = tempItr.next();

                pepSb.append("S\t");
                pepSb.append(each.isUnique() ? "U" : "");
                pepSb.append("\t");
                pepSb.append(each.getSequence());
                pepSb.append("\t");

                long[] larr = each.getTotalIntArr();

                for (long l : larr) {
                  pepSb.append(l).append("\t");
                }

                for (int i = 1; i < larr.length; i++) {
                  pepSb.append(CensusHelper.format.format((double) larr[0] / larr[i])).append("\t");
                }

                                /*
                                 for(int jj=0;jj<pepRatioArr[tmpCount].length;jj++)
                                 {
                                 if(-1 == pepRatioArr[tmpCount][jj])
                                 pepSb.append("OL").append("\t");
                                 else
                                 {
                                 pepSb.append(  CensusHelper.format.format(pepRatioArr[tmpCount][jj]) ).append("\t");
                                 }
                                 }
                                 */
                pepSb.append(each.getFileName()).append("\n");

                tmpCount++;
              }


              proteinSb.append(protein.getDescription());
              proteinSb.append("\n");

              percent += eachSeg;
              aJProgressBar.setValue((int) percent);

              if (pepSb.length() <= 0) {
                continue;
              }

              result.append(proteinSb.toString());
              result.append(pepSb.toString());

            }
          } catch (Exception e) {
            isSuccessful = false;
            errorMessage = e.getMessage();
            e.printStackTrace();
          }

          SwingUtilities.invokeLater(new Runnable() {
                                       public void run() {
                                         progress.setVisible(false);

                                         if (isSuccessful) {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             "Report file was successfully created.",
                                             "Report file Creation",
                                             JOptionPane.PLAIN_MESSAGE);
                                         } else {
                                           JOptionPane.showMessageDialog(
                                             chroPanel,
                                             errorMessage,
                                             "Report file Creation",
                                             JOptionPane.ERROR_MESSAGE);
                                         }

                                       }

                                     }
          );

          p.print(result.toString());

          if (null != p) {
            p.close();
          }

        }

      };

      t.start();

    } catch (IOException e) {
      System.out.println("Failed to write file" + e);
      e.printStackTrace();

    }

  }

  //use this with dtaselect -t 0
  public void exportITRAQReport(double intensityThreshold) {

    JFileChooser choose = new JFileChooser(this.currentDirectory);
    choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Save CenSus Report File");
    choose.setApproveButtonText("Write");
    choose.addChoosableFileFilter(new SimpleFileNameFilter("txt", "CenSus Report File (*.txt)"));

    if (currentDirectory != null && !"".equals(currentDirectory)) {
      choose.setCurrentDirectory(new File(currentDirectory));
    }

    File f = new File(this.currentDirectory + File.separator + "census-out.txt");
    choose.setSelectedFile(f);

    int returnVal = choose.showOpenDialog(chroPanel);

    if (returnVal == choose.CANCEL_OPTION) {
      return;
    }

    final File file = choose.getSelectedFile();
    //isotopeFileField.setText(file.getAbsolutePath());

    if (cr.getExpType() == CensusConstants.MSMS_SPECIFIC_SINGLE_MASS) {
      exportITRAQSingleReport(file, intensityThreshold);
    } else {
      exportITRAQMultipleReport(file);
    }
  }

  //label free report  //alignment
  public void exportReport() {

    JFileChooser choose = new JFileChooser(this.currentDirectory);
    choose.setMultiSelectionEnabled(false);
    choose.setDialogTitle("Save Census Report File");
    choose.setApproveButtonText("Write");
    choose.addChoosableFileFilter(new SimpleFileNameFilter("txt", "CenSus Report File (*.txt)"));

    if (currentDirectory != null && !"".equals(currentDirectory)) {
      choose.setCurrentDirectory(new File(currentDirectory));
    }

    File f = new File(this.currentDirectory + File.separator + "census-out.txt");
    choose.setSelectedFile(f);

    int returnVal = choose.showOpenDialog(chroPanel);

    if (returnVal == choose.CANCEL_OPTION) {
      return;
    }

    final File file = choose.getSelectedFile();
    //isotopeFileField.setText(file.getAbsolutePath());

    try {
      //display progress bar
      final JDialog progress = new JDialog(this);
      final String tempFolder = this.currentDirectory + File.separator;
      final ChroPeptide tempPeptide = this.currentPeptide;

      //final JPanel tempPepPanel = this.peptidePanel;
      final JProgressBar aJProgressBar = new JProgressBar(0, 100);
      //aJProgressBar.setIndeterminate(true);
      aJProgressBar.setStringPainted(true);

      Container cp = progress.getContentPane();

      //cp.setSize(1000, 500);
      JLabel jb = new JLabel("Generating report...",
        SwingConstants.CENTER);
      cp.add(jb, BorderLayout.SOUTH);
      cp.add(aJProgressBar, BorderLayout.NORTH);
      progress.setSize(500, 100);
      progress.setLocationRelativeTo(this);
      progress.pack();

      progress.setResizable(false);
      progress.setVisible(true);

      Thread t = null;

      if (!isLabeled() && quantLevel == 1) {

        t = new Thread() {
          private boolean isSuccessful = true;
          private String errorMessage = "";
          private PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));

          public void run() {
            printHeader(p);
            p.println("H\tCensus nonlabeling analysis");

            ArrayList sampleList = cr.getSampleList();
            ArrayList fileList = cr.getFileList();

            p.print("H\tPLINE\tLOCUS\t");

                        /*
                         for(Iterator<String> itr=sampleList.iterator(); itr.hasNext(); )
                         {
                         String each = itr.next();
                         p.print(each);
                         p.print("_WEIGHTED_AVG_INT\t");
                         }
                         */
            for (int i = 1; i < sampleList.size(); i++) {
              String str = sampleList.get(i).toString();

              p.print("AVG(");
              p.print(sampleList.get(0));
              p.print("_INT");
              p.print("/" + str);
              p.print("_INT");
              p.print(")\t");

              p.print("STDEV(");
              p.print(sampleList.get(0));
              p.print("_INT");
              p.print("/" + str);
              p.print("_INT");
              p.print(")\t");

              p.print("SPEC_COUNT_NO_NORMALIZATION(");
              p.print(sampleList.get(0));
              p.print("/" + str);
              p.print(")\t");

              if (conf.getSpecCountNormal() == conf.SPEC_COUNT_NORMALIZATION1) {
                p.print("SPEC_COUNT_NORMALIZED(");
                p.print(sampleList.get(0));
                p.print("/" + str);
                p.print(")\t");
              }
            }

            p.println("DESCRIPTION");

            p.print("H\tSLINE\tUNIQUE\tSEQUENCE\t"); //RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tSAM_INT\tREF_INT\tSIGNAL_TO_NOISE_RATIO\tFILE_NAME");

            int fileSize = fileList.size();
            for (int i = 1; i < sampleList.size(); i++) {
              String str = sampleList.get(i).toString();

              p.print(sampleList.get(0));
              p.print("_INT");
              p.print("/" + str);
              p.print("_INT\t");

//                            p.print(sampleList.get(0) + "/" + str);
              //                          p.print("_STDEV\t");
            }

            //for(int i=0;i<sampleList.size();i++)
            for (Iterator<String> itr = sampleList.iterator(); itr.hasNext(); ) {
              String each = itr.next(); //sampleList.get(i).toString();

              p.print("INT(");
              p.print(each);
              p.print(")\t");

              p.print("SPEC_COUNT(");
              p.print(each);
              p.print(")\t");
            }

            p.print("FILE_NAME");
            p.println();

            StringBuffer result = new StringBuffer();

            int totalCount = 0;
            int quantifiedCount = 0;

            try {
              double eachSeg = (double) 100 / proteinList.size();
              double percent = 0;

              StringBuffer proteinSb = new StringBuffer();
              boolean isSameGroup = false;

              ArrayList<String> fList = cr.getFileList();

              /**
               * ** Calculate Spec Count Normalization ***
               */
              // two dimensional array for normalization
              double[] specSumArr = new double[fList.size()];
              for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
                ChroProtein protein = proItr.next();

                String[] specArr = protein.getSpectrumCount().split(",");
                int pLength = Integer.parseInt(protein.getLength());

                for (int i = 0; i < specArr.length; i++) {
                  specSumArr[i] += (double) Integer.parseInt(specArr[i]) / pLength;

                }
              }

              for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
                ChroProtein protein = proItr.next();

                //**  for normalization **//
                double[] normSpecC = new double[fList.size()];

                int proLength = -1;
                if (!"N/A".equals(protein.getLength())) {
                  proLength = Integer.parseInt(protein.getLength());
                }

                String[] specArr = protein.getSpectrumCount().split(",");

                //ArrayList<Double> sampleSpecCountNorm = new ArrayList<Double>();
                ArrayList<ChroXmlReader.Sample> samExpList = cr.getSampleExpList();

                if (conf.getSpecCountNormal() == conf.SPEC_COUNT_NORMALIZATION1) {
                  //double denominator = 0;
                  for (int i = 0; i < normSpecC.length; i++) {
                    if (proLength > 0) {
                      normSpecC[i] = Double.parseDouble(specArr[i]) / proLength / specSumArr[i];
                    }
                  }

                  int tmpCount = 0;

                  for (Iterator<ChroXmlReader.Sample> sItr = samExpList.iterator(); sItr.hasNext(); ) {
                    ChroXmlReader.Sample eachSam = sItr.next();
                    ArrayList<String> eList = eachSam.getExpList();

                    for (Iterator<String> eItr = eList.iterator(); eItr.hasNext(); ) {
                      String eEach = eItr.next();

                      eachSam.addNormalizedSpecC(normSpecC[tmpCount]);
                      tmpCount++;
                    }

                    //sampleSpecCountNorm.add( eachSam.getNormalizedSpecC() );
                  }
                }

                ArrayList<Double> sampleSpecCount = new ArrayList<Double>();

                String tmpSam = cr.getSampleName(fList.get(0));
                double spcSum = 0;

                int tCount = 0;
                for (int i = 0; i < specArr.length; i++) {
                  double d = Double.parseDouble(specArr[i]);

                  String sampleName = cr.getSampleName(fList.get(i));

                  if (!tmpSam.equals(sampleName)) {
                    sampleSpecCount.add(spcSum / tCount);
                    spcSum = 0;
                    tCount = 1;
                  } else {
                    tCount++;
                  }

                  spcSum += d;
                  tmpSam = sampleName;

                }

                sampleSpecCount.add(spcSum);

                if (!isSameGroup) {
                  proteinSb = new StringBuffer();
                }

                if (protein.isRedundant()) {
                  isSameGroup = true;
                  proteinSb.append("P\t");
                  proteinSb.append(protein.getLocus());
                  proteinSb.append("\n");

                  continue;
                } else {
                  isSameGroup = false;
                }

                List<ChroPeptide> peptideList = protein.getPeptideList();
                List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

                long[] sampleIntSumArr = new long[sampleList.size()];

                for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
                  ChroPeptide peptide = pepItr.next();

                  //peptide.getSpecCount();
                  totalCount++;

                  List dataList = peptide.getDataList();

                  long[] intensitySumArr = new long[fileSize];
                  //long[] samArr = new long[dataList.size()];
                  //long[] refArr = new long[samArr.length];

                  //double samIntSum = 0;
                  //double rpepSb.append(each.getSpScore()).append("\t");efIntSum = 0;
                  //int index=0;
                  //int startIndex=0;
                  //int endIndex=samArr.length-1;
                  int startRange = Integer.parseInt(peptide.getStartRange());
                  int endRange = Integer.parseInt(peptide.getEndRange());

                  //LinearRegression reg = null;
                  for (Iterator<ChroNonLabelData> dataItr = dataList.iterator(); dataItr.hasNext(); ) {
                    ChroNonLabelData eachData = dataItr.next();

                    int[] scanNumArr = eachData.getScanNumArr();
                    long[] intenArr = eachData.getIntensityArr();

                    for (int i = 0; i < scanNumArr.length; i++) {
                      if (scanNumArr[0] >= startRange && scanNumArr[0] <= endRange) {
                        intensitySumArr[i] += intenArr[i];
                      }
                    }

                  }

                  peptide.setTotalIntArr(intensitySumArr);

                  tempPepList.add(peptide);
                }

                int peptideCount = 0;
                //double ratioSum=0;

                for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
                  ChroPeptide each = tempItr.next();

                  //fix me when we use filtering
                  //if(!noFilter && each.isFilterOut())
                  //    continue;
                  //ratioSum += each.getSlope();
                  peptideCount++;
                  quantifiedCount++;

                }

                proteinSb.append("P\t");
                proteinSb.append(protein.getLocus());
                proteinSb.append("\t");

                double devSum = 0;

                StringBuffer pepSb = new StringBuffer();

                //int tmpc=0;
                long[][] sampleIntWeightStdevArr = new long[sampleList.size()][tempPepList.size()]; //standard dev for weighted avg
                long[][] sampleIntStdevArr = new long[sampleList.size()][tempPepList.size()]; //standard dev for simple avg

                //remove outlier
				/*
                                 for(Iterator<ChroPeptide> tempItr=tempPepList.iterator(); tempItr.hasNext(); )
                                 {
                                 ChroPeptide each = tempItr.next();

                                 List<String> samNameList = cr.getSampleList();
                                 String refSamName = samNameList.get(0);
                                 ChroXmlReader.Sample refSam = sampleHt.get(refSamName);

                                 for(int i=1;i<samNameList.size();i++)
                                 {
                                 //String sampleName = itr.next();
                                 String sampleName = samNameList.get(i);
                                 ChroXmlReader.Sample eachSam = sampleHt.get(sampleName);
                                 //sampleIntSum += eachSam.getSumIntensity();

                                 if(eachSam.getSumIntensity()>0)
                                 {
                                 pepSb.append( CensusHelper.format.format( (double)refSam.getAverage()/eachSam.getAverage()) ).append("\t");
                                 } else
                                 {
                                 //   pepSb.append( 0.0 ).append("\t");
                                 pepSb.append( 0.0 ).append("\t");
                                 }
                                 }

                                 }
                                 */
                //peptide ratio array before removing outlier
                double[][] pepRatioArr = new double[tempPepList.size()][sampleList.size() - 1];

                int pepCount = 0;

                for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
                  ChroPeptide each = tempItr.next();

                  //fix me when we use filtering
                  //if(noFilter || !each.isFilterOut())
                  //{
//                                        devSum += dev*dev;
/*
                                     pepSb.append("S\t");
                                     pepSb.append(each.isUnique()?"U":"");
                                     pepSb.append("\t");
                                     pepSb.append(each.getSequence());
                                     pepSb.append("\t");
                                     */
                  long[] totalIntArr = each.getTotalIntArr();
                  Hashtable<String, ChroXmlReader.Sample> sampleHt = new Hashtable<String, ChroXmlReader.Sample>();

                  for (int i = 0; i < totalIntArr.length; i++) {
                    String sampleName = cr.getSampleName(i);

                    ChroXmlReader.Sample sample = sampleHt.get(sampleName);

                    if (null == sample) {
                      sample = new ChroXmlReader.Sample(sampleName);
                      sample.addIntensity(totalIntArr[i]);
                      sampleHt.put(sampleName, sample);
                    } else {
                      sample.addIntensity(totalIntArr[i]);
                    }

                  }

                  int samIndex = 0;
                  List<String> samNameList = cr.getSampleList();
                  String refSamName = samNameList.get(0);
                  ChroXmlReader.Sample refSam = sampleHt.get(refSamName);

                  for (int i = 1; i < samNameList.size(); i++) {
                    //String sampleName = itr.next();
                    String sampleName = samNameList.get(i);
                    ChroXmlReader.Sample eachSam = sampleHt.get(sampleName);
                    //sampleIntSum += eachSam.getSumIntensity();

                    sampleIntStdevArr[samIndex][pepCount] = eachSam.getAverage() > 0 ? eachSam.getAverage() : 0;
                    samIndex++;

                    if (eachSam.getSumIntensity() > 0) {
                      double dValue = (double) refSam.getAverage() / eachSam.getAverage();
                      //pepSb.append( CensusHelper.format.format(dValue) ).append("\t");
                      pepRatioArr[pepCount][i - 1] = dValue;
                    } else {
                      //pepSb.append( 0.0 ).append("\t");
                      pepRatioArr[pepCount][i - 1] = 0;
                    }

                    int count = 0;
                    long[] intArr = new long[eachSam.getIntensityArr().size()];

                    for (Iterator<Long> itrInt = eachSam.getIntensityArr().iterator(); itrInt.hasNext(); ) {
                      //itrInt.next();
                      intArr[count++] = itrInt.next().longValue();
                    }

                  }

                  //for(int i=0;i<totalIntArr.length;i++)
                  //  pepSb.append( totalIntArr[i]>0?scientificFormat.format(totalIntArr[i]):0.0 ).append("\t");
//                                        pepSb.append(each.getFileName()).append("\n");
                  //}
                  pepCount++;
                }

                double[] averageRatioArr = new double[pepRatioArr[0].length];
                double[] stdevRatioArr = new double[pepRatioArr[0].length];

                for (int ii = 0; ii < pepRatioArr[0].length; ii++) {
                  double[] tmpArr = new double[pepRatioArr.length];
                  double[] tmpArrForStdev = new double[pepRatioArr.length];
                  double[] dArr = null;

                  for (int jj = 0; jj < pepRatioArr.length; jj++) {
                    tmpArr[jj] = pepRatioArr[jj][ii];
                    tmpArrForStdev[jj] = pepRatioArr[jj][ii];
                  }

                  if (tmpArr.length >= 3) //iterate until there is no outliers
                  {
                    while (true) {
                      dArr = edu.scripps.pms.stats.GrubbsTest.filterExcludingNegative(tmpArr, conf.getOutlierPValue());

                      int num1 = 0;
                      int num2 = 0;
                      for (double each : dArr) {
                        if (each == -1) {
                          num1++;
                        }
                      }

                      for (double each : tmpArr) {
                        if (each == -1) {
                          num2++;
                        }
                      }

                      tmpArr = dArr;

                      if (num1 == num2) {
                        break;
                      }

                    }
                  } else {
                    dArr = tmpArr;
                  }

                  int pepTmpCount = 0;
                  double tmpSum = 0;

                  for (int jj = 0; jj < dArr.length; jj++) {
                    if (-1 == dArr[jj]) {
                      pepRatioArr[jj][ii] = -1;
                      tmpArrForStdev[jj] = -1;
                    } else {
                      pepTmpCount++;
                      tmpSum += pepRatioArr[jj][ii];
                    }
                  }

                  averageRatioArr[ii] = tmpSum / pepTmpCount;
                  stdevRatioArr[ii] = STDev.getStdevWithoutNegative(tmpArrForStdev);

                }

                //if( protein.getLocus().equals("YMR303C") )
                //  System.exit(0);
                int tmpCount = 0;

                for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
                  ChroPeptide each = tempItr.next();

                  pepSb.append("S\t");
                  pepSb.append(each.isUnique() ? "U" : "");
                  pepSb.append("\t");
                  pepSb.append(each.getSequence());
                  pepSb.append("\t");

                  for (int jj = 0; jj < pepRatioArr[tmpCount].length; jj++) {
                    if (-1 == pepRatioArr[tmpCount][jj]) {
                      pepSb.append("OL").append("\t");
                    } else {
                      pepSb.append(CensusHelper.format.format(pepRatioArr[tmpCount][jj])).append("\t");
                    }
                  }

                  long[] pepIntArr = each.getTotalIntArr();
                  //each.get
                  for (int jj = 0; jj < pepIntArr.length; jj++) {
                    pepSb.append(pepIntArr[jj]).append("\t");
                  }

                  //pepSb.append(each.getSpScore()).append("\t");
                  pepSb.append(each.getFileName()).append("\n");

                  tmpCount++;
                }

//                                if(!noFilter && tempPepList.size()>3 && pValueSelect)
                //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);
                //double refSpec = Double.parseDouble( specArr[0] );
                //double specEach = Double.parseDouble( specArr[i+1] );
                Double refSpecC = sampleSpecCount.get(0);

                double refNormSpecC = samExpList.get(0).getNormalizedSpecC();

                for (int i = 0; i < averageRatioArr.length; i++) {
                  proteinSb.append(CensusHelper.format.format(averageRatioArr[i])).append("\t");
                  proteinSb.append(CensusHelper.format.format(stdevRatioArr[i])).append("\t");

                  double eachSpecC = sampleSpecCount.get(i + 1);

                  if (0 == refSpecC || 0 == eachSpecC) {
                    proteinSb.append("N/A").append("\t");
                  } else {
                    proteinSb.append(CensusHelper.format.format(refSpecC / eachSpecC)).append("\t");
                  }

                  if (conf.getSpecCountNormal() == conf.SPEC_COUNT_NORMALIZATION1) {
                    double eachNormSpecC = samExpList.get(i + 1).getNormalizedSpecC();

                    if (0 == eachNormSpecC || 0 == refNormSpecC) {
                      proteinSb.append("N/A").append("\t");
                    } else {
                      proteinSb.append(CensusHelper.format.format(refNormSpecC / eachNormSpecC)).append("\t");
                    }
                  }

                }

                                /*
                                 //average intensity
                                 for(int i=0;i<sampleIntSumArr.length;i++)
                                 {
                                 proteinSb.append("\t");
                                 //proteinSb.append( scientificFormat.format(peptideCount>0?(sampleIntSumArr[i]/peptideCount):0) );
                                 proteinSb.append( scientificFormat.format(peptideCount>0?(sampleIntSumArr[i]/peptideCount):0) );
                                 }
                                 */
                //proteinSb.append(peptideCount>0?peptideCount:"");
                //proteinSb.append("\t");
                //proteinSb.append(protein.getSpectrumCount());
                //proteinSb.append("\t");

                proteinSb.append(protein.getDescription());
                proteinSb.append("\n");

                percent += eachSeg;
                aJProgressBar.setValue((int) percent);

                if (pepSb.length() <= 0) {
                  continue;
                }

                result.append(proteinSb.toString());
                result.append(pepSb.toString());

              }

            } catch (Exception e) {
              isSuccessful = false;
              errorMessage = e.getMessage();
              e.printStackTrace();
            }

            final int finalTotalCount = totalCount;
            final int finalQuantifiedCount = quantifiedCount;

            SwingUtilities.invokeLater(new Runnable() {
                                         public void run() {
                                           progress.setVisible(false);
                                           progress.hide();

                                           if (isSuccessful) {
                                             JOptionPane.showMessageDialog(
                                               chroPanel,
                                               "Report file was successfully created."
                                                 + "\n\nTotal peptides : " + finalTotalCount
                                                 + "\nQuantified peptides : " + finalQuantifiedCount
                                                 + "\nQuantification efficiency : " + CensusHelper.format.format((double) finalQuantifiedCount / finalTotalCount * 100) + " %",
                                               "Report file Creation",
                                               JOptionPane.PLAIN_MESSAGE);
                                           } else {
                                             JOptionPane.showMessageDialog(
                                               chroPanel,
                                               errorMessage,
                                               "Report file Creation",
                                               JOptionPane.ERROR_MESSAGE);
                                           }

                                         }

                                       }
            );

            p.print("H\t");
            p.print("Total peptides\t");
            p.println(finalTotalCount);
            p.print("H\t");
            p.print("Quantified peptides\t");
            p.print(finalQuantifiedCount);
            p.print("\n");
            p.print("H\t");
            p.print("Quantification efficiency\t");
            p.print(CensusHelper.format.format((double) finalQuantifiedCount / finalTotalCount * 100));
            p.print(" %\n");
            //p.print("H\t"); p.print("Correction Factor (Ln)\t"); p.println(correctFactorValue);

            p.print(result.toString());

            if (null != p) {
              p.close();
            }

          }

        };

      }

      t.start();

    } catch (IOException e) {
      System.out.println("Failed to write file" + e);
      e.printStackTrace();

    }
  }

  //labeled report with all-none
  public void exportReportANGUI(
    final boolean noFilter,
    final boolean detSelect,
    final boolean pValueSelect,
    final boolean filterFragmentIons,
    final double detValue,
    final double pValue,
    final double correctFactorValue,
    final boolean isUniquePeptide,
    final boolean removeNegative,
    final boolean discardAN,
    final double allNoneLowerBound,
    final double allNoneUpperBound,
    final double allNoneCompositeScore,
    final int allNoneMinPeptide,
    final boolean discardUnlabeledPeptide,
    final boolean discardReverseProtein
  ) {

    exportReportANGUI(
      noFilter,
      detSelect,
      pValueSelect,
      filterFragmentIons,
      detValue,
      pValue,
      correctFactorValue,
      isUniquePeptide,
      removeNegative,
      discardAN,
      allNoneLowerBound,
      allNoneUpperBound,
      allNoneCompositeScore,
      allNoneMinPeptide,
      true,
      discardUnlabeledPeptide,
      false
    );
  }

  //labeled

  public static ReportResult runReportMultipleMs1labeling(ReportParam param,
                                                          PrintStream p
  ) throws Exception {

    ReportResult rResult = new ReportResult();
    StringBuffer result = new StringBuffer();
    StringBuffer singletonResult = new StringBuffer();
    int totalCount = 0;
    int quantifiedCount = 0;
    int quantifiedCountWithSingleton = 0;
    int redunProteinCount = 0;
    int uniqueProteinCount = 0;
    int proteinGroupCount = 0;
    HashMap<String, Double> pepValratioLM = new HashMap<>();
    HashMap<String, Double> pepValratioLH = new HashMap<>();
    HashMap<String, Double> pepValratioMH = new HashMap<>();
    ArrayList<ChroProtein> proteinList = param.getProteinList();
    Configuration conf = param.getConf();
    boolean discardAN = param.isDiscardAN();
    boolean noFilter = param.isNoFilter();

    boolean isGui = param.isIsGui();
    JProgressBar aJProgressBar = param.getAJProgressBar();
    double detValue = param.getDetValue();
    boolean filterFragmentIons = param.isFilterFragmentIons();
    double correctFactorValue = param.getCorrectFactorValue();
    boolean discardUnlabeledPeptide = param.isDiscardUnlabeledPeptide();
    boolean discardReverseProtein = param.isDiscardReverseProtein();
    boolean removeNegative = param.isRemoveNegative();
    boolean detSelect = param.isDetSelect();
    boolean isUniquePeptide = param.isIsUniquePeptide();
    boolean pValueSelect = param.isPValueSelect();
    double pValue = param.getPValue();
    double allNoneLowerBound = param.getAllNoneLowerBound();
    double allNoneUpperBound = param.getAllNoneUpperBound();
    double allNoneCompositeScore = param.getAllNoneCompositeScore();
    int allNoneMinPeptide = param.getAllNoneMinPeptide();
    double profileScore = param.getProfileScore();

    int maxSpectrumShift = param.getMaxSpectrumShift();

    printHeader(p);

    if (detSelect) {
      p.print("H\tDeterminant Factor : ");
      p.println(detValue);
    } else {
      p.println("H\tNo Determinant Factor");
    }


    if (pValueSelect) {
      p.print("H\tIterate Outlier: ");
      p.println(param.isIterateOutlier());
      p.print("H\tOutlier pValue: ");
      p.println(pValue);
    } else {
      p.println("H\tNo Outlier pValue");
    }

    if (filterFragmentIons) {
      p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
    }

    p.print("H\tDiscard reverse proteins : ");
    p.println(discardReverseProtein);
    p.print("H\tSmoothing peaks: ");
    p.println(param.isSmoothingPeaks());

    p.print("H\tCorrection Factor Value : ");
    p.println(correctFactorValue);
    p.print("H\tallNoneLowerBound : ");
    p.println(allNoneLowerBound);
    p.print("H\tallNoneUpperBound : ");
    p.println(allNoneUpperBound);
    p.print("H\tallNoneCompositeScore : ");
    p.println(allNoneCompositeScore);
    p.print("H\tallNoneMinPeptideNum: ");
    p.println(allNoneMinPeptide);
    p.print("H\tprofileScore : ");
    p.println(profileScore);
    p.print("H\tmaxScanShift : ");
    p.println(maxSpectrumShift);
    p.print("H\tUnique Peptide only : ");
    p.println(isUniquePeptide ? "true" : "false");

    //p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
    p.println("H\tPLINE\tProtein line");
    p.println("H\tSLINE\tPeptide line");

    //p.println("H\tPLINE\tLOCUS\tMEDIAN_RATIO_L_M\tMEDIAN_RATIO_L_H\tMEDIAN_RATIO_M_H\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tCOMPOSITE_RATIO\tCOMPOSITE_RATIO_STANDARD_DEVIATION\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tDESCRIPTION");
    //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO_L_M\tREV_RATIO_L_M\tDET_FACTOR_L_M\tRATIO_L_H\tREV_RATIO_L_H\tDET_FACTOR_L_H\tRATIO_M_H\tREV_RATIO_M_H\tDET_FACTOR_M_H\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
    p.println("H\tPLINE\tLOCUS\tMEDIAN_RATIO_L_M\tSTDEV_L_M\tMEDIAN_NORM_RATIO_L_M\tNORM_STDEV_L_M\tMEDIAN_RATIO_L_H\tSTDEV_L_H\tMEDIAN_NORM_RATIO_L_H\tNORM_STDEV_L_H\tMEDIAN_RATIO_M_H\tSTDEV_M_H\tMEDIAN_NORM_RATIO_M_H\tNORM_STDEV_M_H\tMEDIAN_AREA_RATIO_L_M\tMEDIAN_AREA_RATIO_L_H\tMEDIAN_AREA_RATIO_M_H\tCOMPOSITE_RATIO_L_M\tNORM_COMPOSITE_RATIO_L_M\tCOMPOSITE_RATIO_L_H\tNORM_COMPOSITE_RATIO_L_H\tCOMPOSITE_RATIO_M_H\tNORM_COMPOSITE_RATIO_M_H\tCOMPOSITE_RATIO_STDEV_L_M\tNORM_COMPOSITE_RATIO_STDEV_L_M\tCOMPOSITE_RATIO_STDEV_L_H\tNORM_COMPOSITE_RATIO_STDEV_L_H\tCOMPOSITE_RATIO_STDEV_M_H\tNORM_COMPOSITE_RATIO_STDEV_M_H\tlist_COMPOSITE_RATIO_L_M\tlist_NORM_COMPOSITE_RATIO_L_M\tlist_COMPOSITE_RATIO_L_H\tlist_NORM_COMPOSITE_RATIO_L_H\tlist_COMPOSITE_RATIO_M_H\tlist_NORM_COMPOSITE_RATIO_M_H\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
    p.print("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO_L_M\tREV_RATIO_L_M\tDET_FACTOR_L_M\tNORM_RATIO_L_M\tRATIO_L_H\tREV_RATIO_L_H\tDET_FACTOR_L_H\tNORM_RATIO_L_H\tRATIO_M_H\tREV_RATIO_M_H\tDET_FACTOR_M_H\tNORM_RATIO_M_H\t");
    p.print("PEAK_AREA_L\tPEAK_AREA_M\tPEAK_AREA_H\t");
    p.print("AREA_RATIO_L_M\tAREA_RATIO_L_H\tAREA_RATIO_M_H\t");
    p.println("XCorr\tdeltaCN\tFILE_NAME\tSCAN\tCS\tPROFILE_SCORE");

    double eachSeg = (double) 100 / proteinList.size();
    double percent = 0;

    StringBuffer proteinSb = new StringBuffer();
    boolean isSameGroup = false;

    List<ChroProtein> proteinResultList = new ArrayList<ChroProtein>();
    List<ChroProtein> proteinANList = new ArrayList<ChroProtein>();
    Hashtable<String, TDoubleArrayList> highScoreProteinHt = new Hashtable<String, TDoubleArrayList>();

    if (!discardAN) {
      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();

        List<ChroPeptide> peptideList = protein.getPeptideList();

        if (discardReverseProtein && protein.getLocus().startsWith("Rev")) // || protein.getLocus().startsWith("cont"))
        {
          continue;
        }

        ChroProtein newChroProtein = new ChroProtein();
        newChroProtein.setLocus(protein.getLocus());
        newChroProtein.setDescription(protein.getDescription());

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
          ChroPeptide peptide = pepItr.next();
          newChroProtein.addPeptide(peptide);

        }

        proteinResultList.add(newChroProtein);
      }
    }

    //Hashtable<String, ChroProtein> anProteinHt = new Hashtable<String, ChroProtein>();
    //************************************  analyze singleton peptides  *********************************************/
    for (Iterator<ChroProtein> itr = proteinResultList.iterator(); itr.hasNext(); ) {

      ChroProtein proAN = itr.next();
      //String desc = proAN.getDescription();
      if (discardReverseProtein && proAN.getLocus().startsWith("Rev")) // || proAN.getLocus().startsWith("cont"))
      {
        continue;
      }

      //anProteinHt.put(proAN.getLocus(), proAN);
      for (Iterator<ChroPeptide> itrp = proAN.getPeptideList().iterator(); itrp.hasNext(); ) {
        ChroPeptide pepAN = itrp.next();

        int peakStart = Integer.parseInt(pepAN.getStartRange());
        int peakEnd = Integer.parseInt(pepAN.getEndRange());

        List dataList = pepAN.getDataList();

        AllNoneUtil.getANScore(pepAN, dataList, peakStart, peakEnd);
      }
    }

    ///////***********   End of Analyzing singleton peptides   ******************************/
    for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
      HashSet<String> discardPeptideSet = new HashSet<>();
      ChroProtein protein = proItr.next();

      if (!noFilter && discardReverseProtein && (protein.getLocus().startsWith("Rev"))) // || protein.getLocus().startsWith("cont")))
      {
        continue;
      }
      redunProteinCount++;

      if (!isSameGroup) {
        proteinSb = new StringBuffer();
        proteinGroupCount = 0;
      }

      proteinGroupCount++;

      if (protein.isRedundant()) {
        isSameGroup = true;
        proteinSb.append("P\t");
        proteinSb.append(protein.getLocus());
        proteinSb.append("\n");

        continue;
      } else {
        isSameGroup = false;
      }

      List<ChroPeptide> peptideList = protein.getPeptideList();
      List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

      DescriptiveStatistics statIntensityValues = new DescriptiveStatistics();      //find intensity threshold outlier
      DescriptiveStatistics peakAreaRatioList = new DescriptiveStatistics();

      double INTENSITY_COMPOSITE_RATIO = 0.2;  //outlier value

//		for(Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); )
      for (int j = 0; j < peptideList.size(); j++) {
        totalCount++;
        statIntensityValues = new DescriptiveStatistics();
        for (int i = 0; i < peptideList.size(); i++) {
          if (i == j) {
            continue;
          }
          ChroPeptide tempPeptide = peptideList.get(i);
          statIntensityValues.addValue(tempPeptide.getMaxIntensity());
        }

        double median = statIntensityValues.getPercentile(50);
        double intensityThresholdComposite = INTENSITY_COMPOSITE_RATIO * median;

        ChroPeptide peptide = peptideList.get(j);

        //System.out.println(peptide.getScanNum() + "\t" + peptide.getSequence() + "\t" + protein.getLocus()) ;
        //if(param.getIntensityThreshold()>0 && (peptide.getMaxIntensity() < intensityThresholdComposite))
        if ((peptide.getMaxIntensity() < intensityThresholdComposite)) {
          peptide.setFilterOut(true);
          String key = peptide.getSequence() + "_" + peptide.getScanNum() + "_" + peptide.getFileName() + "_" + peptide.getChargeState();
          discardPeptideSet.add(key);

        } else {
          peptide.setFilterOut(false);
        }
        // System.out.println("====\t" + param.getIntensityThreshold() + "\t"+ (peptide.getMaxIntensity() < intensityThresholdComposite) + "\t" + peptide.getMaxIntensity() + "\t" + intensityThresholdComposite);


        List<ChroData> l = peptide.getDataList();

        int startRange = Integer.parseInt(peptide.getStartRange());
        int endRange = Integer.parseInt(peptide.getEndRange());
        //AllNoneUtil.getANScore(peptide, l, startRange, endRange);

        //long[] samArr = new long[l.size()];
        //long[] refArr = new long[samArr.length];
        long identifiedSamIntensity = 0;
        long identifiedRefIntensity = 0;

        double samIntSum = 0;
        double refIntSum = 0;

        int index = 0;
        int startIndex = 0;
        int endIndex = l.size() - 1;
        //int endIndex=samArr.length-1;

        int identifiedScanIndex = -1;
        double[][] dIntensityArr = new double[l.size()][];
        double[] retentionArr = new double[l.size()];

        for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext(); ) {
          ChroData data = dataItr.next();

          //samArr[index] = data.getSampleIntensity();
          //refArr[index] = data.getRefIntensity();
          dIntensityArr[index] = data.getdIntensityArr();
          retentionArr[index] = data.getRetentionTime();

          int scanTemp = data.getScanNum();
          if (identifiedScanIndex < 0 && data.getScanNum() >= peptide.getScanNum()) {
            identifiedScanIndex = index - 1;

            identifiedSamIntensity = data.getSampleIntensity();
            identifiedRefIntensity = data.getRefIntensity();

          }

          if (startRange >= scanTemp) {
            startIndex = index;
          }
          if (endRange >= scanTemp) {
            endIndex = index;
          }

          index++;
        }

        peptide.setPrecursorLightIntensity(identifiedSamIntensity);
        peptide.setPrecursorHeavyIntensity(identifiedRefIntensity);

        long[][] intensityArrToCompare = new long[dIntensityArr[0].length][dIntensityArr.length];
        for (int ii = 0; ii < dIntensityArr.length; ii++) {
          for (int iii = 0; iii < dIntensityArr[ii].length; iii++) {
            intensityArrToCompare[iii][ii] = (long) dIntensityArr[ii][iii];
          }
        }

        double[] intensitySumArr = new double[dIntensityArr[0].length];
        //samIntSum = 0;
        //refIntSum = 0;

        endIndex = (endIndex != 0) ? endIndex : (l.size() - 1);
        //int peakSize = endIndex-startIndex+1;
        //double[] retentionPeakArr = new double[peakSize];
        //double[][] peakIntensityArr = new double[dIntensityArr[0].length][peakSize];

        for (int ii = startIndex; ii <= endIndex; ii++) {
          for (int iii = 0; iii < dIntensityArr[ii].length; iii++) {
            intensitySumArr[iii] += dIntensityArr[ii][iii];
            //      peakIntensityArr[iii][ii-startIndex] = dIntensityArr[ii][iii];
            //      retentionPeakArr[ii-startIndex] = retentionArr[ii];
          }

          //samIntSum += samArr[ii];
          //refIntSum += refArr[ii];
        }

        //List<Double> regList = new ArrayList<Double>();
        LinearRegression reg = null;
        LinearRegression regrev = null;
        double bestCorr = 0;
        List<LinearRegression> regList = new ArrayList<LinearRegression>();

        //   List<Double> revRegList = new ArrayList<Double>();
        for (int ii = 0; ii < intensityArrToCompare.length; ii++) {
          for (int iii = ii + 1; iii < intensityArrToCompare.length; iii++) {
            if (ii == iii) {
              continue;
            }

            //apply smoothing
            reg = new LinearRegression(intensityArrToCompare[ii], intensityArrToCompare[iii], startIndex, endIndex, param.getMaxSpectrumShift(), true);
            reg.setName("" + (ii + 1) + "_" + (iii + 1));
            regrev = new LinearRegression(intensityArrToCompare[iii], intensityArrToCompare[ii], startIndex, endIndex, param.getMaxSpectrumShift(), true);
            regrev.setName("" + (ii + 1) + "_" + (iii + 1));

            if (reg.getCorr() > regrev.getCorr()) {
              regList.add(reg);
            } else {
              regrev.reverseRatio();
              regList.add(regrev);
            }
            //System.out.println(ii + " " + iii);
          }
        }
        peptide.setSpectraDataPoints(endIndex - startIndex + 1);

        // double slope = reg.getSlope();
        // double intercept = reg.getIntercept();
        // double slopeRev = regrev.getSlope();
        // double interceptRev = regrev.getIntercept();
        //slope = Math.exp(Math.log(slope) + correctFactorValue);
        //slopeRev = Math.exp(Math.log(slopeRev) - correctFactorValue);
        peptide.setRegList(regList);
        //peptide.setSlope(slope);
        //peptide.setCorr(reg.getCorr());
        //peptide.setSlopeRev(slopeRev);
        //peptide.setCorrRev(regrev.getCorr());

        peptide.setIntensitySum(intensitySumArr);

        //peptide.setSamIntensity(samIntSum);
        //peptide.setRefIntensity(refIntSum);
        //if(param.isUseProfileScore() && peptide.getAnCompositeScore()<param.getProfileScore()) continue;
        AllNoneUtil.getANScoreTriple(peptide);
        //   System.out.println("--------------->>" + peptide.getProfileScore());

        if (!noFilter) {
          if (detSelect) {
            for (Iterator<LinearRegression> itr = regList.iterator(); itr.hasNext(); ) {
              LinearRegression each = itr.next();

              if (each.getCorr() * each.getCorr() < detValue || each.getCorr() <= 0) {
                //      System.out.println(each.getCorr()*each.getCorr() + "\t" + detValue);
                //     remove = true;
                peptide.setSingleton(true);
                break;
              }
            }
          }
        }

        if (peptide.getProfileScore() < param.getProfileScore())
          continue;

        tempPepList.add(peptide);

      }

      if (tempPepList.size() <= 0) {
        continue;
      }
      //no outlier, and use median value

      int peptideCount = 0;

      if (tempPepList.size() <= 0) {
        continue;
      }

      double[][] ratioArr = new double[tempPepList.get(0).getRegList().size()][tempPepList.size()];
      int count = 0;
      DescriptiveStatistics[] descStatArr = new DescriptiveStatistics[3];
      for (int i = 0; i < descStatArr.length; i++) {
        descStatArr[i] = new DescriptiveStatistics();
      }

      for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
        ChroPeptide each = tempItr.next();
        if (each.isSingleton()) {
          peptideCount++;
          quantifiedCountWithSingleton++;
          continue;
        }

        // if(!noFilter && each.isFilterOut())
        //     continue;
        // if(!noFilter && (( Double.compare(each.getSlope(), Double.NaN) == 0) || each.getSlope()==0) ) continue;
        List<LinearRegression> tmpList = each.getRegList();
        //System.out.println(each.getSequence() + " " +tmpList.size());
        for (int i = 0; i < tmpList.size(); i++) {
          LinearRegression reg = tmpList.get(i);

          double slope = reg.getSlope();
          ratioArr[i][count] = Math.log(slope) / Math.log(2);

          if (slope >= 0) {
            descStatArr[i].addValue(ratioArr[i][count]);
          }
        }

        count++;

                /*
                 ratioSum += each.getSlope();
                 ratioSumRev += each.getSlopeRev();
                 double logRatio = Math.log(each.getSlope());
                 double logRatioRev = Math.log(each.getSlopeRev());
                 //System.out.println("==" + ratioSum + "\t" + logRatio + "\t"  + each.getSlope());
                 logRatioSum += logRatio;
                 logRatioSumRev += logRatioRev;
                 */
//		    System.out.println("==\t" + each.getSlope() + "\t" + logRatio);
        peptideCount++;
        quantifiedCount++;
      }

            /*
             double invLogRatio = 0;
             double invLogRatioRev = 0;
             if(peptideCount>0) {
             averageRatio = ratioSum/peptideCount;
             averageRatioRev = ratioSumRev/peptideCount;
             //		    System.out.println("---->>\t" + peptideCount + "\t" + ratioSum + "\t" + averageRatio );

             logAverageRatio = logRatioSum/peptideCount;
             logAverageRatioRev = logRatioSumRev/peptideCount;
             invLogRatio = Math.exp(logAverageRatio);
             invLogRatioRev = Math.exp(logAverageRatioRev);
             }
             */
      proteinSb.append("P\t");
      proteinSb.append(protein.getLocus());
      proteinSb.append("\t");

      StringBuffer commonProteinSb = new StringBuffer();

      for (DescriptiveStatistics each : descStatArr) {
        double median = each.getPercentile(50);
        median = Math.exp(median * Math.log(2));
        commonProteinSb.append(median).append("\t");
        double stdev = each.getStandardDeviation();
        stdev = Math.exp(stdev * Math.log(2));
        commonProteinSb.append(stdev).append("\t");
        //proteinSb.append(median).append("\t");
      }

      /////////////////////continue here
      //let's use log value
      //proteinSb.append( averageRatio>0?CensusHelper.d5format.format(averageRatio):"NA" );
      //proteinSb.append("\t");
      //proteinSb.append( averageRatioRev>0?CensusHelper.d5format.format(averageRatioRev):"NA" );
      //proteinSb.append("\t");
      //double devSum=0;
      //double devSumRev=0;
      StringBuffer pepSb = new StringBuffer();
      StringBuffer singlePepSb = new StringBuffer();

      //WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();
      //	double totalPepIntensity=0;
      //	double areaRatioSum=0;
      for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
        ChroPeptide each = tempItr.next();

        if (each.isSingleton()) {
          singlePepSb.append("&S\t");
          singlePepSb.append(each.isUnique() ? "U" : "");
          singlePepSb.append("\t");
          singlePepSb.append(each.getSequence());
          singlePepSb.append("\t");

          List<LinearRegression> tmpList = each.getRegList();
          for (int i = 0; i < tmpList.size(); i++) {
            LinearRegression reg = tmpList.get(i);

            double slope = reg.getSlope();
            //System.out.println(i + "== " + slope + " " + reg.getCorr());
            if (reg.getCorr() >= 0)
              singlePepSb.append(slope).append("\t").append(1 / slope).append("\t");
            else
              singlePepSb.append(-1).append("\t").append(-1).append("\t");

            singlePepSb.append(reg.getCorr());
            singlePepSb.append("\t");


            if (reg.getCorr() * reg.getCorr() >= 0.5) {
              if (i == 0) {
                pepValratioLM.put(each.getSequence() + slope, slope);
              } else if (i == 1) {
                pepValratioLH.put(each.getSequence() + slope, slope);
              } else {
                pepValratioMH.put(each.getSequence() + slope, slope);
              }
            }

          }

          double[] intArr = each.getIntensitySum();
          for (double d : intArr)
            singlePepSb.append(d).append("\t");

          if (intArr[0] <= 0 && intArr[1] <= 0)
            singlePepSb.append("NaN\t");
          else if (intArr[0] <= 0)
            singlePepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[1] <= 0)
            singlePepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            singlePepSb.append(intArr[0] / intArr[1]).append("\t");

          if (intArr[0] <= 0 && intArr[2] <= 0)
            singlePepSb.append("NaN\t");
          else if (intArr[0] <= 0)
            singlePepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[2] <= 0)
            singlePepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            singlePepSb.append(intArr[0] / intArr[2]).append("\t");

          if (intArr[1] <= 0 && intArr[2] <= 0)
            singlePepSb.append("NaN\t");
          else if (intArr[1] <= 0)
            singlePepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[2] <= 0)
            singlePepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            singlePepSb.append(intArr[1] / intArr[2]).append("\t");

          singlePepSb.append((null == each.getXCorr()) ? "" : each.getXCorr());
          singlePepSb.append("\t");
          singlePepSb.append((null == each.getDeltCN()) ? "" : each.getDeltCN());
          singlePepSb.append("\t");

          singlePepSb.append(each.getFileName()).append("\t");
          singlePepSb.append(each.getScanNum()).append("\t");
          singlePepSb.append(each.getChargeState());

          singlePepSb.append("\n");


        } else {
          pepSb.append("S\t");
          //    System.out.println("---" + each.getCorr());
          pepSb.append(each.isUnique() ? "U" : "");
          pepSb.append("\t");
          pepSb.append(each.getSequence());
          pepSb.append("\t");

          List<LinearRegression> tmpList = each.getRegList();
          for (int i = 0; i < tmpList.size(); i++) {
            LinearRegression reg = tmpList.get(i);

            double slope = reg.getSlope();
            //System.out.println(i + "== " + slope + " " + reg.getCorr());
            if (reg.getCorr() >= 0)
              pepSb.append(slope).append("\t").append(1 / slope).append("\t");
            else
              pepSb.append(-1).append("\t").append(-1).append("\t");

            pepSb.append(reg.getCorr());
            pepSb.append("\t");


            if (reg.getCorr() * reg.getCorr() >= 0.5) {
              if (i == 0) {
                pepValratioLM.put(each.getSequence() + slope, slope);
              } else if (i == 1) {
                pepValratioLH.put(each.getSequence() + slope, slope);
              } else {
                pepValratioMH.put(each.getSequence() + slope, slope);
              }
            }

          }

          double[] intArr = each.getIntensitySum();
          for (double d : intArr)
            pepSb.append(d).append("\t");

          if (intArr[0] <= 0 && intArr[1] <= 0)
            pepSb.append("NaN\t");
          else if (intArr[0] <= 0)
            pepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[1] <= 0)
            pepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            pepSb.append(intArr[0] / intArr[1]).append("\t");

          if (intArr[0] <= 0 && intArr[2] <= 0)
            pepSb.append("NaN\t");
          else if (intArr[0] <= 0)
            pepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[2] <= 0)
            pepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            pepSb.append(intArr[0] / intArr[2]).append("\t");

          if (intArr[1] <= 0 && intArr[2] <= 0)
            pepSb.append("NaN\t");
          else if (intArr[1] <= 0)
            pepSb.append(RelExMainFrame.MIN_RATIO_VALUE).append("\t");
          else if (intArr[2] <= 0)
            pepSb.append(RelExMainFrame.MAX_RATIO_VALUE).append("\t");
          else
            pepSb.append(intArr[1] / intArr[2]).append("\t");

          pepSb.append((null == each.getXCorr()) ? "" : each.getXCorr());
          pepSb.append("\t");
          pepSb.append((null == each.getDeltCN()) ? "" : each.getDeltCN());
          pepSb.append("\t");

          pepSb.append(each.getFileName()).append("\t");
          pepSb.append(each.getScanNum()).append("\t");
          pepSb.append(each.getChargeState()).append("\t");
          pepSb.append(each.getProfileScore());

          pepSb.append("\n");

        }
      }

      commonProteinSb.append(peptideCount).append("\t");
      commonProteinSb.append(protein.getSpectrumCount()).append("\t");
      proteinSb.append(commonProteinSb.toString());

      //System.out.println(proteinSb.toString());

      proteinSb.append(protein.getDescription());
      proteinSb.append("\n");

      //here
      List redProList = protein.getRedunList();
      for (Iterator<ChroProtein> redPItr = redProList.iterator(); redPItr.hasNext(); ) {
        ChroProtein redPro = redPItr.next();

        proteinSb.append("P\t");
        proteinSb.append(redPro.getLocus());
        proteinSb.append("\t");
        proteinSb.append(commonProteinSb.toString());
        proteinSb.append(redPro.getDescription());
        proteinSb.append("\n");

      }

      percent += eachSeg;

      if (isGui && null != aJProgressBar) {
        aJProgressBar.setValue((int) percent);
      }

      result.append(proteinSb.toString());
      result.append(pepSb.toString());
      result.append(singlePepSb.toString());

      uniqueProteinCount++;
    }

    //calculating median of ratio LM LH MH
    double[] arr1 = new double[pepValratioLM.size()];
    double[] arr2 = new double[pepValratioLH.size()];
    double[] arr3 = new double[pepValratioMH.size()];
    int a = 0, b = 0, c = 0;
    for (String key : pepValratioLM.keySet()) {
      double temp = Math.log(pepValratioLM.get(key)) / Math.log(2);
      if (Double.isNaN(temp)) {
        continue;
      }
      arr1[a] = temp;
      a++;
    }
    for (String key : pepValratioLH.keySet()) {
      double temp = Math.log(pepValratioLH.get(key)) / Math.log(2);
      if (Double.isNaN(temp)) {
        continue;
      }
      arr2[b] = temp;
      b++;
    }
    for (String key : pepValratioMH.keySet()) {
      double temp = Math.log(pepValratioMH.get(key)) / Math.log(2);
      if (Double.isNaN(temp)) {
        continue;
      }
      arr3[c] = temp;
      c++;
    }
    double median1 = CommonStat.getMedianValue(arr1);
    double median2 = CommonStat.getMedianValue(arr2);
    double median3 = CommonStat.getMedianValue(arr3);
    double[] medians = new double[3];
    medians[0] = median1;
    medians[1] = median2;
    medians[2] = median3;
    double mid = WeightedProtein.ProteinModel.getWeightedAverage(medians);

    //double logm1 = Math.log(med7ian1) / Math.log(2);
    //double logm2 = Math.log(median2) / Math.log(2);
    //double logm3 = Math.log(median3) / Math.log(2);
    // double avglog = (median1 + median2 + median3) / 3;
    double corr1 = mid - median1;
    double corr2 = mid - median2;
    double corr3 = mid - median3;
    //System.out.println("corrections "+corr1+ " "+corr2+ " "+corr3);

    StringBuffer newresult = new StringBuffer();
    int flag = 0;
    List<String> proList = new ArrayList<>();
    HashMap<String, List<Double>> normLM = new HashMap<>();
    HashMap<String, List<Double>> normLH = new HashMap<>();
    HashMap<String, List<Double>> normMH = new HashMap<>();
    List<Double> normL_M = new ArrayList<>();
    List<Double> normL_H = new ArrayList<>();
    List<Double> normM_H = new ArrayList<>();

    HashMap<String, List<Double>> peakLM = new HashMap<>();
    HashMap<String, List<Double>> peakLH = new HashMap<>();
    HashMap<String, List<Double>> peakMH = new HashMap<>();
    List<Double> peakL_M = new ArrayList<>();
    List<Double> peakL_H = new ArrayList<>();
    List<Double> peakM_H = new ArrayList<>();
    HashMap<String, List<Double>> compLM = new HashMap<>();
    HashMap<String, List<Double>> compLH = new HashMap<>();
    HashMap<String, List<Double>> compMH = new HashMap<>();
    HashMap<String, List<Double>> normcompLM = new HashMap<>();
    HashMap<String, List<Double>> normcompLH = new HashMap<>();
    HashMap<String, List<Double>> normcompMH = new HashMap<>();
    List<Double> compL_M = new ArrayList<>();
    List<Double> compL_H = new ArrayList<>();
    List<Double> compM_H = new ArrayList<>();
    List<Double> normcompL_M = new ArrayList<>();
    List<Double> normcompL_H = new ArrayList<>();
    List<Double> normcompM_H = new ArrayList<>();

    String[] temp = result.toString().split("\n");

    for (int g = 0; g < temp.length; g++) {

      if (temp[g].startsWith("S\t") || temp[g].startsWith("&S\t")) {

        StringBuffer pepSB = new StringBuffer();
        String[] words = temp[g].split("\t");
        pepSB.append(words[0] + "\t" + words[1] + "\t" + words[2] + "\t" + words[3] + "\t" + words[4] + "\t" + words[5] + "\t");
        double corrRatio1 = Math.log(Double.parseDouble(words[3])) / Math.log(2) - corr1;
        double val1 = Math.exp(corrRatio1 * Math.log(2));
        pepSB.append(Math.exp(corrRatio1 * Math.log(2)) + "\t");
        pepSB.append(words[6] + "\t" + words[7] + "\t" + words[8] + "\t");
        double corrRatio2 = Math.log(Double.parseDouble(words[6])) / Math.log(2) - corr2;
        double val2 = Math.exp(corrRatio2 * Math.log(2));
        pepSB.append(Math.exp(corrRatio2 * Math.log(2)) + "\t");
        pepSB.append(words[9] + "\t" + words[10] + "\t" + words[11] + "\t");
        double corrRatio3 = Math.log(Double.parseDouble(words[9])) / Math.log(2) - corr3;
        double val3 = Math.exp(corrRatio3 * Math.log(2));
        pepSB.append(Math.exp(corrRatio3 * Math.log(2)) + "\t");

        for (int i = 12; i < words.length - 1; i++)
          pepSB.append(words[i]).append("\t");
        pepSB.append(words[words.length - 1]).append("\n");
        //pepSB.append(words[12] + "\t" + words[13] + "\t" + words[14] + "\t" + words[15] + "\t" + words[16] + "\n");
        newresult.append(pepSB);

        normcompL_M.add(val1);
        normcompL_H.add(val2);
        normcompM_H.add(val3);

        if (words[0].startsWith("S")) {
          normL_M.add(Math.exp(corrRatio1 * Math.log(2)));
          normL_H.add(Math.exp(corrRatio2 * Math.log(2)));
          normM_H.add(Math.exp(corrRatio3 * Math.log(2)));
        }
        //put values in list to calculate median of peak area ratio of LM LH MH
        double dp1 = Double.parseDouble(words[15]);
        double dp2 = Double.parseDouble(words[16]);
        double dp3 = Double.parseDouble(words[17]);
        if (Double.isInfinite(dp1))
          dp1 = RelExMainFrame.MAX_RATIO_VALUE;
        else if (dp1 == 0)
          dp1 = RelExMainFrame.MIN_RATIO_VALUE;

        if (Double.isInfinite(dp2))
          dp2 = RelExMainFrame.MAX_RATIO_VALUE;
        else if (dp2 == 0)
          dp2 = RelExMainFrame.MIN_RATIO_VALUE;

        if (Double.isInfinite(dp3))
          dp3 = RelExMainFrame.MAX_RATIO_VALUE;
        else if (dp3 == 0)
          dp3 = RelExMainFrame.MIN_RATIO_VALUE;


        //        System.out.println("============\t" + dp1);
        peakL_M.add(Math.log(dp1) / Math.log(2));
        peakL_H.add(Math.log(dp2) / Math.log(2));
        peakM_H.add(Math.log(dp3) / Math.log(2));

        //for composite ratio calculations
        if (words[0].startsWith("S")) {
          double d1 = Double.parseDouble(words[3]);
          double d2 = Double.parseDouble(words[6]);
          double d3 = Double.parseDouble(words[9]);
          if (Double.isInfinite(d1))
            d1 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d1 == 0)
            d1 = RelExMainFrame.MIN_RATIO_VALUE;

          if (Double.isInfinite(d2))
            d2 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d2 == 0)
            d2 = RelExMainFrame.MIN_RATIO_VALUE;

          if (Double.isInfinite(d3))
            d3 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d3 == 0)
            d3 = RelExMainFrame.MIN_RATIO_VALUE;


          //      System.out.println("S==\t" + d1 + "\t" + d2 + "\t" + d3);

          compL_M.add(d1);
          compL_H.add(d2);
          compM_H.add(d3);
        } else {
          double d1 = Double.parseDouble(words[15]);
          double d2 = Double.parseDouble(words[16]);
          double d3 = Double.parseDouble(words[17]);
          if (Double.isInfinite(d1))
            d1 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d1 == 0)
            d1 = RelExMainFrame.MIN_RATIO_VALUE;

          if (Double.isInfinite(d2))
            d2 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d2 == 0)
            d2 = RelExMainFrame.MIN_RATIO_VALUE;

          if (Double.isInfinite(d3))
            d3 = RelExMainFrame.MAX_RATIO_VALUE;
          else if (d3 == 0)
            d3 = RelExMainFrame.MIN_RATIO_VALUE;


          //       System.out.println("&S==\t" + d1 + "\t" + d2 + "\t" + d3);

          compL_M.add(d1);
          compL_H.add(d2);
          compM_H.add(d3);

        }


      } else if (temp[g].startsWith("P\t")) {
        if ((g - 1 > 0)) {
          if (temp[g - 1].startsWith("S\t") || temp[g - 1].startsWith("&S\t")) {
            for (String pro : proList) {
              normLM.put(pro, normL_M);
              normLH.put(pro, normL_H);
              normMH.put(pro, normM_H);
              peakLM.put(pro, peakL_M);
              peakLH.put(pro, peakL_H);
              peakMH.put(pro, peakM_H);
              compLM.put(pro, compL_M);
              compLH.put(pro, compL_H);
              compMH.put(pro, compM_H);
              normcompLM.put(pro, normcompL_M);
              normcompLH.put(pro, normcompL_H);
              normcompMH.put(pro, normcompM_H);

            }
            normL_H = new ArrayList<>();
            normL_M = new ArrayList<>();
            normM_H = new ArrayList<>();
            peakL_H = new ArrayList<>();
            peakL_M = new ArrayList<>();
            peakM_H = new ArrayList<>();
            compL_H = new ArrayList<>();
            compL_M = new ArrayList<>();
            compM_H = new ArrayList<>();
            normcompL_H = new ArrayList<>();
            normcompL_M = new ArrayList<>();
            normcompM_H = new ArrayList<>();
            proList = new ArrayList<>();
            flag = 0;
          }

        }
        StringBuffer proSB = new StringBuffer();

        String[] words = temp[g].split("\t");
        proList.add(words[1]);
        for (int h = 0; h < words.length; h++) {
          if (h == words.length - 1) {
            proSB.append(words[h] + "\n");
          } else {
            proSB.append(words[h] + "\t");
          }

        }


        newresult.append(proSB);
      }
    }
    //to add last protein details in HashMap
    for (String pro : proList) {
      normLM.put(pro, normL_M);
      normLH.put(pro, normL_H);
      normMH.put(pro, normM_H);
      peakLM.put(pro, peakL_M);
      peakLH.put(pro, peakL_H);
      peakMH.put(pro, peakM_H);
      compLM.put(pro, compL_M);
      compLH.put(pro, compL_H);
      compMH.put(pro, compM_H);
      normcompLM.put(pro, normcompL_M);
      normcompLH.put(pro, normcompL_H);
      normcompMH.put(pro, normcompM_H);
    }
    HashMap<String, Double> normprocompLM = new HashMap<>();
    HashMap<String, Double> normprocompLH = new HashMap<>();
    HashMap<String, Double> normprocompMH = new HashMap<>();
    HashMap<String, Double> procompLM = new HashMap<>();
    HashMap<String, Double> procompLH = new HashMap<>();
    HashMap<String, Double> procompMH = new HashMap<>();
    HashMap<String, Double> procompSTDLM = new HashMap<>();
    HashMap<String, Double> procompSTDLH = new HashMap<>();
    HashMap<String, Double> procompSTDMH = new HashMap<>();
    HashMap<String, Double> pronormcompSTDLM = new HashMap<>();
    HashMap<String, Double> pronormcompSTDLH = new HashMap<>();
    HashMap<String, Double> pronormcompSTDMH = new HashMap<>();
    HashMap<String, Double> pronormLM = new HashMap<>();
    HashMap<String, Double> pronormLH = new HashMap<>();
    HashMap<String, Double> pronormMH = new HashMap<>();
    HashMap<String, Double> pronormSTDLM = new HashMap<>();
    HashMap<String, Double> pronormSTDLH = new HashMap<>();
    HashMap<String, Double> pronormSTDMH = new HashMap<>();
    HashMap<String, Double> propeakLM = new HashMap<>();
    HashMap<String, Double> propeakLH = new HashMap<>();
    HashMap<String, Double> propeakMH = new HashMap<>();
    for (String key : normLM.keySet()) {
      List<Double> temparr = normLM.get(key);
      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {


          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double median_std = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      pronormLM.put(key, median);
      pronormSTDLM.put(key, median_std);
    }
    for (String key : normLH.keySet()) {
      List<Double> temparr = normLH.get(key);
      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {
          System.out.println("invalid==\t" + temparr.get(i));
          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double median_std = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      pronormLH.put(key, median);
      pronormSTDLH.put(key, median_std);
    }
    for (String key : normMH.keySet()) {
      List<Double> temparr = normMH.get(key);
      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {
          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double median_std = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      pronormMH.put(key, median);
      pronormSTDMH.put(key, median_std);
    }

    for (String key : peakLM.keySet()) {
      List<Double> temparr = peakLM.get(key);

      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {
          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = median = Math.exp(median * Math.log(2));
      propeakLM.put(key, fmed);
    }
    for (String key : peakLH.keySet()) {
      List<Double> temparr = peakLH.get(key);
      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {
          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = median = Math.exp(median * Math.log(2));
      propeakLH.put(key, fmed);
    }
    for (String key : peakMH.keySet()) {
      List<Double> temparr = peakMH.get(key);
      List<Double> tempArr = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity")) {
          continue;
        }
        tempArr.add(temparr.get(i));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = median = Math.exp(median * Math.log(2));
      propeakMH.put(key, fmed);
    }
    for (String key : compLM.keySet()) {
      List<Double> temparr = compLM.get(key);
      List<Double> temparr2 = normcompLM.get(key);
      List<Double> tempArr = new ArrayList<>();
      List<Double> tempArr2 = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }

        tempArr.add(Math.log(temparr.get(i)) / Math.log(2));
      }
      for (int i = 0; i < temparr2.size(); i++) {
        if (String.valueOf(temparr2.get(i)).equals("NaN") || String.valueOf(temparr2.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }
        tempArr2.add(Math.log(temparr2.get(i)) / Math.log(2));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);


      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = Math.exp(median * Math.log(2));
      procompLM.put(key, fmed);
      double value = Math.exp((median - corr1) * Math.log(2));
      normprocompLM.put(key, value);
      double[] tArr = {-2, -100, 3};
      double med = CommonStat.getMedianValue(tArr);
      double aaa = Math.exp(med * Math.log(2));


      double stdval = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      double stdval2 = CommonStat.getRatioStdevValueWithMax(tempArr2, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);


      procompSTDLM.put(key, stdval);
      pronormcompSTDLM.put(key, stdval2);


    }
    for (String key : compLH.keySet()) {
      List<Double> temparr = compLH.get(key);
      List<Double> temparr2 = normcompLH.get(key);
      List<Double> tempArr = new ArrayList<>();
      List<Double> tempArr2 = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }
        tempArr.add(Math.log(temparr.get(i)) / Math.log(2));
      }
      for (int i = 0; i < temparr2.size(); i++) {
        if (String.valueOf(temparr2.get(i)).equals("NaN") || String.valueOf(temparr2.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }
        tempArr2.add(Math.log(temparr2.get(i)) / Math.log(2));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = Math.exp(median * Math.log(2));
      procompLH.put(key, fmed);
      double value = Math.exp((median - corr2) * Math.log(2));
      normprocompLH.put(key, value);

      double stdval = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      double stdval2 = CommonStat.getRatioStdevValueWithMax(tempArr2, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);


      procompSTDLH.put(key, stdval);
      pronormcompSTDLH.put(key, stdval2);

    }
    for (String key : compMH.keySet()) {
      List<Double> temparr = compMH.get(key);
      List<Double> temparr2 = normcompMH.get(key);
      List<Double> tempArr = new ArrayList<>();
      List<Double> tempArr2 = new ArrayList<>();
      for (int i = 0; i < temparr.size(); i++) {
        if (String.valueOf(temparr.get(i)).equals("NaN") || String.valueOf(temparr.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }
        tempArr.add(Math.log(temparr.get(i)) / Math.log(2));
      }
      for (int i = 0; i < temparr2.size(); i++) {
        if (String.valueOf(temparr2.get(i)).equals("NaN") || String.valueOf(temparr2.get(i)).equalsIgnoreCase("Infinity") || temparr.get(i) <= 0) {

          continue;
        }
        tempArr2.add(Math.log(temparr2.get(i)) / Math.log(2));
      }
      double[] arr = new double[tempArr.size()];
      for (int i = 0; i < tempArr.size(); i++) {
        arr[i] = tempArr.get(i);
      }
      double median = CommonStat.getMedianValue(arr);
      double fmed = Math.exp(median * Math.log(2));
      procompMH.put(key, fmed);
      double value = Math.exp((median - corr3) * Math.log(2));
      normprocompMH.put(key, value);

      double stdval = CommonStat.getRatioStdevValueWithMax(tempArr, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);
      double stdval2 = CommonStat.getRatioStdevValueWithMax(tempArr2, LOWER_BOUND_RATIO_THRESHOLD, UPPER_BOUND_RATIO_THRESHOLD);


      procompSTDMH.put(key, stdval);
      pronormcompSTDMH.put(key, stdval2);
    }
    StringBuffer finalresult = new StringBuffer();
    finalresult.append("H\tLog_Normalization_LM\t" + corr1 + "\n");
    finalresult.append("H\tLog_Normalization_LH\t" + corr2 + "\n");
    finalresult.append("H\tLog_Normalization_MH\t" + corr3 + "\n");
    finalresult.append("H\tLog_Median_LM\t" + median1 + "\n");
    finalresult.append("H\tLog_Median_LH\t" + median2 + "\n");
    finalresult.append("H\tLog_Median_MH\t" + median3 + "\n");
    String[] ftemp = newresult.toString().split("\n");
    for (int g = 0; g < ftemp.length; g++) {
      if (ftemp[g].startsWith("P\t")) {
        StringBuffer proSB = new StringBuffer();
        String[] words = ftemp[g].split("\t");
        proSB.append(words[0] + "\t" + words[1] + "\t" + words[2] + "\t" + words[3] + "\t");
        proSB.append(pronormLM.get(words[1]) + "\t");
        proSB.append(pronormSTDLM.get(words[1]) + "\t");
        proSB.append(words[4] + "\t" + words[5] + "\t");
        proSB.append(pronormLH.get(words[1]) + "\t");
        proSB.append(pronormSTDLH.get(words[1]) + "\t");
        proSB.append(words[6] + "\t" + words[7] + "\t");
        proSB.append(pronormMH.get(words[1]) + "\t");
        proSB.append(pronormSTDMH.get(words[1]) + "\t");
        proSB.append(propeakLM.get(words[1]) + "\t");
        proSB.append(propeakLH.get(words[1]) + "\t");
        proSB.append(propeakMH.get(words[1]) + "\t");
        proSB.append(procompLM.get(words[1]) + "\t");
        proSB.append(normprocompLM.get(words[1]) + "\t");
        proSB.append(procompLH.get(words[1]) + "\t");
        proSB.append(normprocompLH.get(words[1]) + "\t");
        proSB.append(procompMH.get(words[1]) + "\t");
        proSB.append(normprocompMH.get(words[1]) + "\t");
        proSB.append(procompSTDLM.get(words[1]) + "\t");
        proSB.append(pronormcompSTDLM.get(words[1]) + "\t");
        proSB.append(procompSTDLH.get(words[1]) + "\t");
        proSB.append(pronormcompSTDLH.get(words[1]) + "\t");
        proSB.append(procompSTDMH.get(words[1]) + "\t");
        proSB.append(pronormcompSTDMH.get(words[1]) + "\t");
        proSB.append(compLM.get(words[1]) + "\t");
        proSB.append(normcompLM.get(words[1]) + "\t");
        proSB.append(compLH.get(words[1]) + "\t");
        proSB.append(normcompLH.get(words[1]) + "\t");
        proSB.append(compMH.get(words[1]) + "\t");
        proSB.append(normcompMH.get(words[1]) + "\t");
        proSB.append(words[8] + "\t" + words[9] + "\t" + words[10] + "\n");
        finalresult.append(proSB);

      } else if (ftemp.length > 0 && !"".equals(ftemp[0])) {
        StringBuffer pepSB = new StringBuffer();
        String[] words = ftemp[g].split("\t");
        proList.add(words[1]);
        for (int h = 0; h < words.length; h++) {
          if (h == words.length - 1) {
            pepSB.append(words[h] + "\n");
          } else {
            pepSB.append(words[h] + "\t");
          }

        }
        finalresult.append(pepSB);
      }
    }

    rResult.setResult(finalresult);
    rResult.setSingletonResult(singletonResult);
    rResult.setTotalCount(totalCount);
    rResult.setQuantifiedCount(quantifiedCount);
    rResult.setRedunProteinCount(redunProteinCount);
    rResult.setUniqueProteinCount(uniqueProteinCount);
    rResult.setProteinGroupCount(proteinGroupCount);
    rResult.setQuantifiedCountWithSingleton(quantifiedCountWithSingleton);

    p.print("H\t");
    p.print("Total Redundant Proteins\t");
    p.println(redunProteinCount);
    p.print("H\t");
    p.print("Total Unique Proteins\t");
    p.println(uniqueProteinCount);
    p.print("H\t");
    p.print("Total peptides\t");
    p.println(totalCount);
    p.print("H\tQuantified peptides\t");
    p.print(quantifiedCount);
    p.print("\n");
    p.print("H\tQuantified peptides\t");
    p.print(quantifiedCountWithSingleton);
    p.print("\n");
    p.print("H\t");
    p.print("Quantification efficiency\t");
    p.print(CensusHelper.format.format((double) quantifiedCount / totalCount * 100));
    p.print(" %\n");
    p.print("H\tQuantification efficiency with singleton\t");
    //System.out.println(totalCount);
    // System.out.println(quantifiedCount);
    // System.out.println(quantifiedCountWithSingleton);

    p.print(CensusHelper.format.format((double) (quantifiedCount + quantifiedCountWithSingleton) / totalCount * 100));
    p.print(" %\n");
    p.print("H\t");
    p.print("Correction Factor (Ln)\t");
    p.println(correctFactorValue);

    p.print(finalresult.toString());
    //  singleP.print(singletonResult.toString());


    return rResult;

  }

  public static ReportResult runReport(ReportParam param,
                                       PrintStream p,
                                       PrintStream singleP
  ) throws Exception {

    ReportResult rResult = new ReportResult();
    StringBuffer result = new StringBuffer();
    StringBuffer singletonResult = new StringBuffer();
    int totalCount = 0;
    int quantifiedCount = 0;
    int totalSingletonCount = 0;
    int uniqueQuantifiedCount = 0;
    int uniqueTtalCount = 0;
    int redunProteinCount = 0;
    int uniqueProteinCount = 0;
    int proteinGroupCount = 0;
    //int singletonquantified_count=0;
    int uniqueSingleTonCount = 0;
    Set<String> peptideContainer = new HashSet<>();
    Set<String> uniqueTotalContainer = new HashSet<>();
    ArrayList<ChroProtein> proteinList = param.getProteinList();
    Configuration conf = param.getConf();
    boolean discardAN = param.isDiscardAN();
    boolean noFilter = param.isNoFilter();

    boolean isGui = param.isIsGui();
    JProgressBar aJProgressBar = param.getAJProgressBar();
    double detValue = param.getDetValue();
    boolean filterFragmentIons = param.isFilterFragmentIons();
    double correctFactorValue = param.getCorrectFactorValue();
    boolean discardUnlabeledPeptide = param.isDiscardUnlabeledPeptide();
    boolean discardReverseProtein = param.isDiscardReverseProtein();
    boolean removeNegative = param.isRemoveNegative();
    boolean detSelect = param.isDetSelect();
    boolean isUniquePeptide = param.isIsUniquePeptide();
    boolean pValueSelect = param.isPValueSelect();
    double pValue = param.getPValue();
    double allNoneLowerBound = param.getAllNoneLowerBound();
    double allNoneUpperBound = param.getAllNoneUpperBound();
    double allNoneCompositeScore = param.getAllNoneCompositeScore();
    int allNoneMinPeptide = param.getAllNoneMinPeptide();
    double profileScore = param.getProfileScore();

    int maxSpectrumShift = param.getMaxSpectrumShift();
    //boolean n15param = param.ge
        /*
         System.out.println(discardAN);
         System.out.println(noFilter);
         System.out.println(detValue);
         System.out.println(filterFragmentIons);
         System.out.println(correctFactorValue);
         System.out.println(discardUnlabeledPeptide);
         System.out.println(removeNegative);
         System.out.println(detSelect);
         System.out.println(isUniquePeptide);
         System.out.println(pValueSelect);
         System.out.println(pValue);
         System.out.println(allNoneLowerBound);
         System.out.println(allNoneUpperBound);
         System.out.println(allNoneCompositeScore);
         System.out.println(allNoneMinPeptide);
         */
    boolean isN15param = false;
    if (conf.getQuantType().startsWith("15N")) {
      isN15param = true;
    }
    printHeader(p);
    printHeader(singleP);

    if (detSelect) {
      p.print("H\tDeterminant Factor : ");
      p.println(detValue);
    } else {
      p.println("H\tNo Determinant Factor");
    }

    if (pValueSelect) {
      p.print("H\tIterate Outlier: ");
      p.println(param.isIterateOutlier());
      p.print("H\tOutlier pValue: ");
      p.println(pValue);
    } else {
      p.println("H\tNo Outlier pValue");
    }

    if (filterFragmentIons) {
      p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
    }

    p.print("H\tDiscard reverse proteins : ");
    p.println(discardReverseProtein);
    p.print("H\tminimum peptides per protein : ");
    p.println(param.getMinimumPeptidePerProtein());

    p.print("H\tSmoothing peaks: ");
    p.println(param.isSmoothingPeaks());

    p.print("H\tCorrection Factor Value : ");
    p.println(correctFactorValue);
    p.print("H\tallNoneLowerBound : ");
    p.println(allNoneLowerBound);
    p.print("H\tallNoneUpperBound : ");
    p.println(allNoneUpperBound);
    p.print("H\tallNoneCompositeScore : ");
    p.println(allNoneCompositeScore);
    p.print("H\tallNoneMinPeptideNum: ");
    p.println(allNoneMinPeptide);
    p.print("H\tprofileScore : ");
    p.println(profileScore);
    p.print("H\tmaxScanShift : ");
    p.println(maxSpectrumShift);

    singleP.print("H\tallNoneLowerBound : ");
    singleP.println(allNoneLowerBound);
    singleP.print("H\tallNoneUpperBound : ");
    singleP.println(allNoneUpperBound);
    singleP.print("H\tallNoneCompositeScore : ");
    singleP.println(allNoneCompositeScore);
    singleP.print("H\tallNoneMinPeptideNum: ");
    singleP.println(allNoneMinPeptide);
    singleP.print("H\tmaxScanShift : ");
    singleP.println(maxSpectrumShift);
    singleP.println("H\tPLINE\tProtein line");
    singleP.println("H\t&SLINE\tSingleton Peptide line");
    singleP.println("H\tPLINE\tLOCUS\tPEPTIDE_NUM\tSPEC_COUNT\tLIGHT_SPEC_COUNT\tHEAVY_SPEC_COUNT\tSINGLETON_STATUS\tDESCRIPTION");
    singleP.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

    p.print("H\tUnique Peptide only : ");
    p.println(isUniquePeptide ? "true" : "false");

    //p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
    p.println("H\tPLINE\tProtein line");
    p.println("H\tSLINE\tPeptide line");
    p.println("H\t&SLINE\tSingleton Peptide line");
    if (isN15param == true) {
      p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tCOMPOSITE_RATIO\tCOMPOSITE_RATIO_STANDARD_DEVIATION\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tTOTAL_PEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tFRACTIONAL_ABUND\tDESCRIPTION");
    } else {
      p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tCOMPOSITE_RATIO\tCOMPOSITE_RATIO_STANDARD_DEVIATION\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tTOTAL_PEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tDESCRIPTION");
    }
    //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
    p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREV_SLOPE_RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tPROBABILITY_SCORE\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
    p.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

    double eachSeg = (double) 100 / proteinList.size();
    double percent = 0;

    StringBuffer proteinSb = new StringBuffer();
    StringBuffer sproteinSb = new StringBuffer();
    boolean isSameGroup = false;

    List<ChroProtein> proteinResultList = new ArrayList<ChroProtein>();
    List<ChroProtein> proteinANList = new ArrayList<ChroProtein>();
    Hashtable<String, TDoubleArrayList> highScoreProteinHt = new Hashtable<String, TDoubleArrayList>();

    if (!discardAN) {
      for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
        ChroProtein protein = proItr.next();

        List<ChroPeptide> peptideList = protein.getPeptideList();

        if (discardReverseProtein && protein.getLocus().startsWith("Rev")) // || protein.getLocus().startsWith("cont"))
        {
          continue;
        }

        ChroProtein newChroProtein = new ChroProtein();
        newChroProtein.setLocus(protein.getLocus());
        newChroProtein.setDescription(protein.getDescription());

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); ) {
          ChroPeptide peptide = pepItr.next();
          newChroProtein.addPeptide(peptide);

        }

        proteinResultList.add(newChroProtein);

      }
    }

    Hashtable<String, ChroProtein> anProteinHt = new Hashtable<String, ChroProtein>();

    //************************************  analyze singleton peptides  *********************************************/
    for (Iterator<ChroProtein> itr = proteinResultList.iterator(); itr.hasNext(); ) {

      ChroProtein proAN = itr.next();
      //String desc = proAN.getDescription();
      if (discardReverseProtein && proAN.getLocus().startsWith("Rev")) // || proAN.getLocus().startsWith("cont"))
      {
        continue;
      }

      anProteinHt.put(proAN.getLocus(), proAN);

      for (Iterator<ChroPeptide> itrp = proAN.getPeptideList().iterator(); itrp.hasNext(); ) {
        ChroPeptide pepAN = itrp.next();

        int peakStart = Integer.parseInt(pepAN.getStartRange());
        int peakEnd = Integer.parseInt(pepAN.getEndRange());

        List dataList = pepAN.getDataList();

        AllNoneUtil.getANScore(pepAN, dataList, peakStart, peakEnd);

      }
    }

    ///////***********   End of Analyzing singleton peptides   ******************************/
    for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
      HashSet<String> discardPeptideSet = new HashSet<>();
      ChroProtein protein = proItr.next();

      if (!noFilter && discardReverseProtein && (protein.getLocus().startsWith("Rev"))) // || protein.getLocus().startsWith("cont")))
      {
        continue;
      }
      redunProteinCount++;

      if (!isSameGroup) {
        proteinSb = new StringBuffer();
        sproteinSb = new StringBuffer();
        proteinGroupCount = 0;
      }

      proteinGroupCount++;

      if (protein.isRedundant()) {
        isSameGroup = true;
        proteinSb.append("P\t");
        proteinSb.append(protein.getLocus());
        proteinSb.append("\n");
        sproteinSb.append("P\t");
        sproteinSb.append(protein.getLocus());
        sproteinSb.append("\n");

        continue;
      } else {
        isSameGroup = false;
      }

      List<ChroPeptide> peptideList = protein.getPeptideList();
      List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

      DescriptiveStatistics statIntensityValues = new DescriptiveStatistics();      //find intensity threshold outlier
      DescriptiveStatistics peakAreaRatioList = new DescriptiveStatistics();
      DescriptiveStatistics anPeakAreaRatioList = new DescriptiveStatistics();

      double INTENSITY_COMPOSITE_RATIO = 0.2;  //outlier value

            /*
             //  System.out.println("-----------" + compositeList);
             System.out.println("----------->" + statInLog.getN());
             System.out.println("----------->" + statInLog);

             System.out.println("-----------" + median);
             double intensityThresholdComposite = INTENSITY_COMPOSITE_RATIO * median;


             */
//		for(Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext(); )
      Set<String> sequenceSet = new HashSet<>();
      for (int j = 0; j < peptideList.size(); j++) {
        totalCount++;
        statIntensityValues = new DescriptiveStatistics();
        for (int i = 0; i < peptideList.size(); i++) {
          if (i == j) {
            continue;
          }
          ChroPeptide tempPeptide = peptideList.get(i);
          statIntensityValues.addValue(tempPeptide.getMaxIntensity());
        }

        double median = statIntensityValues.getPercentile(50);
        double intensityThresholdComposite = INTENSITY_COMPOSITE_RATIO * median;

//		    ChroPeptide peptide = pepItr.next();
        ChroPeptide peptide = peptideList.get(j);
        sequenceSet.add(peptide.getSequence());

//                    System.out.println(peptide.getScanNum() + "=========" + peptide.getAreaRatio() + " " + (peptide.getMaxIntensity() < intensityThresholdComposite) + " " + peptide.getMaxIntensity() + " " + INTENSITY_COMPOSITE_RATIO + " " + intensityThresholdComposite);
        //if(param.getIntensityThreshold()>0 && (peptide.getMaxIntensity() < intensityThresholdComposite))

        //robin revisit
                /*if ((peptide.getMaxIntensity() < intensityThresholdComposite)) {
                    peptide.setFilterOut(true);
                    String key = peptide.getSequence() + "_" + peptide.getScanNum() + "_" + peptide.getFileName() + "_" + peptide.getChargeState();
                    discardPeptideSet.add(key);



                } else {
                    peptide.setFilterOut(false);

                }
                */
        // System.out.println("====\t" + param.getIntensityThreshold() + "\t"+ (peptide.getMaxIntensity() < intensityThresholdComposite) + "\t" + peptide.getMaxIntensity() + "\t" + intensityThresholdComposite);
               /* String pepVal=peptide.getSequence()+peptide.getScanNum()+peptide.getChargeState();
                if(!uniqueTotalContainer.contains(pepVal)){
                    uniqueTtalCount++;
                    uniqueTotalContainer.add(pepVal);
                }*/

        List<ChroData> l = peptide.getDataList();

        int startRange = Integer.parseInt(peptide.getStartRange());
        int endRange = Integer.parseInt(peptide.getEndRange());
        AllNoneUtil.getANScore(peptide, l, startRange, endRange);

        long[] samArr = new long[l.size()];
        long[] refArr = new long[samArr.length];

        double samIntSum = 0;
        double refIntSum = 0;

        int index = 0;
        int startIndex = 0;
        int endIndex = samArr.length - 1;

        LinearRegression reg = null;
        LinearRegression regrev = null;

        if (conf.getQuantLevel() == 2 && filterFragmentIons) {
          int pepLength = l.get(0).getResidueLength();

          long[][] bsTempArr = new long[pepLength][samArr.length];
          long[][] ysTempArr = new long[pepLength][samArr.length];
          long[][] brTempArr = new long[pepLength][samArr.length];
          long[][] yrTempArr = new long[pepLength][samArr.length];

          int[] scanNumArr = new int[samArr.length];

          for (Iterator<ChroData> itr = l.iterator(); itr.hasNext(); ) {
            ChroData data = itr.next();

            long bsArr[] = data.getBsIntensity();
            long ysArr[] = data.getYsIntensity();
            long brArr[] = data.getBrIntensity();
            long yrArr[] = data.getYrIntensity();

            for (int i = 0; i < bsArr.length; i++) {
              bsTempArr[i][index] = bsArr[i];
              ysTempArr[i][index] = ysArr[i];
              brTempArr[i][index] = brArr[i];
              yrTempArr[i][index] = yrArr[i];
            }

            scanNumArr[index] = data.getScanNum();


            int scanTemp = data.getScanNum();
            if (startRange >= scanTemp) {
              startIndex = index;
            }
            if (endRange >= scanTemp) {
              endIndex = index;
            }

            index++;

          }

          FragIonList ionList = CalcUtil.getBestFragIons(bsTempArr, ysTempArr, brTempArr, yrTempArr, startIndex, endIndex, maxSpectrumShift);

          int tempIndex = 0;

          for (Iterator<FragIon> itr = ionList.iterator(); itr.hasNext(); ) {
            FragIon ion = itr.next();

            long[] tempSArr = ion.getSArr();
            long[] tempRArr = ion.getRArr();

            for (int i = 0; i < tempSArr.length; i++) {
              samArr[i] += tempSArr[i];
              refArr[i] += tempRArr[i];
            }

            if (tempIndex == ionList.getBestIndex()) {
              break;
            }

            tempIndex++;

          }

          samIntSum = 0;
          refIntSum = 0;

          endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);
          for (int ii = startIndex; ii <= endIndex; ii++) {
            samIntSum += samArr[ii];
            refIntSum += refArr[ii];
          }

          if (param.isSmoothingPeaks()) {
            reg = new LinearRegression(samArr, refArr, startIndex, endIndex, maxSpectrumShift, true);
            //	regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, conf.getMaxSpectrumShift(), true);
          } else {
            reg = new LinearRegression(samArr, refArr, startIndex, endIndex, maxSpectrumShift);
            //	regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, conf.getMaxSpectrumShift());
          }

          peptide.setSpectraDataPoints(endIndex - startIndex + 1);
        } else {

          for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext(); ) {
            ChroData data = dataItr.next();

            samArr[index] = data.getSampleIntensity();
            refArr[index] = data.getRefIntensity();

            int scanTemp = data.getScanNum();
            if (startRange >= scanTemp) {
              startIndex = index;
            }
            if (endRange >= scanTemp) {
              endIndex = index;
            }

            index++;
          }

          samIntSum = 0;
          refIntSum = 0;

          endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);
          for (int ii = startIndex; ii <= endIndex; ii++) {
            samIntSum += samArr[ii];
            refIntSum += refArr[ii];

          }

          if (param.isSmoothingPeaks()) {
            reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
            regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
          } else {
            reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift());
            regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift());
          }

          peptide.setSpectraDataPoints(endIndex - startIndex + 1);
        }

        double slope = reg.getSlope()  ;
        double intercept = reg.getIntercept();
        double slopeRev = regrev.getSlope()   ;
        double interceptRev = regrev.getIntercept();


        double logSlope =  slope > 0 ? (Math.log(slope) + correctFactorValue) : 0;
        double logSlopeRev =  slopeRev > 0 ? (Math.log(slopeRev) - correctFactorValue): 0;
        slope = Math.exp( logSlope);
        slopeRev = Math.exp(logSlopeRev);
        peptide.setSlope(slope);
        peptide.setCorr(reg.getCorr());
        peptide.setSlopeRev(slopeRev);
        peptide.setCorrRev(regrev.getCorr());

        peptide.setSamIntensity(samIntSum);
        peptide.setRefIntensity(refIntSum);
        if (!noFilter && param.isUseProfileScore() && peptide.getAnCompositeScore() < param.getProfileScore()) {
          continue;
        }

        //peptide.setSnRatio(snRatio);
        if (!noFilter) {
          if (discardUnlabeledPeptide) {
            double d1 = peptide.getLightMass();
            double d2 = peptide.getHeavyMass();

            if (d1 > 0 && d2 > 0 && (d1 == d2)) {
              continue;
            }
          }

          if (reg.getCorr() < 0 && removeNegative) {
            continue;
          }

          if (detSelect && detValue > reg.getCorr() * reg.getCorr()) {
            continue;
          }

          if (isUniquePeptide && !peptide.isUnique()) {
            continue;
          }

          if (peptide.getAnCompositeScore() < 0.5) {
            continue;
          }
        }

        tempPepList.add(peptide);


      }

      if (!noFilter && tempPepList.size() >= 3 && pValueSelect) {
        if (param.isIterateOutlier()) {
          //keep outlier iteration now
          //iterate until there is no outliers
          int currentSize;
          int prevSize = tempPepList.size();
          while (true) {
            //dArr = edu.scripps.pms.stats.GrubbsTest.filter(tmpArr, 0.1);
            //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);
            edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);

            currentSize = tempPepList.size();

            if (prevSize <= currentSize) {
              break;
            }

            prevSize = currentSize;
          }

        } else {
          edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);
          //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);

        }
      }

      int peptideCount = 0;
      double averageRatio = 0;
      double ratioSum = 0;
      double logRatioSum = 0;
      double logAverageRatio = 0;

      double averageRatioRev = 0;
      double ratioSumRev = 0;
      double logRatioSumRev = 0;
      double logAverageRatioRev = 0;


      for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
        ChroPeptide each = tempItr.next();

        if (!noFilter && each.isFilterOut()) {
          continue;
        }
        if (!noFilter && ((Double.compare(each.getSlope(), Double.NaN) == 0) || each.getSlope() == 0)) {
          continue;
        }

        ratioSum += each.getSlope();
        ratioSumRev += each.getSlopeRev();
        double logRatio = Math.log(each.getSlope());
        double logRatioRev = Math.log(each.getSlopeRev());
//System.out.println("==" + ratioSum + "\t" + logRatio + "\t"  + each.getSlope());
        if ((Double.compare(logRatio, Double.NaN) != 0))
          logRatioSum += logRatio;

        if ((Double.compare(logRatioRev, Double.NaN) != 0))
          logRatioSumRev += logRatioRev;

//		    System.out.println("==\t" + each.getSlope() + "\t" + logRatio);
        peptideCount++;
        //  quantifiedCount++;
        String peptideVal = each.getSequence() + each.getChargeState() + each.getScanNum() + each.getFileName();
        if (!peptideContainer.contains(peptideVal)) {
          //        uniqueQuantifiedCount++;
          peptideContainer.add(peptideVal);
        }

      }


      double invLogRatio = 0;
      double invLogRatioRev = 0;
      if (peptideCount > 0) {
        averageRatio = ratioSum / peptideCount;
        averageRatioRev = ratioSumRev / peptideCount;
//		    System.out.println("---->>\t" + peptideCount + "\t" + ratioSum + "\t" + averageRatio );

        logAverageRatio = logRatioSum / peptideCount;
        logAverageRatioRev = logRatioSumRev / peptideCount;
        invLogRatio = Math.exp(logAverageRatio);
        invLogRatioRev = Math.exp(logAverageRatioRev);
      }

      proteinSb.append("P\t");
      proteinSb.append(protein.getLocus());
      proteinSb.append("\t");
      proteinSb.append(invLogRatio > 0 ? CensusHelper.d5format.format(invLogRatio) : "NA");
      proteinSb.append("\t");
      proteinSb.append(invLogRatioRev > 0 ? CensusHelper.d5format.format(invLogRatioRev) : "NA");
      proteinSb.append("\t");
      //let's use log value
      //proteinSb.append( averageRatio>0?CensusHelper.d5format.format(averageRatio):"NA" );
      //proteinSb.append("\t");
      //proteinSb.append( averageRatioRev>0?CensusHelper.d5format.format(averageRatioRev):"NA" );
      //proteinSb.append("\t");
      sproteinSb.append("P\t");
      sproteinSb.append(protein.getLocus());
      sproteinSb.append("\t");

      double devSum = 0;
      double devSumRev = 0;

      StringBuffer pepSb = new StringBuffer();

      WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();

      double totalPepIntensity = 0;
      double areaRatioSum = 0;

      for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext(); ) {
        ChroPeptide each = tempItr.next();

        if (!noFilter && each.isFilterOut()) {
          continue;
        }

        //discard peptides with low profile score
        if (!noFilter && (each.getAnCompositeScore() < 0.5)) {
          continue;
        }

        //if(noFilter || !each.isFilterOut())
        double dev = each.getSlope() - averageRatio;
        double devRev = each.getSlopeRev() - averageRatioRev;

        if ((Double.compare(each.getSlope(), Double.NaN) != 0)) {
          devSum += dev * dev;
          devSumRev += devRev * devRev;
        }

        quantifiedCount++;

        pepSb.append("S\t");
        pepSb.append(each.isUnique() ? "U" : "");
        pepSb.append("\t");
        pepSb.append(each.getSequence());
        pepSb.append("\t");

        if (each.getCorr() < 0) {

          if ((Double.compare(each.getSlope(), Double.NaN) == 0)) {
            pepSb.append("NA\tNA\t");
          } else {
            pepSb.append(CensusHelper.d5format.format(each.getSlope()));
            pepSb.append("\t");
            pepSb.append(CensusHelper.d5format.format(each.getSlopeRev()));
            pepSb.append("\t");
          }

          pepSb.append(CensusHelper.d5format.format(each.getCorr()));
          pepSb.append("\t");
          pepSb.append(CensusHelper.d5format.format(each.getDetValue()));
        } else {
          pepSb.append(CensusHelper.d5format.format(each.getSlope()));
          pepSb.append("\t");
          pepSb.append(CensusHelper.d5format.format(each.getSlopeRev()));
          pepSb.append("\t");
          pepSb.append(CensusHelper.d5format.format(each.getCorr()));
          pepSb.append("\t");
          pepSb.append(CensusHelper.d5format.format(each.getDetValue()));
        }

//			System.out.println("====" + detSelect + "\t" + detValue + "\t" + each.getCorr() + "\t" +  (each.getCorr()*each.getCorr()));
        double corr = each.getCorr();
        int dataNum = each.getSpectraDataPoints();
        double tvalue = corr / Math.sqrt((1 - corr * corr) / (dataNum - 2));
        int df = dataNum - 1;

        // if(df<=0) continue;
        //new filters
        //if(each.getAnCompositeScore()<0.5) continue;
        //if(each.getEnrichment()<0.8) continue;
        // System.out.println(CensusHelper.d5format.format(each.getAnCompositeScore()) +" " + each.getFileName());
        double d = 1;
        if (df > 0) {
          TDistribution t = new TDistribution(df);
          d = 1 - t.cumulativeProbability(tvalue);
        }

        pepSb.append("\t");
        pepSb.append(d);
        pepSb.append("\t");
        pepSb.append((null == each.getXCorr()) ? "" : each.getXCorr());
        pepSb.append("\t");
        pepSb.append((null == each.getDeltCN()) ? "" : each.getDeltCN());
        pepSb.append("\t");
        pepSb.append(each.getSamIntensity());
        pepSb.append("\t");
        pepSb.append(each.getRefIntensity());
        pepSb.append("\t");
        pepSb.append(each.getMaxIntensity());
        pepSb.append("\t");
        double intRatio = (0 == each.getRefIntensity()) ? -1 : (each.getSamIntensity() / each.getRefIntensity());

        if (correctFactorValue != 0) {
          intRatio = (Math.log(intRatio) + correctFactorValue);
          intRatio = Math.exp(intRatio);
        }

        pepSb.append((intRatio >= 0) ? CensusHelper.d5format.format(intRatio) : "INF");

        if (each.getSamIntensity() <= 0) {
          peakAreaRatioList.addValue(-MAX_RATIO_LOG2); //statInLog.addValue(11.0);  //roughly 2000 folder based on log2
        } else if (each.getRefIntensity() <= 0) {
          peakAreaRatioList.addValue(MAX_RATIO_LOG2);
        } else {

//                        double intRatio = each.getSamIntensity()/each.getRefIntensity();
          peakAreaRatioList.addValue(Math.log(intRatio) / Math.log(2));
        }

        pepSb.append("\t");

        //15N remaining
        double enrichRatio = each.getRefIntensity() / (each.getSamIntensity() + each.getRefIntensity());
        //       pepSb.append( CensusHelper.format.format(each.getSnRatio()) );
        //       pepSb.append("\t");
        pepSb.append(CensusHelper.d5format.format(each.getAnCompositeScore()));

        pepSb.append("\t");
        pepSb.append(each.getFileName()).append("\t");
        pepSb.append(each.getScanNum()).append("\t");
        pepSb.append(each.getChargeState());

        String key = each.getSequence() + "_" + each.getScanNum() + "_" + each.getFileName() + "_" + each.getChargeState();
        discardPeptideSet.add(key);

/*
                if (each.getBestEnrichCorr() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getBestEnrichCorr());
                }

                if (each.getBestEnrichDelCN() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getBestEnrichDelCN());
                }

                if (each.getCorrOnePlus() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getCorrOnePlus());
                }

                if (each.getCorrOneMinus() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getCorrOneMinus());
                }
                */

        double enrichment = each.getEnrichment();

//System.out.print("\t" + enrichment);
        if (enrichment > 0) {
          pepSb.append("\t").append(enrichment);
        } else {
          pepSb.append("\t").append("N/A");
        }


        //each.getEnrichment()
        //System.out.println("===>>" + each.getSpecCount());
        pepSb.append("\n");

        double rsqrtLog = Math.log(each.getDetValue());
        double stdevLog = -0.84 * rsqrtLog + 0.43;
        double invStdev = Math.exp(stdevLog);

        //grap only valid ratio
        if (each.getSlope() > 0) {
          pModel.add(invStdev, each.getSlope());
          // double dev = each.getSlope()-averageRatio;
        }

        //Hashtable<String, ArrayList> highScoreProteinHt = new Hashtable<String, ArrayList>();
        TDoubleArrayList al = highScoreProteinHt.get(protein.getLocus());

        if (null == al) {
          al = new TDoubleArrayList();
          al.add(each.getAreaRatio());
          highScoreProteinHt.put(protein.getLocus(), al);
        } else {
          al.add(each.getAreaRatio());
        }

        totalPepIntensity += each.getSamIntensity();
        totalPepIntensity += each.getRefIntensity();
        areaRatioSum += each.getAreaRatio();
      }

      ////////////add all none peptides here
      ChroProtein anProtein = anProteinHt.get(protein.getLocus());
      StringBuffer singlePepSb = new StringBuffer();
      //System.out.println(anProtein.getLocus());
      //anProtein.getPeptideList()
      int singletonPeptideCount = 0;
      int singletonUpCount = 0;
      int singletonDownCount = 0;

      //  System.out.println("composite score.........................");
      // List<String> cscoreList = new ArrayList<String>();
      // Get a DescriptiveStatistics instance
      //List<CompositeScoreModel> compositeList = new ArrayList<CompositeScoreModel>();
// Compute some statistics
//double mean = stats.getMean();
//double std = stats.getStandardDeviation();
//double median = stats.getPercentile(50);

      if (!discardAN && null != anProtein) {

        for (Iterator<ChroPeptide> anPItr = anProtein.getPeptideList().iterator(); anPItr.hasNext(); ) {
          ChroPeptide anPep = anPItr.next();
          String key = anPep.getSequence() + "_" + anPep.getScanNum() + "_" + anPep.getFileName() + "_" + anPep.getChargeState();
          sequenceSet.add(anPep.getSequence());

          if (discardPeptideSet.contains(key)) {
//                            System.out.println("Peptide discarded....."+ key);
            continue;
          }
          //anPep.getDataList()
          //For singleton and regular peptide combined composite score
//                        System.out.println("==========" + anPep.getScanNum());

                  /*
                    if (noFilter || anPep.getAnCompositeScore() >= allNoneCompositeScore) {


                        if (anPep.getRefIntensity() <= 0) { //no ratio.  assign biggest possible intensity number
//                                compositeList.add(new CompositeScoreModel(11.0));
                        } else if (anPep.getSamIntensity() <= 0) { //no ratio
                            //                              compositeList.add(new CompositeScoreModel(-11.0));
                        } else {
                            //double intRatio = (0==anPep.getRefIntensity())?-1:(anPep.getSamIntensity()/anPep.getRefIntensity());
                            double intRatio = anPep.getSamIntensity() / anPep.getRefIntensity();
                            double intLogRatio = Math.log(intRatio) / Math.log(2);
                                //anPep.getMaxIntensity();

                                //statInLog.addValue(intLogRatio);
                            //                            compositeList.add(new CompositeScoreModel(intLogRatio));
                            //             System.out.println("2======\t" + intRatio + "\t" + Math.log(intRatio)/Math.log(2)  + "\t" + (anPep.getSamIntensity()+anPep.getRefIntensity()));
                        }
                    }

                    */
          // System.out.println("==========" + anPep.getScanNum() + " " + anPep.getAreaRatio()+ " " + allNoneUpperBound + " " +  anPep.getAreaRatio() + " " + allNoneLowerBound);

          if (!noFilter && (anPep.getAreaRatio() < allNoneUpperBound && anPep.getAreaRatio() > allNoneLowerBound)) {
            continue;
          }

          //System.out.println("==========" + anPep.getScanNum());
          //if( detSelect && detValue<anPep.getDetValue() )
          //    continue;
          //  System.out.println("aa");

          if (noFilter || anPep.getAnCompositeScore() >= allNoneCompositeScore) {
            //    System.out.println("bb");

            StringBuffer tmpSb = new StringBuffer();


            totalSingletonCount++;
            tmpSb.append("&S").append("\t").append(anPep.isUnique() ? "U" : "").append("\t").append(anPep.getSequence()).append("\t");


            if (anPep.getCorr() < 0) {
              tmpSb.append("0.0");
              tmpSb.append("\t");
              tmpSb.append(CensusHelper.d5format.format(anPep.getCorr()));
              tmpSb.append("\t");
              tmpSb.append("0.0");
            } else {
              tmpSb.append(CensusHelper.d5format.format(anPep.getSlope()));
              tmpSb.append("\t");
              tmpSb.append(CensusHelper.d5format.format(anPep.getCorr()));
              tmpSb.append("\t");
              tmpSb.append(CensusHelper.d5format.format(anPep.getDetValue()));
            }

            tmpSb.append("\t");
            tmpSb.append((null == anPep.getXCorr()) ? "" : anPep.getXCorr());
            tmpSb.append("\t");
            tmpSb.append((null == anPep.getDeltCN()) ? "" : anPep.getDeltCN());
            tmpSb.append("\t");
            tmpSb.append(anPep.getSamIntensity());
            tmpSb.append("\t");
            tmpSb.append(anPep.getRefIntensity());
            tmpSb.append("\t");
            tmpSb.append(anPep.getMaxIntensity());
            tmpSb.append("\t");

            double intRatio = (0 == anPep.getRefIntensity()) ? -1 : (anPep.getSamIntensity() / anPep.getRefIntensity());
            // System.out.println("============" + intRatio);

            if (correctFactorValue != 0) {
              intRatio = (Math.log(intRatio) + correctFactorValue);
              intRatio = Math.exp(intRatio);
            }

            if (anPep.getSamIntensity() <= 0) {
              anPeakAreaRatioList.addValue(-MAX_RATIO_LOG2); //statInLog.addValue(11.0);  //roughly 2000 folder based on log2
            } else if (anPep.getRefIntensity() <= 0) {
              anPeakAreaRatioList.addValue(MAX_RATIO_LOG2);
            } else {

              anPeakAreaRatioList.addValue(Math.log(intRatio) / Math.log(2));

              //  System.out.println("==========\t" + tmpIntRatio + "\t" + Math.log(tmpIntRatio)/Math.log(2));
            }

            if (anPep.getSamIntensity() >= anPep.getRefIntensity()) {
              singletonUpCount++;

            } else {
              singletonDownCount++;

            }

            if (anPep.getRefIntensity() <= 0) {
              tmpSb.append("INF");
              //  statInLog.addValue(11.0);  //roughly 2000 folder based on log2
              //     System.out.println("1======\t" + intRatio  + "\t" + (anPep.getSamIntensity()+anPep.getRefIntensity()));
            } else {
              tmpSb.append(CensusHelper.d5format.format(intRatio));
              //        statInLog.addValue( Math.log(intRatio)/Math.log(2) );

              //        System.out.println("2======\t" + intRatio + "\t" + Math.log(intRatio)/Math.log(2)  + "\t" + (anPep.getSamIntensity()+anPep.getRefIntensity()));
            }

            //                                        tmpSb.append( CensusHelper.format.format(anPep.getSamIntensity()/anPep.getRefIntensity()) );
            tmpSb.append("\t");
            tmpSb.append(CensusHelper.d5format.format(anPep.getAnCompositeScore()));
            tmpSb.append("\t");
            tmpSb.append(anPep.getFileName()).append("\t");
            tmpSb.append(anPep.getScanNum()).append("\t");
            tmpSb.append(anPep.getChargeState()).append("\n");

            //cscoreList.add((intRatio>=0)?CensusHelper.format.format(intRatio):"2000" );
            //if(intRatio>)
            //pepSb.append(tmpSb);
            singletonPeptideCount++;
            // singletonquantified_count++;
            singlePepSb.append(tmpSb);

          }
        }

      }

      // for(CompositeScoreModel c:compositeList) {
      //    System.out.println(c.getIntRatioLog());
      // }
            /*
             Collections.sort(compositeList);
             for(CompositeScoreModel c:compositeList) {
             //        System.out.println(c.getIntRatioLog());

             }*/
      //////////End of adding an peptides
      if (peptideCount > 1 && Double.compare(devSum, Double.NaN) != 0) {
        proteinSb.append(CensusHelper.d5format.format(Math.sqrt(devSum / (peptideCount - 1))));
        proteinSb.append("\t");
        double value = Math.sqrt(devSumRev / (peptideCount - 1));
        String toPrint;
        if(Double.isNaN(value) || Double.isInfinite(value)) {
          toPrint = "NA";
        }
        else {
          toPrint = CensusHelper.d5format.format(value);
        }
        proteinSb.append(toPrint);
      } else {
        proteinSb.append("NA\tNA");
      }

      proteinSb.append("\t");

      //composite score
      //proteinSb.append( Math.pow(2, peakAreaRatioList.getMean()) );
      if (singletonPeptideCount >= allNoneMinPeptide) {
        double[] anArr = anPeakAreaRatioList.getValues();
        for (double d : anArr) {
          peakAreaRatioList.addValue(d);
        }
      }

      proteinSb.append(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))));
      proteinSb.append("\t");
      proteinSb.append(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getStandardDeviation())));
      proteinSb.append("\t");

      if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
        proteinSb.append("NA\t");

      } else {
        proteinSb.append(CensusHelper.d5format.format(pModel.getStandardWeightedAverage()));
        proteinSb.append("\t");
      }

      //     System.out.println(statInLog.getMean());
      int totalPeptideNumber = 0;

      if (singletonPeptideCount < allNoneMinPeptide) {
        totalPeptideNumber = peptideCount;
      } else {
        totalPeptideNumber = peptideCount + singletonPeptideCount;
      }

      int singletonFinal = singletonPeptideCount;
      if (singletonPeptideCount < allNoneMinPeptide)
        singletonFinal = 0;


      //System.out.println("total peptide=======" + protein.getLocus() + "\t" + totalPeptideNumber + "\t" + peptideCount + "\t" + singletonPeptideCount + "\t" + singletonFinal + "\t" + allNoneMinPeptide);
      //if()

      //if (totalPeptideNumber < param.getMinimumPeptidePerProtein()) continue;
      if (sequenceSet.size() < param.getMinimumPeptidePerProtein()) continue;
      //System.out.println("total peptide=======" + protein.getLocus() + "\t" + totalPeptideNumber + "\t" + peptideCount + "\t" + singletonPeptideCount + "\t" + singletonFinal + "\t" + allNoneMinPeptide);


      proteinSb.append(invLogRatio > 0 ? CensusHelper.d5format.format(invLogRatio) : "NA");
      proteinSb.append("\t");
      proteinSb.append(invLogRatioRev > 0 ? CensusHelper.d5format.format(invLogRatioRev) : "NA");
      proteinSb.append("\t");
      proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
      proteinSb.append("\t");
      proteinSb.append(totalPeptideNumber > 0 ? totalPeptideNumber : "NA");
      proteinSb.append("\t");
      proteinSb.append(protein.getSpectrumCount());
      proteinSb.append("\t");
      proteinSb.append(protein.getLspectrumCount());
      proteinSb.append("\t");
      proteinSb.append(protein.getHspectrumCount());
      proteinSb.append("\t");
      //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
      proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
      proteinSb.append("\t");
      if (isN15param == true) {
        if (CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))).equals("�")) {
          proteinSb.append("NA\t");
        } else {
          double compositeVal = Double.parseDouble(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))));
          double fa = 100 * (1 / (1 + compositeVal));
          proteinSb.append(fa + "\t");
        }
      }
      proteinSb.append(protein.getDescription());
      proteinSb.append("\n");

      sproteinSb.append(singletonPeptideCount);
      sproteinSb.append("\t");
      sproteinSb.append(protein.getSpectrumCount());
      sproteinSb.append("\t");
      sproteinSb.append(protein.getLspectrumCount());
      sproteinSb.append("\t");
      sproteinSb.append(protein.getHspectrumCount());
      sproteinSb.append("\t");
      sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
      sproteinSb.append("\t");
      sproteinSb.append(protein.getDescription());
      sproteinSb.append("\n");

      //here
      List redProList = protein.getRedunList();
      for (Iterator<ChroProtein> redPItr = redProList.iterator(); redPItr.hasNext(); ) {
        ChroProtein redPro = redPItr.next();

        proteinSb.append("P\t");
        proteinSb.append(redPro.getLocus());
        proteinSb.append("\t");
        proteinSb.append(invLogRatio > 0 ? CensusHelper.d5format.format(invLogRatio) : "NA");
        proteinSb.append("\t");
        proteinSb.append(invLogRatioRev > 0 ? CensusHelper.d5format.format(invLogRatioRev) : "NA");
        proteinSb.append("\t");
        //proteinSb.append( averageRatio>0?CensusHelper.d5format.format(averageRatio):"NA" );
        //proteinSb.append("\t");
        //proteinSb.append( averageRatioRev>0?CensusHelper.d5format.format(averageRatioRev):"NA" );

        if (peptideCount > 1) {
          proteinSb.append(CensusHelper.d5format.format(Math.sqrt(devSum / (peptideCount - 1))));
          proteinSb.append("\t");
          proteinSb.append(CensusHelper.d5format.format(Math.sqrt(devSumRev / (peptideCount - 1))));
        } else {
          proteinSb.append("NA\tNA");
        }

        proteinSb.append("\t");
        proteinSb.append(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))));
        proteinSb.append("\t");
        proteinSb.append(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getStandardDeviation())));
        proteinSb.append("\t");

        if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
          proteinSb.append("NA\t");

        } else {
          proteinSb.append(CensusHelper.d5format.format(pModel.getStandardWeightedAverage()));
          proteinSb.append("\t");
        }

        //proteinSb.append( CensusHelper.format.format(pModel.getStandardWeightedAverage()) );
        proteinSb.append(invLogRatio > 0 ? CensusHelper.d5format.format(invLogRatio) : "NA");
        proteinSb.append("\t");
        proteinSb.append(invLogRatioRev > 0 ? CensusHelper.d5format.format(invLogRatioRev) : "NA");
        proteinSb.append("\t");
        proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
        proteinSb.append("\t");
        proteinSb.append(totalPeptideNumber > 0 ? totalPeptideNumber : "NA");
        proteinSb.append("\t");
        proteinSb.append(redPro.getSpectrumCount());
        proteinSb.append("\t");
        proteinSb.append(redPro.getLspectrumCount());
        proteinSb.append("\t");
        proteinSb.append(redPro.getHspectrumCount());
        proteinSb.append("\t");
        //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
        proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
        proteinSb.append("\t");
        if (isN15param == true) {
          if (CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))).equals("�")) {
            proteinSb.append("NA\t");
          } else {
            double compositeVal = Double.parseDouble(CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50))));
            double fa = 100 * (1 / (1 + compositeVal));
            proteinSb.append(fa + "\t");
          }
        }

        //CensusHelper.d5format.format(Math.pow(2, peakAreaRatioList.getPercentile(50)))
        proteinSb.append(redPro.getDescription());
        proteinSb.append("\n");

        sproteinSb.append("P\t");
        sproteinSb.append(redPro.getLocus());
        sproteinSb.append("\t");
        sproteinSb.append(singletonPeptideCount);
        sproteinSb.append("\t");
        sproteinSb.append(protein.getSpectrumCount());
        sproteinSb.append("\t");
        sproteinSb.append(protein.getLspectrumCount());
        sproteinSb.append("\t");
        sproteinSb.append(protein.getHspectrumCount());
        sproteinSb.append("\t");
        sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
        sproteinSb.append("\t");
        sproteinSb.append(protein.getDescription());
        sproteinSb.append("\n");

      }

      percent += eachSeg;

      if (isGui && null != aJProgressBar) {
        aJProgressBar.setValue((int) percent);
      }

      if (singletonPeptideCount < allNoneMinPeptide) {
        singlePepSb = new StringBuffer();
      }

      if (pepSb.length() <= 0 && singlePepSb.length() <= 0) {
        redunProteinCount -= proteinGroupCount;
        continue;
      }

      //if( pepSb.length()<=0 && singlePepSb.length()<=0) {
      //if( pepSb.length()<=0 && singlePepSb.length()>0) {
      if (singlePepSb.length() > 0) {
        singletonResult.append(sproteinSb.toString());
        singletonResult.append(singlePepSb.toString());
        pepSb.append(singlePepSb.toString());

        String[] uniqueOneS = singlePepSb.toString().split("\n");
        for (int u = 0; u < uniqueOneS.length; u++) {
          String[] words = uniqueOneS[u].split("\t");
          String pepVal = words[2] + words[14] + words[15];
          if (!uniqueTotalContainer.contains(pepVal)) {
            uniqueTtalCount++;
            uniqueTotalContainer.add(pepVal);
          }
        }


      }

      result.append(proteinSb.toString());
      result.append(pepSb.toString());

      uniqueProteinCount++;
    }

    rResult.setResult(result);
    rResult.setSingletonResult(singletonResult);
    rResult.setTotalCount(totalCount);
    rResult.setQuantifiedCount(quantifiedCount);
    rResult.setQuantifiedCountWithSingleton(quantifiedCount + totalSingletonCount);
    rResult.setRedunProteinCount(redunProteinCount);
    rResult.setUniqueProteinCount(uniqueProteinCount);
    rResult.setProteinGroupCount(proteinGroupCount);

    p.print("H\t");
    p.print("Total Redundant Proteins\t");
    p.println(redunProteinCount);
    p.print("H\t");
    p.print("Total Unique Proteins\t");
    p.println(uniqueProteinCount);
    p.print("H\t");
    p.print("Total peptides\t");
    p.println(totalCount);
    p.print("H\t");
    p.print("Quantified peptides\t");
    p.print(quantifiedCount);
    p.print("\n");
    p.print("H\t");
    p.print("Quantification efficiency\t");
    p.print(CensusHelper.format.format((double) quantifiedCount / totalCount * 100));
    p.print(" %\n");
    p.print("H\t");
    p.print("Quantification efficiency including singletons\t");

    //   System.out.println("==" + quantifiedCount + " " + totalSingletonCount+ " "  + totalCount);

    p.print(CensusHelper.format.format(((double) (totalSingletonCount + quantifiedCount) / totalCount) * 100));
    p.print(" %\n");
    p.print("H\t");
    p.print("Correction Factor (Ln)\t");
    p.println(correctFactorValue);

    singleP.print("H\t");
    singleP.print("Total Redundant Proteins\t");
    singleP.println(redunProteinCount);
    singleP.print("H\t");
    singleP.print("Total Unique Proteins\t");
    singleP.println(uniqueProteinCount);
    singleP.print("H\t");
    singleP.print("Total peptides\t");
    singleP.println(totalCount);

    p.print(result.toString());
    singleP.print(singletonResult.toString());

    return rResult;

  }

  public static ReportResult runReportDIALF(ReportParam param,
                                            PrintStream p,
                                            PrintStream singleP
  ) throws Exception {

    ReportResult rResult = new ReportResult();
    StringBuffer result = new StringBuffer();
    StringBuffer singletonResult = new StringBuffer();
    int totalCount = 0;
    int quantifiedCount = 0;
    int totalSingletonCount = 0;
    int uniqueQuantifiedCount = 0;
    int uniqueTtalCount = 0;
    int redunProteinCount = 0;
    int uniqueProteinCount = 0;
    int proteinGroupCount = 0;
    //int singletonquantified_count=0;
    int uniqueSingleTonCount = 0;
    Set<String> peptideContainer = new HashSet<>();
    Set<String> uniqueTotalContainer = new HashSet<>();
    ArrayList<ChroProtein> proteinList = param.getProteinList();
    Configuration conf = param.getConf();
    boolean discardAN = param.isDiscardAN();
    boolean noFilter = param.isNoFilter();

    boolean isGui = param.isIsGui();
    JProgressBar aJProgressBar = param.getAJProgressBar();
    double detValue = param.getDetValue();
    boolean filterFragmentIons = param.isFilterFragmentIons();
    double correctFactorValue = param.getCorrectFactorValue();
    boolean discardUnlabeledPeptide = param.isDiscardUnlabeledPeptide();
    boolean discardReverseProtein = param.isDiscardReverseProtein();
    boolean removeNegative = param.isRemoveNegative();
    boolean detSelect = param.isDetSelect();
    boolean isUniquePeptide = param.isIsUniquePeptide();
    boolean pValueSelect = param.isPValueSelect();
    double pValue = param.getPValue();
    double allNoneLowerBound = param.getAllNoneLowerBound();
    double allNoneUpperBound = param.getAllNoneUpperBound();
    double allNoneCompositeScore = param.getAllNoneCompositeScore();
    int allNoneMinPeptide = param.getAllNoneMinPeptide();
    double profileScore = param.getProfileScore();

    int maxSpectrumShift = param.getMaxSpectrumShift();
    //boolean n15param = param.ge
        /*
         System.out.println(discardAN);
         System.out.println(noFilter);
         System.out.println(detValue);
         System.out.println(filterFragmentIons);
         System.out.println(correctFactorValue);
         System.out.println(discardUnlabeledPeptide);
         System.out.println(removeNegative);
         System.out.println(detSelect);
         System.out.println(isUniquePeptide);
         System.out.println(pValueSelect);
         System.out.println(pValue);
         System.out.println(allNoneLowerBound);
         System.out.println(allNoneUpperBound);
         System.out.println(allNoneCompositeScore);
         System.out.println(allNoneMinPeptide);
         */
    boolean isN15param = false;
    if (conf.getQuantType().startsWith("15N")) {
      isN15param = true;
    }
    printHeader(p);
    printHeader(singleP);

    if (detSelect) {
      p.print("H\tDeterminant Factor : ");
      p.println(detValue);
    } else {
      p.println("H\tNo Determinant Factor");
    }

    if (pValueSelect) {
      p.print("H\tIterate Outlier: ");
      p.println(param.isIterateOutlier());
      p.print("H\tOutlier pValue: ");
      p.println(pValue);
    } else {
      p.println("H\tNo Outlier pValue");
    }

    if (filterFragmentIons) {
      p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
    }

    p.print("H\tDiscard reverse proteins : ");
    p.println(discardReverseProtein);
    p.print("H\tminimum peptides per protein : ");
    p.println(param.getMinimumPeptidePerProtein());

    p.print("H\tSmoothing peaks: ");
    p.println(param.isSmoothingPeaks());

    p.print("H\tCorrection Factor Value : ");
    p.println(correctFactorValue);
    p.print("H\tallNoneLowerBound : ");
    p.println(allNoneLowerBound);
    p.print("H\tallNoneUpperBound : ");
    p.println(allNoneUpperBound);
    p.print("H\tallNoneCompositeScore : ");
    p.println(allNoneCompositeScore);
    p.print("H\tallNoneMinPeptideNum: ");
    p.println(allNoneMinPeptide);
    p.print("H\tprofileScore : ");
    p.println(profileScore);
    p.print("H\tmaxScanShift : ");
    p.println(maxSpectrumShift);

    singleP.print("H\tallNoneLowerBound : ");
    singleP.println(allNoneLowerBound);
    singleP.print("H\tallNoneUpperBound : ");
    singleP.println(allNoneUpperBound);
    singleP.print("H\tallNoneCompositeScore : ");
    singleP.println(allNoneCompositeScore);
    singleP.print("H\tallNoneMinPeptideNum: ");
    singleP.println(allNoneMinPeptide);
    singleP.print("H\tmaxScanShift : ");
    singleP.println(maxSpectrumShift);
    singleP.println("H\tPLINE\tProtein line");
    singleP.println("H\t&SLINE\tSingleton Peptide line");
    singleP.println("H\tPLINE\tLOCUS\tPEPTIDE_NUM\tSPEC_COUNT\tLIGHT_SPEC_COUNT\tHEAVY_SPEC_COUNT\tSINGLETON_STATUS\tDESCRIPTION");
    singleP.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

    p.print("H\tUnique Peptide only : ");
    p.println(isUniquePeptide ? "true" : "false");

    //p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
    p.println("H\tPLINE\tProtein line");
    p.println("H\tSLINE\tPeptide line");
    p.println("H\t&SLINE\tSingleton Peptide line");
    if (isN15param == true) {
      p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tCOMPOSITE_RATIO\tCOMPOSITE_RATIO_STANDARD_DEVIATION\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tTOTAL_PEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tFRACTIONAL_ABUND\tDESCRIPTION");
    } else {
      p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tCOMPOSITE_RATIO\tCOMPOSITE_RATIO_STANDARD_DEVIATION\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tTOTAL_PEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tDESCRIPTION");
    }
    //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
    p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREV_SLOPE_RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tPROBABILITY_SCORE\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
    p.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tPEAK_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

    double eachSeg = (double) 100 / proteinList.size();
    double percent = 0;

    StringBuffer proteinSb = new StringBuffer();
    StringBuffer sproteinSb = new StringBuffer();
    boolean isSameGroup = false;

    List<ChroProtein> proteinResultList = new ArrayList<ChroProtein>();
    List<ChroProtein> proteinANList = new ArrayList<ChroProtein>();
    Hashtable<String, TDoubleArrayList> highScoreProteinHt = new Hashtable<String, TDoubleArrayList>();



    ///////***********   End of Analyzing singleton peptides   ******************************/
    for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext(); ) {
      HashSet<String> discardPeptideSet = new HashSet<>();
      ChroProtein protein = proItr.next();

      if (!noFilter && discardReverseProtein && (protein.getLocus().startsWith("Rev"))) // || protein.getLocus().startsWith("cont")))
      {
        continue;
      }
      redunProteinCount++;

      if (!isSameGroup) {
        proteinSb = new StringBuffer();
        sproteinSb = new StringBuffer();
        proteinGroupCount = 0;
      }

      proteinGroupCount++;

      if (protein.isRedundant()) {
        isSameGroup = true;
        proteinSb.append("P\t");
        proteinSb.append(protein.getLocus());
        proteinSb.append("\n");
        sproteinSb.append("P\t");
        sproteinSb.append(protein.getLocus());
        sproteinSb.append("\n");

        continue;
      } else {
        isSameGroup = false;
      }

      List<ChroPeptide> peptideList = protein.getPeptideList();
      List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

      DescriptiveStatistics statIntensityValues = new DescriptiveStatistics();      //find intensity threshold outlier
      DescriptiveStatistics peakAreaRatioList = new DescriptiveStatistics();
      DescriptiveStatistics anPeakAreaRatioList = new DescriptiveStatistics();

      proteinSb.append("P\t");
      proteinSb.append(protein.getLocus());
      proteinSb.append("\t");

      StringBuffer pepSb = new StringBuffer();

      double proteinIntensity=0;
      for (int j = 0; j < peptideList.size(); j++) {

        totalCount++;

//		    ChroPeptide peptide = pepItr.next();
        ChroPeptide peptide = peptideList.get(j);
//                    System.out.println(peptide.getScanNum() + "=========" + peptide.getAreaRatio() + " " + (peptide.getMaxIntensity() < intensityThresholdComposite) + " " + peptide.getMaxIntensity() + " " + INTENSITY_COMPOSITE_RATIO + " " + intensityThresholdComposite);
        //if(param.getIntensityThreshold()>0 && (peptide.getMaxIntensity() < intensityThresholdComposite))



     //   List<ChroData> l = peptide.getDataList();


        double[] regArr = peptide.getDiaFragRegressionArr();
        List<GaussianPeakModel> gaussinList = peptide.getGaussianPeakModelList();


        double totalPepIntensity = 0;
        double areaRatioSum = 0;
        double regSum = 0;

        List<GaussianPeakModel> gaussianPeakModelList = peptide.getGaussianPeakModelList();
        for(int k=0;k<regArr.length;k++) {
          regSum += regArr[k];
          if(regArr[k]>0.9) {
            GaussianPeakModel gmodel = gaussianPeakModelList.get(k);
            if(null == gmodel)
              continue;

            totalPepIntensity += gmodel.getPeakArea();
          }

        }

        pepSb.append("S\t");
        pepSb.append(peptide.isUnique() ? "U" : "");
        pepSb.append("\t");
        pepSb.append(peptide.getSequence());
        pepSb.append("\t");
        pepSb.append(totalPepIntensity);
        pepSb.append("\t");
        pepSb.append(regSum);
        pepSb.append("\t");
        pepSb.append(regSum/peptide.getSequence().length());



        pepSb.append("\t");
        pepSb.append(peptide.getFileName()).append("\t");
        pepSb.append(peptide.getScanNum()).append("\t");
        pepSb.append(peptide.getChargeState());

        pepSb.append("\n");

        proteinIntensity += totalPepIntensity;

      }


      proteinSb.append(protein.getPepCount());
      proteinSb.append("\t");
      proteinSb.append(protein.getSpectrumCount());
      proteinSb.append("\t");
      proteinSb.append(proteinIntensity);
      proteinSb.append("\t");


      proteinSb.append(protein.getDescription());
      proteinSb.append("\n");

      result.append(proteinSb.toString());
      result.append(pepSb.toString());

      uniqueProteinCount++;
    }


    rResult.setResult(result);
    rResult.setSingletonResult(singletonResult);
    rResult.setTotalCount(totalCount);
    rResult.setQuantifiedCount(quantifiedCount);
    rResult.setQuantifiedCountWithSingleton(quantifiedCount + totalSingletonCount);
    rResult.setRedunProteinCount(redunProteinCount);
    rResult.setUniqueProteinCount(uniqueProteinCount);
    rResult.setProteinGroupCount(proteinGroupCount);

    p.print("H\t");
    p.print("Total Redundant Proteins\t");
    p.println(redunProteinCount);
    p.print("H\t");
    p.print("Total Unique Proteins\t");
    p.println(uniqueProteinCount);
    p.print("H\t");
    p.print("Total peptides\t");
    p.println(totalCount);
    p.print("H\t");
    p.print("Quantified peptides\t");
    p.print(quantifiedCount);
    p.print("\n");
    p.print("H\t");
    p.print("Quantification efficiency\t");
    p.print(CensusHelper.format.format((double) quantifiedCount / totalCount * 100));
    p.print(" %\n");
    p.print("H\t");
    p.print("Quantification efficiency including singletons\t");

    p.print(CensusHelper.format.format(   ( (double)(totalSingletonCount+quantifiedCount) / totalCount) * 100));
    p.print(" %\n");
    p.print("H\t");
    p.print("Correction Factor (Ln)\t");
    p.println(correctFactorValue);


    p.print(result.toString());

    return rResult;
  }


  public static ReportResult runReportN15(ReportParam param,
            PrintStream p,
            PrintStream singleP
    ) throws Exception {

        ReportResult rResult = new ReportResult();
        StringBuffer result = new StringBuffer();
        StringBuffer singletonResult = new StringBuffer();
        int totalCount = 0;
        int quantifiedCount = 0;
        int redunProteinCount = 0;
        int uniqueProteinCount = 0;
        int proteinGroupCount = 0;

      param.setUseProfileScore(true);

        ArrayList<ChroProtein> proteinList = param.getProteinList();
        Configuration conf = param.getConf();
        boolean discardAN = param.isDiscardAN();
        boolean noFilter = param.isNoFilter();

        boolean isGui = param.isIsGui();
        JProgressBar aJProgressBar = param.getAJProgressBar();
        double detValue = param.getDetValue();
        boolean filterFragmentIons = param.isFilterFragmentIons();
        double correctFactorValue = param.getCorrectFactorValue();
        boolean discardUnlabeledPeptide = param.isDiscardUnlabeledPeptide();
        boolean discardReverseProtein = param.isDiscardReverseProtein();
        boolean removeNegative = param.isRemoveNegative();
        boolean detSelect = param.isDetSelect();
        boolean isUniquePeptide = param.isIsUniquePeptide();
        boolean pValueSelect = param.isPValueSelect();
        double pValue = param.getPValue();
        double allNoneLowerBound = param.getAllNoneLowerBound();
        double allNoneUpperBound = param.getAllNoneUpperBound();
        double allNoneCompositeScore = param.getAllNoneCompositeScore();
        int allNoneMinPeptide = param.getAllNoneMinPeptide();
        double profileScore = param.getProfileScore();
        int totalSingletonCount=0;

        int maxSpectrumShift = param.getMaxSpectrumShift();

        /*
         System.out.println(discardAN);
         System.out.println(noFilter);
         System.out.println(detValue);
         System.out.println(filterFragmentIons);
         System.out.println(correctFactorValue);
         System.out.println(discardUnlabeledPeptide);
         System.out.println(removeNegative);
         System.out.println(detSelect);
         System.out.println(isUniquePeptide);
         System.out.println(pValueSelect);
         System.out.println(pValue);
         System.out.println(allNoneLowerBound);
         System.out.println(allNoneUpperBound);
         System.out.println(allNoneCompositeScore);
         System.out.println(allNoneMinPeptide);
         */
        printHeader(p);
        printHeader(singleP);

      p.print("H\tNo filter : ");


      if(noFilter) {
          p.println("true");
        } else {
          p.println("false");

        }

        if (detSelect) {
            p.print("H\tDeterminant Factor : ");
            p.println(detValue);
        } else {
            p.println("H\tNo Determinant Factor");
        }

        if (pValueSelect) {
            p.print("H\tIterate Outlier: ");
            p.println(param.isIterateOutlier());
            p.print("H\tOutlier pValue: ");
            p.println(pValue);
        } else {
            p.println("H\tNo Outlier pValue");
        }

        if (filterFragmentIons) {
            p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
        }

        p.print("H\tDiscard reverse proteins : ");
        p.println(discardReverseProtein);
        p.print("H\tSmoothing peaks: ");
        p.println(param.isSmoothingPeaks());

        p.print("H\tCorrection Factor Value : ");
        p.println(correctFactorValue);
        p.print("H\tallNoneLowerBound : ");
        p.println(allNoneLowerBound);
        p.print("H\tallNoneUpperBound : ");
        p.println(allNoneUpperBound);
        p.print("H\tallNoneCompositeScore : ");
        p.println(allNoneCompositeScore);
        p.print("H\tallNoneMinPeptideNum: ");
        p.println(allNoneMinPeptide);
        p.print("H\tprofileScore : ");
        p.println(profileScore);
        p.print("H\tmaxScanShift : ");
        p.println(maxSpectrumShift);

        singleP.print("H\tallNoneLowerBound : ");
        singleP.println(allNoneLowerBound);
        singleP.print("H\tallNoneUpperBound : ");
        singleP.println(allNoneUpperBound);
        singleP.print("H\tallNoneCompositeScore : ");
        singleP.println(allNoneCompositeScore);
        singleP.print("H\tallNoneMinPeptideNum: ");
        singleP.println(allNoneMinPeptide);
        singleP.print("H\tmaxScanShift : ");
        singleP.println(maxSpectrumShift);
        singleP.println("H\tPLINE\tProtein line");
        singleP.println("H\t&SLINE\tSingleton Peptide line");
        singleP.println("H\tPLINE\tLOCUS\tPEPTIDE_NUM\tSPEC_COUNT\tLIGHT_SPEC_COUNT\tHEAVY_SPEC_COUNT\tSINGLETON_STATUS\tDESCRIPTION");
        singleP.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

        p.print("H\tUnique Peptide only : ");
        p.println(isUniquePeptide ? "true" : "false");


        p.print("H\tMinimum peptides per protein : ");
        p.println(param.getMinimumPeptidePerProtein());

        p.print("H\tN15 Source: ");
        p.println(param.getN15Source());
        p.print("H\tN15 APE Threshold: ");
        p.println(param.getN15ApeThreshold());
        p.print("H\tN15 Profile Threshold: ");
        p.println(param.getN15ProfileThreshold());
        p.print("H\tN15 Outlier Threshold: ");
        p.println(param.getN15OutlierThreshold());

        //p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
        p.println("H\tPLINE\tProtein line");
        p.println("H\tSLINE\tPeptide line");
        p.println("H\t&SLINE\tSingleton Peptide line");
        p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tENRICHMENT\tDESCRIPTION");
	    //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
        //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tBEST_ENRICH_CORR\tBEST_ENRICH_DELCN\tCORR_ONE_PLUS\tCORR_ONE_MINUS\tENRICHMENT\tENRICHMENT_MR");
        p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREV_SLOPE_RATIO\tDETERMINANT_FACTOR\tREGRESSION_FACTOR\tPROBABILITY_SCORE\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tBEST_ENRICH_CORR\tBEST_ENRICH_DELCN\tCORR_ONE_PLUS\tCORR_ONE_MINUS\tENRICHMENT\tENRICHMENT_MR");
        p.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

        double eachSeg = (double) 100 / proteinList.size();
        double percent = 0;

        StringBuffer proteinSb = new StringBuffer();
        StringBuffer sproteinSb = new StringBuffer();
        boolean isSameGroup = false;

        List<ChroProtein> proteinResultList = new ArrayList<ChroProtein>();
        List<ChroProtein> proteinANList = new ArrayList<ChroProtein>();
        Hashtable<String, TDoubleArrayList> highScoreProteinHt = new Hashtable<String, TDoubleArrayList>();

        if (!discardAN) {
            for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext();) {
                ChroProtein protein = proItr.next();

                List<ChroPeptide> peptideList = protein.getPeptideList();

                if (discardReverseProtein && protein.getLocus().startsWith("Rev")) // || protein.getLocus().startsWith("cont"))
                {


                    continue;
                }

                ChroProtein newChroProtein = new ChroProtein();
                newChroProtein.setLocus(protein.getLocus());
                newChroProtein.setDescription(protein.getDescription());

                for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
                    ChroPeptide peptide = pepItr.next();
                    newChroProtein.addPeptide(peptide);

                }

                proteinResultList.add(newChroProtein);
            }
        }

        Hashtable<String, ChroProtein> anProteinHt = new Hashtable<String, ChroProtein>();

        //************************************  analyze singleton peptides  *********************************************/
        for (Iterator<ChroProtein> itr = proteinResultList.iterator(); itr.hasNext();) {

            ChroProtein proAN = itr.next();
            //String desc = proAN.getDescription();
            if (discardReverseProtein && proAN.getLocus().startsWith("Rev")) // || proAN.getLocus().startsWith("cont"))
            {
                continue;
            }

            anProteinHt.put(proAN.getLocus(), proAN);

            for (Iterator<ChroPeptide> itrp = proAN.getPeptideList().iterator(); itrp.hasNext();) {
                ChroPeptide pepAN = itrp.next();

                int peakStart = Integer.parseInt(pepAN.getStartRange());
                int peakEnd = Integer.parseInt(pepAN.getEndRange());

                List dataList = pepAN.getDataList();

                AllNoneUtil.getANScore(pepAN, dataList, peakStart, peakEnd);
            }
        }

	    ///////***********   End of Analyzing singleton peptides   ******************************/
        for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext();) {
            ChroProtein protein = proItr.next();

            if (!noFilter && discardReverseProtein && (protein.getLocus().startsWith("Rev"))) // || protein.getLocus().startsWith("cont")))
            {
                continue;
            }
            redunProteinCount++;

            if (!isSameGroup) {
                proteinSb = new StringBuffer();
                sproteinSb = new StringBuffer();
                proteinGroupCount = 0;
            }

            proteinGroupCount++;

            if (protein.isRedundant()) {
                isSameGroup = true;
                proteinSb.append("P\t");
                proteinSb.append(protein.getLocus());
                proteinSb.append("\n");
                sproteinSb.append("P\t");
                sproteinSb.append(protein.getLocus());
                sproteinSb.append("\n");

                continue;
            } else {
                isSameGroup = false;
            }

            List<ChroPeptide> peptideList = protein.getPeptideList();
            List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();
            Set<String> sequenceSet = new HashSet<>();
            for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
                ChroPeptide peptide = pepItr.next();
                peptide.setFilterOut(false);
                sequenceSet.add(peptide.getSequence());

                totalCount++;
                List<ChroData> l = peptide.getDataList();

                int startRange = Integer.parseInt(peptide.getStartRange());
                int endRange = Integer.parseInt(peptide.getEndRange());
                AllNoneUtil.getANScore(peptide, l, startRange, endRange);

                long[] samArr = new long[l.size()];
                long[] refArr = new long[samArr.length];

                double samIntSum = 0;
                double refIntSum = 0;

                int index = 0;
                int startIndex = 0;
                int endIndex = samArr.length - 1;

                LinearRegression reg = null;
                LinearRegression regrev = null;

                for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext();) {
                    ChroData data = dataItr.next();

                    samArr[index] = data.getSampleIntensity();
                    refArr[index] = data.getRefIntensity();

                    int scanTemp = data.getScanNum();
                    if (startRange >= scanTemp) {
                        startIndex = index;
                    }
                    if (endRange >= scanTemp) {
                        endIndex = index;
                    }

                    index++;
                }

                samIntSum = 0;
                refIntSum = 0;

                endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);
                for (int ii = startIndex; ii <= endIndex; ii++) {
                    samIntSum += samArr[ii];
                    refIntSum += refArr[ii];

                }

                if (param.isSmoothingPeaks()) {
                    reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
                    regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
                } else {
                    reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift());
                    regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift());
                }

                peptide.setSpectraDataPoints(endIndex - startIndex + 1);

                double slope = reg.getSlope();
                double intercept = reg.getIntercept();
                double slopeRev = regrev.getSlope();
                double interceptRev = regrev.getIntercept();


              double logSlope =  slope > 0 ? (Math.log(slope) + correctFactorValue) : 0;
              double logSlopeRev =  slopeRev > 0 ? (Math.log(slopeRev) - correctFactorValue): 0;
              slope = Math.exp( logSlope);
              slopeRev = Math.exp(logSlopeRev);
                peptide.setSlope(slope);
                peptide.setCorr(reg.getCorr());
                peptide.setSlopeRev(slopeRev);
                peptide.setCorrRev(regrev.getCorr());

                peptide.setSamIntensity(samIntSum);
                peptide.setRefIntensity(refIntSum);


                //if (param.isUseProfileScore() && peptide.getAnCompositeScore() < param.getProfileScore()) {
              if (param.isUseProfileScore() && peptide.getAnCompositeScore() < param.getProfileScore()) {
         //       System.out.println("aa--" + " " + peptide.getAnCompositeScore() + " " + param.getProfileScore() + " " + param.isUseProfileScore());

                    continue;
                }

                //peptide.setSnRatio(snRatio);
                if (!noFilter) {
                    if (discardUnlabeledPeptide) {
                        double d1 = peptide.getLightMass();
                        double d2 = peptide.getHeavyMass();

                        if (d1 > 0 && d2 > 0 && (d1 == d2)) {
                            continue;
                        }
                    }

                    if (reg.getCorr() < 0 && removeNegative) {
                        continue;
                    }

                    if (detSelect && detValue > reg.getCorr() * reg.getCorr()) {
                        continue;
                    }

                    if (isUniquePeptide && !peptide.isUnique()) {
                        continue;
                    }
                }

                tempPepList.add(peptide);

            }

            if (!noFilter && tempPepList.size() >= 3 && pValueSelect) {

                if (param.isIterateOutlier()) {
				//keep outlier iteration now
                    //iterate until there is no outliers
                    int currentSize;
                    int prevSize = tempPepList.size();
                    while (true) {
					//dArr = edu.scripps.pms.stats.GrubbsTest.filter(tmpArr, 0.1);
                        //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);
                        //edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);
                        edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, param.getN15OutlierThreshold());

                        currentSize = tempPepList.size();

                        if (prevSize <= currentSize) {
                            break;
                        }

                        prevSize = currentSize;
                    }

                } else {
                    edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);
                    //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);

                }
            }

            int peptideCount = 0;
            double averageRatio = 0;
            double ratioSum = 0;
            double logRatioSum = 0;
            double logAverageRatio = 0;

            double averageRatioRev = 0;
            double ratioSumRev = 0;
            double logRatioSumRev = 0;
            double logAverageRatioRev = 0;

            for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
                ChroPeptide each = tempItr.next();

                if (!noFilter && each.isFilterOut()) {

                    continue;
                }
                if (!noFilter && ((Double.compare(each.getSlope(), Double.NaN) == 0) || each.getSlope() == 0)) {

                    continue;
                }

                ratioSum += each.getSlope();
                ratioSumRev += each.getSlopeRev();
                double logRatio = Math.log(each.getSlope());
                double logRatioRev = Math.log(each.getSlopeRev());
//System.out.println("==" + ratioSum + "\t" + logRatio + "\t"  + each.getSlope());
                logRatioSum += logRatio;
                logRatioSumRev += logRatioRev;

//		    System.out.println("==\t" + each.getSlope() + "\t" + logRatio);
                peptideCount++;
                quantifiedCount++;
            }

            double invLogRatio = 0;
            double invLogRatioRev = 0;
            if (peptideCount > 0) {
                averageRatio = ratioSum / peptideCount;
                averageRatioRev = ratioSumRev / peptideCount;
//		    System.out.println("---->>\t" + peptideCount + "\t" + ratioSum + "\t" + averageRatio );

                logAverageRatio = logRatioSum / peptideCount;
                logAverageRatioRev = logRatioSumRev / peptideCount;
                invLogRatio = Math.exp(logAverageRatio);
                invLogRatioRev = Math.exp(logAverageRatioRev);
            }

            proteinSb.append("P\t");
            proteinSb.append(protein.getLocus());
            proteinSb.append("\t");
            proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
            proteinSb.append("\t");
            proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
            proteinSb.append("\t");
		//let's use log value
            //proteinSb.append( averageRatio>0?CensusHelper.format.format(averageRatio):"NA" );
            //proteinSb.append("\t");
            //proteinSb.append( averageRatioRev>0?CensusHelper.format.format(averageRatioRev):"NA" );
            //proteinSb.append("\t");
            sproteinSb.append("P\t");
            sproteinSb.append(protein.getLocus());
            sproteinSb.append("\t");

            double devSum = 0;
            double devSumRev = 0;

            StringBuffer pepSb = new StringBuffer();

            WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();

            //double totalPepIntensity=0;
            double areaRatioSum = 0;

            for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
                ChroPeptide each = tempItr.next();

                if (!noFilter && each.isFilterOut()) {

                    continue;
                }

                    //discard peptides with low profile score
                //if(each.getAnCompositeScore()<0.5) continue;
                //if(noFilter || !each.isFilterOut())
                double dev = each.getSlope() - averageRatio;
                double devRev = each.getSlopeRev() - averageRatioRev;

                if ((Double.compare(each.getSlope(), Double.NaN) != 0)) {
                    devSum += dev * dev;
                    devSumRev += devRev * devRev;
                }

                double corr = each.getCorr();
                int dataNum = each.getSpectraDataPoints();
                double tvalue = corr / Math.sqrt((1 - corr * corr) / (dataNum - 2));
                int df = dataNum - 1;

                if (df <= 0) {

                    continue;
                }
		//new filters
//		    if(each.getAnCompositeScore()<0.5) continue;
//		    if(each.getEnrichment()<0.8) continue;

                if (each.getAnCompositeScore() < param.getN15ProfileThreshold()) {

             //     System.out.println("6");
                    continue;
                }
                //a.p.e threshold
                if (each.getEnrichment() < param.getN15ApeThreshold()) {

                    continue;
                }

                pepSb.append("S\t");
                pepSb.append(each.isUnique() ? "U" : "");
                pepSb.append("\t");
                pepSb.append(each.getSequence());
                pepSb.append("\t");

                if (each.getCorr() < 0) {

                    if ((Double.compare(each.getSlope(), Double.NaN) == 0)) {
                        pepSb.append("NA\tNA\t");
                    } else {
                        pepSb.append(CensusHelper.format.format(each.getSlope()));
                        pepSb.append("\t");
                        pepSb.append(CensusHelper.format.format(each.getSlopeRev()));
                        pepSb.append("\t");
                    }

                    pepSb.append(CensusHelper.format.format(each.getCorr()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getDetValue()));
                } else {
                    pepSb.append(CensusHelper.format.format(each.getSlope()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getSlopeRev()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getCorr()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getDetValue()));
                }

//			System.out.println("====" + detSelect + "\t" + detValue + "\t" + each.getCorr() + "\t" +  (each.getCorr()*each.getCorr()));
                   // System.out.println(CensusHelper.format.format(each.getAnCompositeScore()) +" " + each.getFileName());
                TDistribution t = new TDistribution(df);
                double d = 1 - t.cumulativeProbability(tvalue);

                pepSb.append("\t");
                pepSb.append(d);
                pepSb.append("\t");
                pepSb.append((null == each.getXCorr()) ? "" : each.getXCorr());
                pepSb.append("\t");
                pepSb.append((null == each.getDeltCN()) ? "" : each.getDeltCN());
                pepSb.append("\t");
                pepSb.append(each.getSamIntensity());
                pepSb.append("\t");
                pepSb.append(each.getRefIntensity());
                pepSb.append("\t");
                double intRatio = (0 == each.getRefIntensity()) ? -1 : (each.getSamIntensity() / each.getRefIntensity());
                pepSb.append((intRatio >= 0) ? CensusHelper.format.format(intRatio) : "INF");
                pepSb.append("\t");

                //15N remaining
                double enrichRatio = each.getRefIntensity() / (each.getSamIntensity() + each.getRefIntensity());
		    //       pepSb.append( CensusHelper.format.format(each.getSnRatio()) );
                //       pepSb.append("\t");
                pepSb.append(CensusHelper.format.format(each.getAnCompositeScore()));

                pepSb.append("\t");
                pepSb.append(each.getFileName()).append("\t");
                pepSb.append(each.getScanNum()).append("\t");
                pepSb.append(each.getChargeState());

                if (each.getBestEnrichCorr() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getBestEnrichCorr());
                }

                if (each.getBestEnrichDelCN() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getBestEnrichDelCN());
                }

                if (each.getCorrOnePlus() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getCorrOnePlus());
                }

                if (each.getCorrOneMinus() <= -1) {
                    pepSb.append("\t").append("N/A");
                } else {
                    pepSb.append("\t").append(each.getCorrOneMinus());
                }

                double enrichment = each.getEnrichment();

                //	    if(enrichment>0)
                pepSb.append("\t").append(enrichment);
	//	    else
                //		pepSb.append("\t").append("N/A");

//robin
//		    //multiple regression analysis for enrichment
// constant: 0.2606391001547484
// ratio: -0.0013865804600325519
// profile: -0.1534407004734698
// correlation: 0.20870496987846332
// delCN: -0.05808530803138859
// plus: 0.16118244323905673
// minus: -0.005428458538558106
// enrich: 0.5221466417403269
                /*
                 if(each.getCorrOneMinus()<=0 && each.getCorrOnePlus()<=0) {
                 enrichment = enrichRatio;
                 } else {
                 //check peak area
                 double diff = Math.abs(enrichment-enrichRatio); //isoEnrich);
                 if(diff>0.1)
                 enrichment = -1;

                 }
                 */
                   //continue
                //  System.out.println(each.getCorr() + " " +each.getSequence() + "\t" + enrichment + " " + each.getEnrichment() + " " + enrichRatio);
                //
                pepSb.append("\t").append(enrichRatio);
                pModel.addEnrichRatio(enrichRatio);
                pepSb.append("\n");
                /*
                 if(each.getCorr()>0.75) {
                 if(enrichment>=0)
                 pepSb.append("\t").append(enrichment);
                 else
                 pepSb.append("\t").append("N/A");
                 } else {

                 if(each.getSamIntensity()>each.getRefIntensity()*4) {//System.out.print("1");
                 pepSb.append("\t").append("0.0"); }
                 else if(each.getSamIntensity()*4<each.getRefIntensity()) {//System.out.print("2");
                 pepSb.append("\t").append("1.0"); }
                 else {//System.out.print("3");
                 pepSb.append("\t").append("N/A");}

                 }
                 */

                double rsqrtLog = Math.log(each.getDetValue());
                double stdevLog = -0.84 * rsqrtLog + 0.43;
                double invStdev = Math.exp(stdevLog);

                //grap only valid ratio
                if (each.getSlope() > 0) {
                    pModel.add(invStdev, each.getSlope());
                    // double dev = each.getSlope()-averageRatio;
                }

                //Hashtable<String, ArrayList> highScoreProteinHt = new Hashtable<String, ArrayList>();
                TDoubleArrayList al = highScoreProteinHt.get(protein.getLocus());

                if (null == al) {
                    al = new TDoubleArrayList();
                    al.add(each.getAreaRatio());
                    highScoreProteinHt.put(protein.getLocus(), al);
                } else {
                    al.add(each.getAreaRatio());
                }

		//    totalPepIntensity += each.getSamIntensity();
                //    totalPepIntensity += each.getRefIntensity();
                areaRatioSum += each.getAreaRatio();
            }

            ////////////add all none peptides here
            ChroProtein anProtein = anProteinHt.get(protein.getLocus());
            StringBuffer singlePepSb = new StringBuffer();
		//System.out.println(anProtein.getLocus());
            //anProtein.getPeptideList()
            int singletonPeptideCount = 0;
            int singletonUpCount = 0;
            int singletonDownCount = 0;

            if (!discardAN && null != anProtein) {
                for (Iterator<ChroPeptide> anPItr = anProtein.getPeptideList().iterator(); anPItr.hasNext();) {
                    ChroPeptide anPep = anPItr.next();
                    /*
                     if(anPep.getSequence().equals("R.VLLVIDEPHTDWAK.Y"))
                     {
                     }*/

                    if ((anPep.getAreaRatio() < allNoneUpperBound && anPep.getAreaRatio() > allNoneLowerBound)) {
                        continue;
                    }

                    if (detSelect && detValue < anPep.getDetValue()) {
                        continue;
                    }

        //            System.out.println("aaaaaaaaaaaaaaa");
                    if (anPep.getAnCompositeScore() >= allNoneCompositeScore) {

                      totalSingletonCount++;
                        StringBuffer tmpSb = new StringBuffer();
                        tmpSb.append("&S").append("\t").append(anPep.isUnique() ? "U" : "").append("\t").append(anPep.getSequence()).append("\t");

                        if (anPep.getCorr() < 0) {
                            tmpSb.append("0.0");
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getCorr()));
                            tmpSb.append("\t");
                            tmpSb.append("0.0");
                        } else {
                            tmpSb.append(CensusHelper.format.format(anPep.getSlope()));
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getCorr()));
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getDetValue()));
                        }

                        tmpSb.append("\t");
                        tmpSb.append((null == anPep.getXCorr()) ? "" : anPep.getXCorr());
                        tmpSb.append("\t");
                        tmpSb.append((null == anPep.getDeltCN()) ? "" : anPep.getDeltCN());
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getSamIntensity());
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getRefIntensity());
                        tmpSb.append("\t");
                        double intRatio = (0 == anPep.getRefIntensity()) ? -1 : (anPep.getSamIntensity() / anPep.getRefIntensity());

                        if (anPep.getSamIntensity() >= anPep.getRefIntensity()) {
                            singletonUpCount++;
                        } else {
                            singletonDownCount++;
                        }

                        tmpSb.append((intRatio >= 0) ? CensusHelper.format.format(intRatio) : "INF");
                        //                                        tmpSb.append( CensusHelper.format.format(anPep.getSamIntensity()/anPep.getRefIntensity()) );
                        tmpSb.append("\t");
                        tmpSb.append(CensusHelper.format.format(anPep.getAnCompositeScore()));
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getFileName()).append("\t");
                        tmpSb.append(anPep.getScanNum()).append("\t");
                        tmpSb.append(anPep.getChargeState());

			    //pepSb.append(tmpSb);
                        tmpSb.append(anPep.getChargeState());

                        if (anPep.getBestEnrichCorr() <= -1) {
                            tmpSb.append("\t").append("N/A");
                        } else {
                            tmpSb.append("\t").append(anPep.getBestEnrichCorr());
                        }

                        if (anPep.getBestEnrichDelCN() <= -1) {
                            tmpSb.append("\t").append("N/A");
                        } else {
                            tmpSb.append("\t").append(anPep.getBestEnrichDelCN());
                        }

                        if (anPep.getCorrOnePlus() <= -1) {
                            tmpSb.append("\t").append("N/A");
                        } else {
                            tmpSb.append("\t").append(anPep.getCorrOnePlus());
                        }

                        if (anPep.getCorrOneMinus() <= -1) {
                            tmpSb.append("\t").append("N/A");
                        } else {
                            tmpSb.append("\t").append(anPep.getCorrOneMinus());
                        }

                        double enrichment = anPep.getEnrichment();

//System.out.print("\t" + enrichment);
                        if (enrichment > 0) {
                            tmpSb.append("\t").append(enrichment);
                        } else {
                            tmpSb.append("\t").append("N/A");
                        }

                        tmpSb.append("\n");

                        singletonPeptideCount++;
                        singlePepSb.append(tmpSb);
                    }
                }
            }


      //      System.out.println(peptideCount + "==========================\t" + param.getMinimumPeptidePerProtein());
      //    if(peptideCount<param.getMinimumPeptidePerProtein()) continue;
          if(sequenceSet.size()<param.getMinimumPeptidePerProtein()) continue;

          //////////End of adding an peptides
            if (peptideCount > 1 && Double.compare(devSum, Double.NaN) != 0) {
                proteinSb.append(CensusHelper.format.format(Math.sqrt(devSum / (peptideCount - 1))));
                proteinSb.append("\t");
                proteinSb.append(CensusHelper.format.format(Math.sqrt(devSumRev / (peptideCount - 1))));
            } else {
                proteinSb.append("NA\tNA");
            }

            proteinSb.append("\t");
              //  System.out.println("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
            // System.out.println("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
            //       + pModel.getAverageEnrichRatio(true, 0.1));
            if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
                proteinSb.append("NA\t");

            } else {
                proteinSb.append(CensusHelper.format.format(pModel.getStandardWeightedAverage()));
                proteinSb.append("\t");
            }

            proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
            proteinSb.append("\t");
            proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
            proteinSb.append("\t");
            proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
            proteinSb.append("\t");
            proteinSb.append(protein.getSpectrumCount());
            proteinSb.append("\t");
            proteinSb.append(protein.getLspectrumCount());
            proteinSb.append("\t");
            proteinSb.append(protein.getHspectrumCount());
            proteinSb.append("\t");
            //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
            proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
            proteinSb.append("\t");
            proteinSb.append(pModel.getAverageEnrichRatio(true, 0.1));
            proteinSb.append("\t");
            proteinSb.append(protein.getDescription());
            proteinSb.append("\n");

            sproteinSb.append(singletonPeptideCount);
            sproteinSb.append("\t");
            sproteinSb.append(protein.getSpectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(protein.getLspectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(protein.getHspectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
            sproteinSb.append("\t");
            sproteinSb.append(protein.getDescription());
            sproteinSb.append("\n");

            //here
            List redProList = protein.getRedunList();
            for (Iterator<ChroProtein> redPItr = redProList.iterator(); redPItr.hasNext();) {
                ChroProtein redPro = redPItr.next();

                proteinSb.append("P\t");
                proteinSb.append(redPro.getLocus());
                proteinSb.append("\t");
                proteinSb.append(averageRatio > 0 ? CensusHelper.format.format(averageRatio) : "NA");
                proteinSb.append("\t");
                proteinSb.append(averageRatioRev > 0 ? CensusHelper.format.format(averageRatioRev) : "NA");
                proteinSb.append("\t");

                if (peptideCount > 1) {
                    proteinSb.append(CensusHelper.format.format(Math.sqrt(devSum / (peptideCount - 1))));
                    proteinSb.append("\t");
                    proteinSb.append(CensusHelper.format.format(Math.sqrt(devSumRev / (peptideCount - 1))));
                } else {
                    proteinSb.append("NA\tNA");
                }

                proteinSb.append("\t");

                if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
                    proteinSb.append("NA\t");

                } else {
                    proteinSb.append(CensusHelper.format.format(pModel.getStandardWeightedAverage()));
                    proteinSb.append("\t");
                }
                //proteinSb.append( CensusHelper.format.format(pModel.getStandardWeightedAverage()) );
                proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
                proteinSb.append("\t");
                proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
                proteinSb.append("\t");
                proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
                proteinSb.append("\t");
                proteinSb.append(redPro.getSpectrumCount());
                proteinSb.append("\t");
                proteinSb.append(redPro.getLspectrumCount());
                proteinSb.append("\t");
                proteinSb.append(redPro.getHspectrumCount());
                proteinSb.append("\t");
                //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
                proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
                proteinSb.append("\t");
                proteinSb.append(pModel.getAverageEnrichRatio(true, 0.1));
                proteinSb.append("\t");
                proteinSb.append(redPro.getDescription());
                proteinSb.append("\n");

                sproteinSb.append("P\t");
                sproteinSb.append(redPro.getLocus());
                sproteinSb.append("\t");
                sproteinSb.append(singletonPeptideCount);
                sproteinSb.append("\t");
                sproteinSb.append(protein.getSpectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(protein.getLspectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(protein.getHspectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
                sproteinSb.append("\t");
                sproteinSb.append(protein.getDescription());
                sproteinSb.append("\n");

            }

            percent += eachSeg;

            if (isGui && null != aJProgressBar) {
                aJProgressBar.setValue((int) percent);
            }

            if (singletonPeptideCount < allNoneMinPeptide) {
                singlePepSb = new StringBuffer();
            }

            if (pepSb.length() <= 0 && singlePepSb.length() <= 0) {
                redunProteinCount -= proteinGroupCount;
                continue;
            }

		//if( pepSb.length()<=0 && singlePepSb.length()<=0) {
            //if( pepSb.length()<=0 && singlePepSb.length()>0) {
            if (singlePepSb.length() > 0) {
                singletonResult.append(sproteinSb.toString());
                singletonResult.append(singlePepSb.toString());
                pepSb.append(singlePepSb.toString());
            }

            result.append(proteinSb.toString());
            result.append(pepSb.toString());

            uniqueProteinCount++;
        }

        rResult.setResult(result);
        rResult.setSingletonResult(singletonResult);
        rResult.setTotalCount(totalCount);
        rResult.setQuantifiedCount(quantifiedCount);
        rResult.setRedunProteinCount(redunProteinCount);
        rResult.setUniqueProteinCount(uniqueProteinCount);
        rResult.setProteinGroupCount(proteinGroupCount);
        rResult.setQuantifiedCountWithSingleton(quantifiedCount);


      p.print("H\t");
        p.print("Total Redundant Proteins\t");
        p.println(redunProteinCount);
        p.print("H\t");
        p.print("Total Unique Proteins\t");
        p.println(uniqueProteinCount);
        p.print("H\t");
        p.print("Total peptides\t");
        p.println(totalCount);
        p.print("H\t");
        p.print("Quantified peptides\t");
        p.print(quantifiedCount);
        p.print("\n");
        p.print("Quantified peptides including singletons\t");
        p.print(quantifiedCount);
        p.print("\n");

        p.print("H\t");
        p.print("Quantification efficiency\t");
        p.print(CensusHelper.format.format((double) quantifiedCount / totalCount * 100));
        p.print(" %\n");

        p.print("H\t");
        p.print("Correction Factor (Ln)\t");
        p.println(correctFactorValue);

        singleP.print("H\t");
        singleP.print("Total Redundant Proteins\t");
        singleP.println(redunProteinCount);
        singleP.print("H\t");
        singleP.print("Total Unique Proteins\t");
        singleP.println(uniqueProteinCount);
        singleP.print("H\t");
        singleP.print("Total peptides\t");
        singleP.println(totalCount);

        p.print(result.toString());
        singleP.print(singletonResult.toString());

        return rResult;

    }

    // added Harshi Shah
    public static void generateRegression(ReportParam param) {

        List<ChroProtein> proteinList = param.getProteinList();
        StringBuffer sb = new StringBuffer();
        sb.append("SlopeConfidenceInterval\tSlope\tR\tN\tSlopeStdErr\tgetSumSquaredErrors()\n");
        try {
            bw.write(sb.toString());
            sb = new StringBuffer();
        } catch (IOException ex) {
            Logger.getLogger(RelExMainFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
        for (ChroProtein protein : proteinList) {
            List<ChroPeptide> peptideList = protein.getPeptideList();
            for (ChroPeptide peptide : peptideList) {
                List<IsoData> isoDataList = peptide.getIsoDataList();
                List<IsoData> isoOrigDataList = peptide.getIsoOrigDataList();

                List lightList = new ArrayList();
                List heavyList = new ArrayList();

//                sb.append(protein.getLocus()).append("\t").append(peptide.getSequence() + "\n");
                for (int i = 0; i < isoDataList.get(0).getLightData().length; i++) {
                    System.out.print("LN\t");
                    double[] light = new double[isoDataList.size()];
                    for (int j = 0; j < isoDataList.size(); j++) {
                        light[j] = isoDataList.get(j).getLightData()[i];
                        System.out.print(light[j] + "\t");  //+ isoOrigDataList.get(j).getLightData()[i] + "\t");
                    }
                    System.out.println("");
                    lightList.add(light);
                }

                for (int i = 0; i < isoDataList.get(0).getHeavyData().length; i++) {
                    System.out.print("HN\t");
                    double[] heavy = new double[isoDataList.size()];
                    for (int j = 0; j < isoDataList.size(); j++) {
                        heavy[j] = isoDataList.get(j).getHeavyData()[i];
                        System.out.print(heavy[j] + "\t"); //  + isoOrigDataList.get(j).getHeavyData()[i] + "\t");
                    }

                    System.out.println("");
                    heavyList.add(heavy);
                }

                for (int i = 0; i < isoDataList.get(0).getLightData().length; i++) {
                    System.out.print("L\t");
                    double[] light = new double[isoDataList.size()];
                    for (int j = 0; j < isoDataList.size(); j++) {
                        light[j] = isoDataList.get(j).getLightData()[i];
                        System.out.print(isoOrigDataList.get(j).getLightData()[i] + "\t");
                    }
                    System.out.println("");
                    lightList.add(light);
                }

                for (int i = 0; i < isoDataList.get(0).getHeavyData().length; i++) {
                    System.out.print("H\t");
                    double[] heavy = new double[isoDataList.size()];
                    for (int j = 0; j < isoDataList.size(); j++) {
                        heavy[j] = isoDataList.get(j).getHeavyData()[i];
                        System.out.print(isoOrigDataList.get(j).getHeavyData()[i] + "\t");
                    }

                    System.out.println("");

                    heavyList.add(heavy);
                }

                runSimpleRegression(lightList, heavyList);

            }
        }

        try {
            System.out.println("regression file created");
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(RelExMainFrame.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    //added harshil Shah

    public static void write(SimpleRegression sr) {
        try {
//            StringBuffer sb =new StringBuffer();
            double val = sr.getSlopeConfidenceInterval();
            StringBuffer sb = new StringBuffer();

            sb.append(val).append("\t");
            sb.append(sr.getSlope()).append("\t");
            sb.append(sr.getR()).append("\t");
            sb.append(sr.getN()).append("\t");
            sb.append(sr.getSlopeStdErr()).append("\t");
            sb.append(sr.getSumSquaredErrors()).append("\n");

            bw.write(sb.toString());
        } catch (Exception ex) {
            Logger.getLogger(RelExMainFrame.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    //added harshil Shah
    public static void runSimpleRegression(List<double[]> lightList, List<double[]> heavyList) {
        SimpleRegression sr = new SimpleRegression();
        for (double[] currentLightData : lightList) {
            for (double[] currentHeavyData : heavyList) {
                try {
                    for (int i = 0; i < currentHeavyData.length; i++) {
                        sr.addData(currentLightData[i], currentHeavyData[i]);
                    }

                    write(sr);

                } catch (Exception ex) {
                    Logger.getLogger(RelExMainFrame.class.getName()).log(Level.SEVERE, null, ex);
                }

            }
        }
    }
//added harshil SHah

    public static ReportResult runReportIsotops(ReportParam param,
            PrintStream p,
            PrintStream singleP
    ) throws Exception {
        //     bw = new BufferedWriter(new FileWriter("e:\\regression.txt"));
        bw = new BufferedWriter(new FileWriter("/home/rpark/pms/Census/build/classes/n15/regression.txt"));
        generateRegression(param);
        ReportResult rResult = new ReportResult();
        StringBuffer result = new StringBuffer();
        StringBuffer singletonResult = new StringBuffer();
        int totalCount = 0;
        int quantifiedCount = 0;
        int redunProteinCount = 0;
        int uniqueProteinCount = 0;
        int proteinGroupCount = 0;

        ArrayList<ChroProtein> proteinList = param.getProteinList();
        Configuration conf = param.getConf();
        boolean discardAN = param.isDiscardAN();
        boolean noFilter = param.isNoFilter();

        boolean isGui = param.isIsGui();
        JProgressBar aJProgressBar = param.getAJProgressBar();
        double detValue = param.getDetValue();
        boolean filterFragmentIons = param.isFilterFragmentIons();
        double correctFactorValue = param.getCorrectFactorValue();
        boolean discardUnlabeledPeptide = param.isDiscardUnlabeledPeptide();
        boolean discardReverseProtein = param.isDiscardReverseProtein();
        boolean removeNegative = param.isRemoveNegative();
        boolean detSelect = param.isDetSelect();
        boolean isUniquePeptide = param.isIsUniquePeptide();
        boolean pValueSelect = param.isPValueSelect();
        double pValue = param.getPValue();
        double allNoneLowerBound = param.getAllNoneLowerBound();
        double allNoneUpperBound = param.getAllNoneUpperBound();
        double allNoneCompositeScore = param.getAllNoneCompositeScore();
        int allNoneMinPeptide = param.getAllNoneMinPeptide();
        double profileScore = param.getProfileScore();

        int maxSpectrumShift = param.getMaxSpectrumShift();

        printHeader(p);
        printHeader(singleP);

        if (detSelect) {
            p.print("H\tDeterminant Factor : ");
            p.println(detValue);
        } else {
            p.println("H\tNo Determinant Factor");
        }

        if (pValueSelect) {
            p.print("H\tIterate Outlier: ");
            p.println(param.isIterateOutlier());
            p.print("H\tOutlier pValue: ");
            p.println(pValue);
        } else {
            p.println("H\tNo Outlier pValue");
        }

        if (filterFragmentIons) {
            p.println("H\tFilter Fragment Ions on MS/MS pValue : true");
        }

        p.print("H\tDiscard reverse proteins : ");
        p.println(discardReverseProtein);
        p.print("H\tSmoothing peaks: ");
        p.println(param.isSmoothingPeaks());

        p.print("H\tCorrection Factor Value : ");
        p.println(correctFactorValue);
        p.print("H\tallNoneLowerBound : ");
        p.println(allNoneLowerBound);
        p.print("H\tallNoneUpperBound : ");
        p.println(allNoneUpperBound);
        p.print("H\tallNoneCompositeScore : ");
        p.println(allNoneCompositeScore);
        p.print("H\tallNoneMinPeptideNum: ");
        p.println(allNoneMinPeptide);
        p.print("H\tprofileScore : ");
        p.println(profileScore);
        p.print("H\tmaxScanShift : ");
        p.println(maxSpectrumShift);

        singleP.print("H\tallNoneLowerBound : ");
        singleP.println(allNoneLowerBound);
        singleP.print("H\tallNoneUpperBound : ");
        singleP.println(allNoneUpperBound);
        singleP.print("H\tallNoneCompositeScore : ");
        singleP.println(allNoneCompositeScore);
        singleP.print("H\tallNoneMinPeptideNum: ");
        singleP.println(allNoneMinPeptide);
        singleP.print("H\tmaxScanShift : ");
        singleP.println(maxSpectrumShift);
        singleP.println("H\tPLINE\tProtein line");
        singleP.println("H\t&SLINE\tSingleton Peptide line");
        singleP.println("H\tPLINE\tLOCUS\tPEPTIDE_NUM\tSPEC_COUNT\tLIGHT_SPEC_COUNT\tHEAVY_SPEC_COUNT\tSINGLETON_STATUS\tDESCRIPTION");
        singleP.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

        p.print("H\tUnique Peptide only : ");
        p.println(isUniquePeptide ? "true" : "false");
        p.print("H\tMinimum peptides per protein : ");
        p.println(param.getMinimumPeptidePerProtein());

        //p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tSTANDARD_DEVIATION\tPEPTIDE_NUM\tSPEC_COUNT\tDESCRIPTION");
        p.println("H\tPLINE\tProtein line");
        p.println("H\tSLINE\tPeptide line");
        p.println("H\t&SLINE\tSingleton Peptide line");
        p.println("H\tPLINE\tLOCUS\tAVERAGE_RATIO\tAVERAGE_RATIO_REV\tSTANDARD_DEVIATION\tSTANDARD_DEVIATION_REV\tWEIGHTED_AVERAGE\tLOG_INV_AVERAGE\tLOG_INV_AVERAGE_REV\tPEPTIDE_NUM\tSPEC_COUNT\tLSPEC_COUNT\tHSPEC_COUNT\tAREA_RATIO\tDESCRIPTION");
        //p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tENRICHMENT");
        p.println("H\tSLINE\tUNIQUE\tSEQUENCE\tRATIO\tREV_SLOPE_RATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tPROBABILITY_SCORE\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tPROFILE_SCORE\tFILE_NAME\tSCAN\tCS\tBEST_ENRICH_CORR\tBEST_ENRICH_DELCN\tCORR_ONE_PLUS\tCORR_ONE_MINUS\tENRICHMENT\tENRICHMENT_MR");

        p.println("H\t&SLINE\tUNIQUE\tSEQUENCE\tRATIO\tREGRESSION_FACTOR\tDETERMINANT_FACTOR\tXCorr\tdeltaCN\tSAM_INT\tREF_INT\tAREA_RATIO\tSINGLETON_SCORE\tFILE_NAME\tSCAN\tCS");

        double eachSeg = (double) 100 / proteinList.size();
        double percent = 0;

        StringBuffer proteinSb = new StringBuffer();
        StringBuffer sproteinSb = new StringBuffer();
        boolean isSameGroup = false;

        List<ChroProtein> proteinResultList = new ArrayList<ChroProtein>();
        List<ChroProtein> proteinANList = new ArrayList<ChroProtein>();
        Hashtable<String, TDoubleArrayList> highScoreProteinHt = new Hashtable<String, TDoubleArrayList>();

        if (!discardAN) {
            for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext();) {
                ChroProtein protein = proItr.next();

                List<ChroPeptide> peptideList = protein.getPeptideList();

                if (discardReverseProtein && protein.getLocus().startsWith("Rev")) // || protein.getLocus().startsWith("cont"))
                {
                    continue;
                }

                ChroProtein newChroProtein = new ChroProtein();
                newChroProtein.setLocus(protein.getLocus());
                newChroProtein.setDescription(protein.getDescription());

                for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
                    ChroPeptide peptide = pepItr.next();
                    newChroProtein.addPeptide(peptide);

                }

                proteinResultList.add(newChroProtein);
            }
        }

        Hashtable<String, ChroProtein> anProteinHt = new Hashtable<String, ChroProtein>();

        //************************************  analyze singleton peptides  *********************************************/
        for (Iterator<ChroProtein> itr = proteinResultList.iterator(); itr.hasNext();) {

            ChroProtein proAN = itr.next();
            //String desc = proAN.getDescription();
            if (discardReverseProtein && proAN.getLocus().startsWith("Rev")) // || proAN.getLocus().startsWith("cont"))
            {
                continue;
            }

            anProteinHt.put(proAN.getLocus(), proAN);

            for (Iterator<ChroPeptide> itrp = proAN.getPeptideList().iterator(); itrp.hasNext();) {
                ChroPeptide pepAN = itrp.next();

                int peakStart = Integer.parseInt(pepAN.getStartRange());
                int peakEnd = Integer.parseInt(pepAN.getEndRange());

                List dataList = pepAN.getDataList();

                AllNoneUtil.getANScore(pepAN, dataList, peakStart, peakEnd);
            }
        }

	    ///////***********   End of Analyzing singleton peptides   ******************************/
        for (Iterator<ChroProtein> proItr = proteinList.iterator(); proItr.hasNext();) {
            ChroProtein protein = proItr.next();

            if (!noFilter && discardReverseProtein && (protein.getLocus().startsWith("Rev"))) // || protein.getLocus().startsWith("cont")))
            {
                continue;
            }
            redunProteinCount++;

            if (!isSameGroup) {
                proteinSb = new StringBuffer();
                sproteinSb = new StringBuffer();
                proteinGroupCount = 0;
            }

            proteinGroupCount++;

            if (protein.isRedundant()) {
                isSameGroup = true;
                proteinSb.append("P\t");
                proteinSb.append(protein.getLocus());
                proteinSb.append("\n");
                sproteinSb.append("P\t");
                sproteinSb.append(protein.getLocus());
                sproteinSb.append("\n");

                continue;
            } else {
                isSameGroup = false;
            }

            List<ChroPeptide> peptideList = protein.getPeptideList();
            List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

            for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
                ChroPeptide peptide = pepItr.next();
                peptide.setFilterOut(false);

                totalCount++;
                List<ChroData> l = peptide.getDataList();
                List<IsoData> lIso = peptide.getIsoDataList();

                int startRange = Integer.parseInt(peptide.getStartRange());
                int endRange = Integer.parseInt(peptide.getEndRange());
                AllNoneUtil.getANScore(peptide, l, startRange, endRange);

                long[] samArr = new long[lIso.size()];
                long[] refArr = new long[samArr.length];

                double samIntSum = 0;
                double refIntSum = 0;

                int index = 0;
                int startIndex = 0;
                int endIndex = samArr.length - 1;

                LinearRegression reg = null;
                LinearRegression regrev = null;

                for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext();) {
                    ChroData data = dataItr.next();

                    samArr[index] = data.getSampleIntensity();
                    refArr[index] = data.getRefIntensity();

                    int scanTemp = data.getScanNum();
                    if (startRange >= scanTemp) {
                        startIndex = index;
                    }
                    if (endRange >= scanTemp) {
                        endIndex = index;
                    }

                    index++;
                }

                samIntSum = 0;
                refIntSum = 0;

                endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);
                for (int ii = startIndex; ii <= endIndex; ii++) {
                    samIntSum += samArr[ii];
                    refIntSum += refArr[ii];

                }

                if (param.isSmoothingPeaks()) {
                    reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
                    regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift(), true);
                } else {
                    reg = new LinearRegression(samArr, refArr, startIndex, endIndex, param.getMaxSpectrumShift());
                    regrev = new LinearRegression(refArr, samArr, startIndex, endIndex, param.getMaxSpectrumShift());
                }

                peptide.setSpectraDataPoints(endIndex - startIndex + 1);

                double slope = reg.getSlope();
                double intercept = reg.getIntercept();
                double slopeRev = regrev.getSlope();
                double interceptRev = regrev.getIntercept();
              double logSlope =  slope > 0 ? (Math.log(slope) + correctFactorValue) : 0;
              double logSlopeRev =  slopeRev > 0 ? (Math.log(slopeRev) - correctFactorValue): 0;
              slope = Math.exp( logSlope);
              slopeRev = Math.exp(logSlopeRev);
                peptide.setSlope(slope);
                peptide.setCorr(reg.getCorr());
                peptide.setSlopeRev(slopeRev);
                peptide.setCorrRev(regrev.getCorr());

                peptide.setSamIntensity(samIntSum);
                peptide.setRefIntensity(refIntSum);
                if (param.isUseProfileScore() && peptide.getAnCompositeScore() < param.getProfileScore()) {
                    continue;
                }

		    //peptide.setSnRatio(snRatio);
                if (!noFilter) {
                    if (discardUnlabeledPeptide) {
                        double d1 = peptide.getLightMass();
                        double d2 = peptide.getHeavyMass();

                        if (d1 > 0 && d2 > 0 && (d1 == d2)) {
                            continue;
                        }
                    }

                    if (reg.getCorr() < 0 && removeNegative) {
                        continue;
                    }

                    if (detSelect && detValue > reg.getCorr() * reg.getCorr()) {
                        continue;
                    }

                    if (isUniquePeptide && !peptide.isUnique()) {
                        continue;
                    }
                }

                tempPepList.add(peptide);

            }

            if (!noFilter && tempPepList.size() >= 3 && pValueSelect) {
                if (param.isIterateOutlier()) {
				//keep outlier iteration now
                    //iterate until there is no outliers
                    int currentSize;
                    int prevSize = tempPepList.size();
                    while (true) {
					//dArr = edu.scripps.pms.stats.GrubbsTest.filter(tmpArr, 0.1);
                        //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);
                        edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);

                        currentSize = tempPepList.size();

                        if (prevSize <= currentSize) {
                            break;
                        }

                        prevSize = currentSize;
                    }

                } else {
                    edu.scripps.pms.stats.GrubbsTest.filterAndRemove(tempPepList, pValue);
                    //edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);

                }
            }

            int peptideCount = 0;
            double averageRatio = 0;
            double ratioSum = 0;
            double logRatioSum = 0;
            double logAverageRatio = 0;

            double averageRatioRev = 0;
            double ratioSumRev = 0;
            double logRatioSumRev = 0;
            double logAverageRatioRev = 0;

            for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
                ChroPeptide each = tempItr.next();

                if (!noFilter && each.isFilterOut()) {
                    continue;
                }
                if (!noFilter && ((Double.compare(each.getSlope(), Double.NaN) == 0) || each.getSlope() == 0)) {
                    continue;
                }

                ratioSum += each.getSlope();
                ratioSumRev += each.getSlopeRev();
                double logRatio = Math.log(each.getSlope());
                double logRatioRev = Math.log(each.getSlopeRev());
//System.out.println("==" + ratioSum + "\t" + logRatio + "\t"  + each.getSlope());
                logRatioSum += logRatio;
                logRatioSumRev += logRatioRev;

//		    System.out.println("==\t" + each.getSlope() + "\t" + logRatio);
                peptideCount++;
                quantifiedCount++;
            }

            double invLogRatio = 0;
            double invLogRatioRev = 0;
            if (peptideCount > 0) {
                averageRatio = ratioSum / peptideCount;
                averageRatioRev = ratioSumRev / peptideCount;
//		    System.out.println("---->>\t" + peptideCount + "\t" + ratioSum + "\t" + averageRatio );

                logAverageRatio = logRatioSum / peptideCount;
                logAverageRatioRev = logRatioSumRev / peptideCount;
                invLogRatio = Math.exp(logAverageRatio);
                invLogRatioRev = Math.exp(logAverageRatioRev);
            }

            proteinSb.append("P\t");
            proteinSb.append(protein.getLocus());
            proteinSb.append("\t");
            proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
            proteinSb.append("\t");
            proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
            proteinSb.append("\t");
		//let's use log value
            //proteinSb.append( averageRatio>0?CensusHelper.format.format(averageRatio):"NA" );
            //proteinSb.append("\t");
            //proteinSb.append( averageRatioRev>0?CensusHelper.format.format(averageRatioRev):"NA" );
            //proteinSb.append("\t");
            sproteinSb.append("P\t");
            sproteinSb.append(protein.getLocus());
            sproteinSb.append("\t");

            double devSum = 0;
            double devSumRev = 0;

            StringBuffer pepSb = new StringBuffer();

            WeightedProtein.ProteinModel pModel = new WeightedProtein.ProteinModel();

            double totalPepIntensity = 0;
            double areaRatioSum = 0;

            for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
                ChroPeptide each = tempItr.next();

                if (!noFilter && each.isFilterOut()) {
                    continue;
                }

                //if(noFilter || !each.isFilterOut())
                double dev = each.getSlope() - averageRatio;
                double devRev = each.getSlopeRev() - averageRatioRev;

                if ((Double.compare(each.getSlope(), Double.NaN) != 0)) {
                    devSum += dev * dev;
                    devSumRev += devRev * devRev;
                }

                pepSb.append("S\t");
                pepSb.append(each.isUnique() ? "U" : "");
                pepSb.append("\t");
                pepSb.append(each.getSequence());
                pepSb.append("\t");

                if (each.getCorr() < 0) {

                    if ((Double.compare(each.getSlope(), Double.NaN) == 0)) {
                        pepSb.append("NA\tNA\t");
                    } else {
                        pepSb.append(CensusHelper.format.format(each.getSlope()));
                        pepSb.append("\t");
                        pepSb.append(CensusHelper.format.format(each.getSlopeRev()));
                        pepSb.append("\t");
                    }

                    pepSb.append(CensusHelper.format.format(each.getCorr()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getDetValue()));
                } else {
                    pepSb.append(CensusHelper.format.format(each.getSlope()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getSlopeRev()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getCorr()));
                    pepSb.append("\t");
                    pepSb.append(CensusHelper.format.format(each.getDetValue()));
                }

//			System.out.println("====" + detSelect + "\t" + detValue + "\t" + each.getCorr() + "\t" +  (each.getCorr()*each.getCorr()));
                double corr = each.getCorr();
                int dataNum = each.getSpectraDataPoints();
                double tvalue = corr / Math.sqrt((1 - corr * corr) / (dataNum - 2));
                int df = dataNum - 1;

                if (df <= 0) {
                    continue;
                }
                TDistribution t = new TDistribution(df);
                double d = 1 - t.cumulativeProbability(tvalue);

                pepSb.append("\t");
                pepSb.append(d);
                pepSb.append("\t");
                pepSb.append((null == each.getXCorr()) ? "" : each.getXCorr());
                pepSb.append("\t");
                pepSb.append((null == each.getDeltCN()) ? "" : each.getDeltCN());
                pepSb.append("\t");
                pepSb.append(each.getSamIntensity());
                pepSb.append("\t");
                pepSb.append(each.getRefIntensity());
                pepSb.append("\t");
                double intRatio = (0 == each.getRefIntensity()) ? -1 : (each.getSamIntensity() / each.getRefIntensity());
                pepSb.append((intRatio >= 0) ? CensusHelper.format.format(intRatio) : "INF");
                pepSb.append("\t");

		    //       pepSb.append( CensusHelper.format.format(each.getSnRatio()) );
                //       pepSb.append("\t");
                pepSb.append(CensusHelper.format.format(each.getAnCompositeScore()));
                pepSb.append("\t");
                pepSb.append(each.getFileName()).append("\t");
                pepSb.append(each.getScanNum()).append("\t");
                pepSb.append(each.getChargeState());

                double enrichment = each.getEnrichment();

                if (enrichment > 0) {
                    pepSb.append("\t").append(enrichment);
                } else {
                    pepSb.append("\t").append("N/A");
                }

		    //each.getEnrichment()
		    //System.out.println("===>>" + each.getSpecCount());
                pepSb.append("\n");

                double rsqrtLog = Math.log(each.getDetValue());
                double stdevLog = -0.84 * rsqrtLog + 0.43;
                double invStdev = Math.exp(stdevLog);

                //grap only valid ratio
                if (each.getSlope() > 0) {
                    pModel.add(invStdev, each.getSlope());
                    // double dev = each.getSlope()-averageRatio;
                }

                //Hashtable<String, ArrayList> highScoreProteinHt = new Hashtable<String, ArrayList>();
                TDoubleArrayList al = highScoreProteinHt.get(protein.getLocus());

                if (null == al) {
                    al = new TDoubleArrayList();
                    al.add(each.getAreaRatio());
                    highScoreProteinHt.put(protein.getLocus(), al);
                } else {
                    al.add(each.getAreaRatio());
                }

                totalPepIntensity += each.getSamIntensity();
                totalPepIntensity += each.getRefIntensity();
                areaRatioSum += each.getAreaRatio();
            }

            ////////////add all none peptides here
            ChroProtein anProtein = anProteinHt.get(protein.getLocus());
            StringBuffer singlePepSb = new StringBuffer();
		//System.out.println(anProtein.getLocus());
            //anProtein.getPeptideList()
            int singletonPeptideCount = 0;
            int singletonUpCount = 0;
            int singletonDownCount = 0;

            if (!discardAN && null != anProtein) {
                for (Iterator<ChroPeptide> anPItr = anProtein.getPeptideList().iterator(); anPItr.hasNext();) {
                    ChroPeptide anPep = anPItr.next();
                    /*
                     if(anPep.getSequence().equals("R.VLLVIDEPHTDWAK.Y"))
                     {
                     }*/

                    if ((anPep.getAreaRatio() < allNoneUpperBound && anPep.getAreaRatio() > allNoneLowerBound)) {
                        continue;
                    }

                    if (detSelect && detValue < anPep.getDetValue()) {
                        continue;
                    }

                    if (anPep.getAnCompositeScore() >= allNoneCompositeScore) {

                        StringBuffer tmpSb = new StringBuffer();
                        tmpSb.append("&S").append("\t").append(anPep.isUnique() ? "U" : "").append("\t").append(anPep.getSequence()).append("\t");

                        if (anPep.getCorr() < 0) {
                            tmpSb.append("0.0");
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getCorr()));
                            tmpSb.append("\t");
                            tmpSb.append("0.0");
                        } else {
                            tmpSb.append(CensusHelper.format.format(anPep.getSlope()));
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getCorr()));
                            tmpSb.append("\t");
                            tmpSb.append(CensusHelper.format.format(anPep.getDetValue()));
                        }

                        tmpSb.append("\t");
                        tmpSb.append((null == anPep.getXCorr()) ? "" : anPep.getXCorr());
                        tmpSb.append("\t");
                        tmpSb.append((null == anPep.getDeltCN()) ? "" : anPep.getDeltCN());
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getSamIntensity());
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getRefIntensity());
                        tmpSb.append("\t");
                        double intRatio = (0 == anPep.getRefIntensity()) ? -1 : (anPep.getSamIntensity() / anPep.getRefIntensity());

                        if (anPep.getSamIntensity() >= anPep.getRefIntensity()) {
                            singletonUpCount++;
                        } else {
                            singletonDownCount++;
                        }

                        tmpSb.append((intRatio >= 0) ? CensusHelper.format.format(intRatio) : "INF");
                        //                                        tmpSb.append( CensusHelper.format.format(anPep.getSamIntensity()/anPep.getRefIntensity()) );
                        tmpSb.append("\t");
                        tmpSb.append(CensusHelper.format.format(anPep.getAnCompositeScore()));
                        tmpSb.append("\t");
                        tmpSb.append(anPep.getFileName()).append("\t");
                        tmpSb.append(anPep.getScanNum()).append("\t");
                        tmpSb.append(anPep.getChargeState()).append("\n");

			    //pepSb.append(tmpSb);
                        singletonPeptideCount++;
                        singlePepSb.append(tmpSb);
                    }
                }
            }

            //////////End of adding an peptides
            if (peptideCount > 1 && Double.compare(devSum, Double.NaN) != 0) {
                proteinSb.append(CensusHelper.format.format(Math.sqrt(devSum / (peptideCount - 1))));
                proteinSb.append("\t");
                proteinSb.append(CensusHelper.format.format(Math.sqrt(devSumRev / (peptideCount - 1))));
            } else {
                proteinSb.append("NA\tNA");
            }

            proteinSb.append("\t");
            if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
                proteinSb.append("NA\t");

            } else {
                proteinSb.append(CensusHelper.format.format(pModel.getStandardWeightedAverage()));
                proteinSb.append("\t");
            }

            proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
            proteinSb.append("\t");
            proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
            proteinSb.append("\t");
            proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
            proteinSb.append("\t");
            proteinSb.append(protein.getSpectrumCount());
            proteinSb.append("\t");
            proteinSb.append(protein.getLspectrumCount());
            proteinSb.append("\t");
            proteinSb.append(protein.getHspectrumCount());
            proteinSb.append("\t");
            //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
            proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
            proteinSb.append("\t");
            proteinSb.append(protein.getDescription());
            proteinSb.append("\n");

            sproteinSb.append(singletonPeptideCount);
            sproteinSb.append("\t");
            sproteinSb.append(protein.getSpectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(protein.getLspectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(protein.getHspectrumCount());
            sproteinSb.append("\t");
            sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
            sproteinSb.append("\t");
            sproteinSb.append(protein.getDescription());
            sproteinSb.append("\n");

            //here
            List redProList = protein.getRedunList();
            for (Iterator<ChroProtein> redPItr = redProList.iterator(); redPItr.hasNext();) {
                ChroProtein redPro = redPItr.next();

                proteinSb.append("P\t");
                proteinSb.append(redPro.getLocus());
                proteinSb.append("\t");
                proteinSb.append(averageRatio > 0 ? CensusHelper.format.format(averageRatio) : "NA");
                proteinSb.append("\t");
                proteinSb.append(averageRatioRev > 0 ? CensusHelper.format.format(averageRatioRev) : "NA");
                proteinSb.append("\t");

                if (peptideCount > 1) {
                    proteinSb.append(CensusHelper.format.format(Math.sqrt(devSum / (peptideCount - 1))));
                    proteinSb.append("\t");
                    proteinSb.append(CensusHelper.format.format(Math.sqrt(devSumRev / (peptideCount - 1))));
                } else {
                    proteinSb.append("NA\tNA");
                }

                proteinSb.append("\t");

                if ((Double.compare(pModel.getStandardWeightedAverage(), Double.NaN) == 0)) {
                    proteinSb.append("NA\t");

                } else {
                    proteinSb.append(CensusHelper.format.format(pModel.getStandardWeightedAverage()));
                    proteinSb.append("\t");
                }
                //proteinSb.append( CensusHelper.format.format(pModel.getStandardWeightedAverage()) );
                proteinSb.append(invLogRatio > 0 ? CensusHelper.format.format(invLogRatio) : "NA");
                proteinSb.append("\t");
                proteinSb.append(invLogRatioRev > 0 ? CensusHelper.format.format(invLogRatioRev) : "NA");
                proteinSb.append("\t");
                proteinSb.append(peptideCount > 0 ? peptideCount : "NA");
                proteinSb.append("\t");
                proteinSb.append(redPro.getSpectrumCount());
                proteinSb.append("\t");
                proteinSb.append(redPro.getLspectrumCount());
                proteinSb.append("\t");
                proteinSb.append(redPro.getHspectrumCount());
                proteinSb.append("\t");
                //proteinSb.append( peptideCount>0?(totalPepIntensity/peptideCount):"");
                proteinSb.append(peptideCount > 0 ? (areaRatioSum / peptideCount) : "NA");
                proteinSb.append("\t");
                proteinSb.append(redPro.getDescription());
                proteinSb.append("\n");

                sproteinSb.append("P\t");
                sproteinSb.append(redPro.getLocus());
                sproteinSb.append("\t");
                sproteinSb.append(singletonPeptideCount);
                sproteinSb.append("\t");
                sproteinSb.append(protein.getSpectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(protein.getLspectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(protein.getHspectrumCount());
                sproteinSb.append("\t");
                sproteinSb.append(singletonUpCount + "/" + singletonDownCount);
                sproteinSb.append("\t");
                sproteinSb.append(protein.getDescription());
                sproteinSb.append("\n");

            }

            percent += eachSeg;

            if (isGui && null != aJProgressBar) {
                aJProgressBar.setValue((int) percent);
            }

            if (singletonPeptideCount < allNoneMinPeptide) {
                singlePepSb = new StringBuffer();
            }

            if (pepSb.length() <= 0 && singlePepSb.length() <= 0) {
                redunProteinCount -= proteinGroupCount;
                continue;
            }

		//if( pepSb.length()<=0 && singlePepSb.length()<=0) {
            //if( pepSb.length()<=0 && singlePepSb.length()>0) {
            if (singlePepSb.length() > 0) {
                singletonResult.append(sproteinSb.toString());
                singletonResult.append(singlePepSb.toString());
                pepSb.append(singlePepSb.toString());
            }

            result.append(proteinSb.toString());
            result.append(pepSb.toString());

            uniqueProteinCount++;
        }

        rResult.setResult(result);
        rResult.setSingletonResult(singletonResult);
        rResult.setTotalCount(totalCount);
        rResult.setQuantifiedCount(quantifiedCount);
        rResult.setRedunProteinCount(redunProteinCount);
        rResult.setUniqueProteinCount(uniqueProteinCount);
        rResult.setProteinGroupCount(proteinGroupCount);

        p.print("H\t");
        p.print("Total Redundant Proteins\t");
        p.println(redunProteinCount);
        p.print("H\t");
        p.print("Total Unique Proteins\t");
        p.println(uniqueProteinCount);
        p.print("H\t");
        p.print("Total peptides\t");
        p.println(totalCount);
        p.print("H\t");
        p.print("Quantified peptides\t");
        p.print(quantifiedCount);
        p.print("\n");
        p.print("H\t");
        p.print("Quantification efficiency\t");
        p.print(CensusHelper.format.format((double) quantifiedCount / totalCount * 100));
        p.print(" %\n");
        p.print("H\t");
        p.print("Correction Factor (Ln)\t");
        p.println(correctFactorValue);

        singleP.print("H\t");
        singleP.print("Total Redundant Proteins\t");
        singleP.println(redunProteinCount);
        singleP.print("H\t");
        singleP.print("Total Unique Proteins\t");
        singleP.println(uniqueProteinCount);
        singleP.print("H\t");
        singleP.print("Total peptides\t");
        singleP.println(totalCount);

        p.print(result.toString());
        singleP.print(singletonResult.toString());

        return rResult;

    }

    //labeled report with all-none
    public void exportReportANGUI(
            final boolean noFilter,
            final boolean detSelect,
            final boolean pValueSelect,
            final boolean filterFragmentIons,
            final double detValue,
            final double pValue,
            final double correctFactorValue,
            final boolean isUniquePeptide,
            final boolean removeNegative,
            final boolean discardAN,
            final double allNoneLowerBound,
            final double allNoneUpperBound,
            final double allNoneCompositeScore,
            final int allNoneMinPeptide,
            final boolean isGui,
            final boolean discardUnlabeledPeptide,
            final boolean discardReverseProteins
    ) {

        File f = new File(this.currentDirectory + File.separator + "census-out.txt");
        File tmpFile = null;

        if (isGui) {
            JFileChooser choose = new JFileChooser(this.currentDirectory);

            choose.setMultiSelectionEnabled(false);
            choose.setDialogTitle("Save CenSus Report File");
            choose.setApproveButtonText("Write");
            choose.addChoosableFileFilter(new SimpleFileNameFilter("txt", "CenSus Report File (*.txt)"));

            if (currentDirectory != null && !"".equals(currentDirectory)) {
                choose.setCurrentDirectory(new File(currentDirectory));
            }

            choose.setSelectedFile(f);

            int returnVal = choose.showOpenDialog(chroPanel);
            tmpFile = choose.getSelectedFile();

            if (returnVal == choose.CANCEL_OPTION) {
                return;
            }

            this.currentDirectory = tmpFile.getAbsolutePath();

            this.currentDirectory = currentDirectory.substring(0, currentDirectory.lastIndexOf(File.separator));
        } else {
            tmpFile = new File(this.currentDirectory + File.separator + "census-out.txt");
        }

        final JDialog progress = new JDialog(this);
        final File file = tmpFile;

        isotopeFileField.setText(file.getAbsolutePath());

        final File singleFile = new File(this.currentDirectory + File.separator + "census-out_singleton.txt");

        final String curDir = this.currentDirectory;

        try {
//            final String tempFolder = this.currentDirectory + File.separator;
//            final ChroPeptide tempPeptide = this.currentPeptide;

            JProgressBar tmpBar = null;
            if (isGui) {
                //final JPanel tempPepPanel = this.peptidePanel;
                tmpBar = new JProgressBar(0, 100);
                tmpBar.setStringPainted(true);
                Container cp = progress.getContentPane();

                //cp.setSize(1000, 500);
                JLabel jb = new JLabel("Generating report...",
                        SwingConstants.CENTER);
                cp.add(jb, BorderLayout.SOUTH);
                cp.add(tmpBar, BorderLayout.NORTH);
                progress.setSize(500, 100);
                progress.setLocationRelativeTo(this);
                progress.pack();

                progress.setResizable(false);
                progress.setVisible(true);

            }
            final JProgressBar aJProgressBar = tmpBar;
            final boolean tempIsdataDependent = this.isDataDependent;

            ////////////  FIX THIS
            //final boolean tempIsFilterFrag = true; //options.isFilterFragmentIons();
            //////////////////
	    //final boolean tempIsFilterFrag = options.isFilterFragmentIons();
            //final float tmpDetValue = options.getDetFactorValue();
            //final double tmpDetValue = detValue;
	    //bar.setValue((int)percent);
            //        progressBar.setStringPainted(true);
            Thread t = new Thread() {
                private boolean isSuccessful = true;
                private String errorMessage = "";
                private PrintStream p = new PrintStream(new BufferedOutputStream(new FileOutputStream(file)));
                private PrintStream singleP = new PrintStream(new BufferedOutputStream(new FileOutputStream(singleFile)));

                public void run() {
                    ReportResult rResult = null;

                    try {
                        /*
                         */
                        ReportParam param = new ReportParam();
                        param.setProteinList(proteinList);
                        param.setConf(conf);
                        param.setDiscardAN(discardAN);
                        param.setNoFilter(noFilter);
                        param.setIsGui(isGui);
                        param.setAJProgressBar(aJProgressBar);
                        param.setDetValue(detValue);
                        param.setFilterFragmentIons(filterFragmentIons);
                        param.setCorrectFactorValue(correctFactorValue);
                        param.setDiscardUnlabeledPeptide(discardUnlabeledPeptide);
                        param.setDiscardReverseProtein(discardReverseProteins);
                        param.setRemoveNegative(removeNegative);
                        param.setDetSelect(detSelect);
                        param.setIsUniquePeptide(isUniquePeptide);
                        param.setPValueSelect(pValueSelect);
                        param.setPValue(pValue);
                        param.setDetValue(detValue);
                        param.setAllNoneLowerBound(allNoneLowerBound);
                        param.setAllNoneUpperBound(allNoneUpperBound);
                        param.setAllNoneCompositeScore(allNoneCompositeScore);
                        param.setAllNoneMinPeptide(allNoneMinPeptide);

                        //rResult = runReport(proteinList, conf, discardAN, noFilter, isGui, aJProgressBar, detValue, filterFragmentIons, correctFactorValue, discardUnlabeledPeptide);
                        rResult = runReport(param, p, singleP);
                    } catch (Exception e) {
                        isSuccessful = false;
                        errorMessage = e.getMessage();
                        e.printStackTrace();
                    }

                    int totalCount = rResult.getTotalCount();
                    int quantifiedCount = rResult.getQuantifiedCount();
                    int redunProteinCount = rResult.getRedunProteinCount();
                    int uniqueProteinCount = rResult.getUniqueProteinCount();
                    int proteinGroupCount = rResult.getProteinGroupCount();

                    final int finalTotalCount = totalCount;
                    final int finalQuantifiedCount = quantifiedCount;

                    SwingUtilities.invokeLater(new Runnable() {
                        public void run() {
                            progress.setVisible(false);
                            progress.hide();

                            if (isSuccessful) {
                                JOptionPane.showMessageDialog(
                                        chroPanel,
                                        "Report file was successfully created."
                                        + "\n\nTotal peptides : " + finalTotalCount
                                        + "\nQuantified peptides : " + finalQuantifiedCount
                                        + "\nQuantification efficiency : " + CensusHelper.format.format((double) finalQuantifiedCount / finalTotalCount * 100) + " %",
                                        "Report file Creation",
                                        JOptionPane.PLAIN_MESSAGE);
                            } else {
                                JOptionPane.showMessageDialog(
                                        chroPanel,
                                        errorMessage,
                                        "Report file Creation",
                                        JOptionPane.ERROR_MESSAGE);
                            }

                        }

                    }
                    );

                    if (null != p) {
                        p.close();
                    }

                    if (null != singleP) {
                        singleP.close();
                    }

                }
            };

            t.start();

        } catch (IOException e) {
            System.out.println("Failed to write file" + e);
            JOptionPane.showMessageDialog(this, e.toString(), e.toString(), JOptionPane.ERROR_MESSAGE);

            if (null != progress) {
                progress.setVisible(false);
                progress.hide();
            }

            e.printStackTrace();

        }
    }

    private void reportItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_reportItemActionPerformed
        // TODO add your handling code here:

        //Configuration conf = Configuration.getInstance();
        //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.
        if (cr.getExpType() > 0) {
            switch (cr.getExpType()) {
                case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
                    ExportTandemTagReport tdialog = new ExportTandemTagReport(this, true);
                    tdialog.pack();
                    tdialog.setLocationRelativeTo(this);
                    tdialog.setVisible(true);
                    tdialog.setResizable(false);
                    break;

                case CensusConstants.MSMS_DATA_INDEPENDENT:
                    ExportDialog dialog = new ExportDialog(this, true);
                    dialog.pack();
                    dialog.setLocationRelativeTo(this);
                    dialog.setVisible(true);
                    dialog.setResizable(false);

                default:
                    break;
            }
        } else if (this.isLabeled()) {
            ExportDialog dialog = new ExportDialog(this, true);
            dialog.pack();
            dialog.setLocationRelativeTo(this);
            dialog.setVisible(true);
            dialog.setResizable(false);
        } else {
            ExportNonlabelDialog dialog = new ExportNonlabelDialog(this, true);
            dialog.pack();
            dialog.setLocationRelativeTo(this);
            dialog.setVisible(true);
            dialog.setResizable(false);

        }

    }//GEN-LAST:event_reportItemActionPerformed

    private void fileMenuActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fileMenuActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_fileMenuActionPerformed

    private void drawPepDist(String data) {
        try {

            pepDistPlot.clear(true);
            pepDistPlot.repaint();
            peptideDistParser.parse(null, data);
            peptideDistParser.postProcessing();

            //pepDistPlot._drawPoints();
        } catch (Exception e) {
            System.out.println("error==>> : " + e);
        }

    }

    private void drawPlot(String data) {
        try {
            //this.plot = new Plot();
            //PlotBoxMLParser chroXmlParser = new PlotMLParser((Plot)plot);
            plot.clear(true);
            plot.repaint();
            chroXmlParser.parse(null, data);
            chroXmlParser.postProcessing();
        } catch (Exception e) {
            System.out.println("error==>> : " + e);
        }
    }

    /**
     * ** Sample DATA ****** xmlData = new StringBuffer();
     * xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>");
     * xmlData.append("</plot>"); xmlData.append("<title>Chromatogram</title>");
     * xmlData.append("<xLabel>Scan Number</xLabel>");
     * xmlData.append("<yLabel>Intensity</yLabel>");
     * xmlData.append("<noGrid/>");
     * xmlData.append("<size width=\"550\" height=\"370\"/>"); - <xTicks>
     * <tick label="50-60" position="50" />
     * <tick label="60-70" position="60" />
     * <tick label="70-80" position="70" />
     * <tick label="80-90" position="80" />
     * <tick label="90-100" position="90" />
     * </xTicks>
     * xmlData.append("<dataset name=\"dmso\" marks=\"none\" connected=\"yes\" stems=\"no\">");
     * xmlData.append("<p x=\"221.724\" y=\"7719.56\"/>");
     * xmlData.append("<p x=\"220.945\" y=\"7616.44\"/>");
     * xmlData.append("<p x=\"220.168\" y=\"7822.69\"/>");
     * xmlData.append("<p x=\"219.392\" y=\"8073.13\"/>");
     * xmlData.append("<p x=\"131.194\" y=\"2519.17\"/>");
     * xmlData.append("</dataset>"); xmlData.append("</plot>");
     *
     *
     ****************************
     */
    public void generateNonLabelData(ChroPeptide peptide) {
        generateNonLabelData(peptide, -1);
    }

    public void generateNonLabelData(ChroPeptide peptide, int selectId) {

        if (null == currentPeptide) {
            return;
        }

        if (null == peptide) {
            peptide = currentPeptide;
        }

        StringBuffer xmlData = new StringBuffer();
        xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>");
        xmlData.append("<plot>");
        xmlData.append("<title>Chromatogram</title>");

        if (conf.isAlign()) {
            xmlData.append("<xLabel>Scan Number</xLabel>");
        } else {
            xmlData.append("<xLabel>Retention Time</xLabel>");
        }

        xmlData.append("<yLabel>Intensity</yLabel>");

        int startRange = 0;
        int endRange = 0;

        xmlData.append("<startRange>").append(peptide.getStartRange()).append("</startRange>");
        xmlData.append("<endRange>").append(peptide.getEndRange()).append("</endRange>");
        xmlData.append("<scanNum>").append(peptide.getScanNum()).append("</scanNum>");
        //xmlData.append("<scanNum>").append("9000").append("</scanNum>");

        startRange = Integer.parseInt(peptide.getStartRange());
        endRange = Integer.parseInt(peptide.getEndRange());

        xmlData.append("<grid/>");

        int fileSize = cr.getFileList().size();

        //System.out.println(cr.getFileList());
        String idFileName = peptide.getFileName();
        idFileName = idFileName.substring(0, idFileName.indexOf("."));
        int idFileIndex = 0;

        StringBuffer ticks = new StringBuffer();
        StringBuffer[] bufferArr = new StringBuffer[fileSize];

        ArrayList grayoutList = new ArrayList();
        for (int i = 0; i < bufferArr.length; i++) {
            bufferArr[i] = new StringBuffer();
            String legendName = cr.getFileName(i);

            if (legendName.contains(idFileName)) {
                idFileIndex = i;
            }

            if (legendName.length() > 10) {
                legendName = legendName.substring(0, 10);
            }
            bufferArr[i].append("<dataset name=\"").append(legendName).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");

            if (selectId >= 0 && selectId != i) {
                grayoutList.add(i);
            }
        }

        this.plot.setGrayOutList(grayoutList);

        ticks.append("<xTicks>");

        List dataList = peptide.getDataList();

        long[] intensitySumArr = new long[fileSize];

        int convertedScan = 0;

        if (cr.getQuantLevel() == 1) {
            for (Iterator<ChroNonLabelData> itr = dataList.iterator(); itr.hasNext();) {
                ChroNonLabelData eachData = itr.next();

                int[] scanNumArr = eachData.getScanNumArr();
                long[] intenArr = eachData.getIntensityArr();

                for (int i = 0; i < scanNumArr.length; i++) {
                    bufferArr[i].append("<p x=\"");
                    bufferArr[i].append(scanNumArr[0]);
                    bufferArr[i].append("\" y=\"");
                    bufferArr[i].append(intenArr[i]);
                    bufferArr[i].append("\"/>");

                    if (scanNumArr[0] >= startRange && scanNumArr[0] <= endRange) {
                        intensitySumArr[i] += intenArr[i];
                    }

                }

                if (peptide.getScanNum() < scanNumArr[idFileIndex] && convertedScan <= 0) {
                    convertedScan = scanNumArr[0];
                }

                //System.out.println(convertedScan + " " + scanNumArr[idFileIndex] + " " + peptide.getScanNum() + " " + idFileIndex);
                ticks.append("<tick label=\"");
                ticks.append(scanNumArr[0]);
                ticks.append("\" position=\"");
                ticks.append(scanNumArr[0]);
                ticks.append("\" />");
            }
        } else if (cr.getQuantLevel() == 2) {
            for (Iterator<ChroNonLabelMSMSData> itr = dataList.iterator(); itr.hasNext();) {
                ChroNonLabelMSMSData eachData = itr.next();

                int[] scanNumArr = eachData.getScanArr();
                long[] intenArr = eachData.getTotalIntArr();

                for (int i = 0; i < scanNumArr.length; i++) {
                    bufferArr[i].append("<p x=\"");
                    bufferArr[i].append(scanNumArr[0]);
                    bufferArr[i].append("\" y=\"");
                    bufferArr[i].append(intenArr[i]);
                    bufferArr[i].append("\"/>");

                    if (scanNumArr[0] >= startRange && scanNumArr[0] <= endRange) {
                        intensitySumArr[i] += intenArr[i];
                    }
                }

                ticks.append("<tick label=\"");
                ticks.append(scanNumArr[0]);
                ticks.append("\" position=\"");
                ticks.append(scanNumArr[0]);
                ticks.append("\" />");
            }

        }

        ticks.append("</xTicks>");

        String prevDir = "";
        String curDir = "";

        int rowCount = nonlabelTableModel.getRowCount();

        for (int i = 0; i < rowCount; i++) {
            nonlabelTableModel.removeRow(0); //note zero here
        }
        rowCount = nonlabelSummaryTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            nonlabelSummaryTableModel.removeRow(0); //note zero here
        }

        Hashtable<String, ChroXmlReader.Sample> sampleHt = new Hashtable<String, ChroXmlReader.Sample>();

        /*
         if(selectId>=0)
         {
         bufferArr[selectId].append("</dataset>");
         xmlData.append(bufferArr[selectId].toString());
         }
         */
        for (int i = 0; i < bufferArr.length; i++) {
          //  if(selectId<0 || selectId!=i)

            bufferArr[i].append("</dataset>");
            xmlData.append(bufferArr[i].toString());

            //JLabel label = new JLabel(cr.getFileList().get(i));
            curDir = cr.getFileList().get(i);
            String sampleName = cr.getSampleName(curDir);
            Vector vec = new Vector();
            vec.add(sampleName);
            vec.add(curDir.substring(curDir.lastIndexOf(File.separator) + 1));

            //we don't know if file path is generated from linux or window.
            if (curDir.startsWith("/")) {
                vec.add(curDir.substring(0, curDir.lastIndexOf("/")));
            } else {
                vec.add(curDir.substring(0, curDir.lastIndexOf("\\")));
            }

            vec.add(CensusHelper.scientificFormat.format(intensitySumArr[i]));

            this.nonlabelTableModel.addRow(vec);

            ChroXmlReader.Sample sample = sampleHt.get(sampleName);

            if (null == sample) {
                sample = new ChroXmlReader.Sample(sampleName);
                sample.addIntensity(intensitySumArr[i]);
                sampleHt.put(sampleName, sample);
            } else {
                sample.addIntensity(intensitySumArr[i]);
            }

            prevDir = curDir;
        }

        xmlData.append("</plot>");

        //this.plot.setScanNum(peptide.getScanNum());
        this.plot.setScanNum(convertedScan);

        this.plot.setDtaStartRange(peptide.getDtaStartRange());
        this.plot.setDtaEndRange(peptide.getDtaEndRange());

        drawPlot(xmlData.toString());

        Hashtable<String, ChroXmlReader.Sample> sampleObjList = cr.getSampleObjList();

        for (Iterator<String> itr = cr.getSampleList().iterator(); itr.hasNext();) {
            String sampleName = itr.next();

            //ChroXmlReader.Sample eachSam = sampleObjList.get(sampleName);
            ChroXmlReader.Sample eachSam = sampleHt.get(sampleName);

            //System.out.println("sum===>>" + sampleName + " " + eachSam.getSumIntensity());
            Vector vec = new Vector();
            vec.add(sampleName);

            vec.add(CensusHelper.scientificFormat.format(eachSam.getSumIntensity()));
            vec.add(CensusHelper.scientificFormat.format(eachSam.getAverage()));
            //System.out.println("====" + eachSam.getIntensityArr().size());

            long[] intArr = new long[eachSam.getIntensityArr().size()];
            int count = 0;
            for (Iterator<Long> itrInt = eachSam.getIntensityArr().iterator(); itrInt.hasNext();) {
                //itrInt.next();
                intArr[count++] = itrInt.next().longValue();
            }

            vec.add(CensusHelper.scientificFormat.format(STDev.getStdev(intArr)));

            this.nonlabelSummaryTableModel.addRow(vec);
        }

    }

    public String chromXMLHeader(ChroPeptide peptide) {
        StringBuffer sb = new StringBuffer();
        sb.append("<?xml version=\"1.0\" standalone=\"no\"?>");
        sb.append("<plot>");
        sb.append("<title>Chromatogram</title>");
        sb.append("<xLabel>Scan Number</xLabel>");
        sb.append("<yLabel>Intensity</yLabel>");

        sb.append("<startRange>").append(peptide.getStartRange()).append("</startRange>");
        sb.append("<endRange>").append(peptide.getEndRange()).append("</endRange>");
        sb.append("<scanNum>").append(peptide.getScanNum()).append("</scanNum>");

        sb.append("<grid/>");

        return sb.toString();

    }

    public void generateITRAQData(ChroPeptide peptide) {
        if (null == currentPeptide) {
            return;
        }

        if (null == peptide) {
            peptide = this.currentPeptide;
        }

        StringBuffer xmlData = new StringBuffer();
        xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>");
        xmlData.append("<plot>");
        xmlData.append("<title>Chromatogram</title>");

        if (conf.isAlign()) {
            xmlData.append("<xLabel>Scan Number</xLabel>");
        } else {
            xmlData.append("<xLabel>Retention Time</xLabel>");
        }

        xmlData.append("<yLabel>Intensity</yLabel>");

        int startRange = 0;
        int endRange = 0;

        xmlData.append("<startRange>").append(peptide.getStartRange()).append("</startRange>");
        xmlData.append("<endRange>").append(peptide.getEndRange()).append("</endRange>");
        xmlData.append("<scanNum>").append(peptide.getScanNum()).append("</scanNum>");

        if (null != peptide.getStartRange() && null != peptide.getEndRange()) {
            startRange = Integer.parseInt(peptide.getStartRange());
            endRange = Integer.parseInt(peptide.getEndRange());
        }

        xmlData.append("<grid/>");

        //xmlData.append("<size width=\"590\" height=\"390\"/>");
        //int fileSize = cr.getFileList().size();
        double[] massMonitorArr = peptide.getMassMonitorArr();

        StringBuffer ticks = new StringBuffer();

        StringBuffer[] bufferArr = null;
        List dataList = peptide.getDataList();

        switch (cr.getExpType()) {
            case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                bufferArr = new StringBuffer[1];

                bufferArr[0] = new StringBuffer();
                bufferArr[0].append("<dataset name=\"m/z ").append("scan num ").append(peptide.getScanNum()).append("\" marks=\"none\" connected=\"no\" stems=\"yes\">");

                ticks.append("<xTicks>");

                ChroiTRAQLabelData data = (ChroiTRAQLabelData) dataList.get(0);

                int tscanNum = data.getScanNum();
                long[] tintenArr = data.getIntensityArr();

                for (int i = 0; i < tintenArr.length; i++) {
                    bufferArr[0].append("<p x=\"");
                    bufferArr[0].append(massMonitorArr[i]);
                    bufferArr[0].append("\" y=\"");
                    bufferArr[0].append(tintenArr[i]);
                    bufferArr[0].append("\"/>");

                    ticks.append("<tick label=\"");
                    ticks.append(massMonitorArr[i]);
                    ticks.append("\" position=\"");
                    ticks.append(massMonitorArr[i]);
                    ticks.append("\" />");

                }

                break;

            case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ

                bufferArr = new StringBuffer[massMonitorArr.length];

                //ArrayList grayoutList = new ArrayList();
                for (int i = 0; i < bufferArr.length; i++) {
                    bufferArr[i] = new StringBuffer();
                    bufferArr[i].append("<dataset name=\"m/z ").append(massMonitorArr[i]).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");

                }

                ticks.append("<xTicks>");

                for (Iterator<ChroiTRAQLabelData> itr = dataList.iterator(); itr.hasNext();) {
                    ChroiTRAQLabelData eachData = itr.next();

                    int scanNum = eachData.getScanNum();
                    long[] intenArr = eachData.getIntensityArr();

                    for (int i = 0; i < intenArr.length; i++) {
                        bufferArr[i].append("<p x=\"");
                        bufferArr[i].append(scanNum);
                        bufferArr[i].append("\" y=\"");
                        bufferArr[i].append(intenArr[i]);
                        bufferArr[i].append("\"/>");
                    }

                    ticks.append("<tick label=\"");
                    ticks.append(scanNum);
                    ticks.append("\" position=\"");
                    ticks.append(scanNum);
                    ticks.append("\" />");
                }

                break;

            default:
                break;
        }

        ticks.append("</xTicks>");

        String prevDir = "";
        String curDir = "";

        /*
         LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex!=0)?endIndex:samArr.length, conf.getMaxSpectrumShift());

         double slope = reg.getSlope();
         double intercept = reg.getIntercept();

         corrPlot.setData(refArr, samArr, startIndex, (endIndex!=0)?endIndex:samArr.length, slope, intercept, reg.getBestShift());
         corrPlot.clear(true);
         corrPlot.repaint();
         correlationPanel.add(corrPlot);

         //proteinRatioDistPanel.add(corrPlot);

         //this.corrCoeffField.setText( CensusHelper.format.format(reg.getCorr()) );
         if(reg.getCorr()<0)
         {
         this.rrField.setText("N/A");
         this.areaRatioLogField.setText("N/A");
         this.regressionRatioField.setText("N/A");
         this.areaRatioField.setText("N/A");
         //this.snRatioField.setText("N/A");
         }
         else
         {
         this.rrField.setText( CensusHelper.format.format(reg.getCorr()*reg.getCorr()));
         this.areaRatioLogField.setText( CensusHelper.format.format( Math.log(slope) ));
         this.regressionRatioField.setText( CensusHelper.format.format(slope) );
         this.areaRatioField.setText( CensusHelper.format.format(reg.getAreaRatio()) );
         //this.snRatioField.setText( CensusHelper.format.format(snRatio) );
         }

         this.shiftField.setText( String.valueOf(reg.getBestShift()) );
         */
        Hashtable<String, ChroXmlReader.Sample> sampleHt = new Hashtable<String, ChroXmlReader.Sample>();

        for (int i = 0; i < bufferArr.length; i++) {
          //  if(selectId<0 || selectId!=i)

            bufferArr[i].append("</dataset>");
            xmlData.append(bufferArr[i].toString());

        }

        xmlData.append("</plot>");

        this.plot.setScanNum(peptide.getScanNum());
        this.plot.setDtaStartRange(peptide.getDtaStartRange());
        this.plot.setDtaEndRange(peptide.getDtaEndRange());

        drawPlot(xmlData.toString());
    }

    public void generateDepData(ChroPeptide peptide) {
        if (null == peptide) {
            peptide = currentPeptide;
        }

        if (null == currentPeptide) {
            return;
        }

        List l = peptide.getDataList();
        StringBuffer xmlData = new StringBuffer();

        /*
         xmlData.append("<?xml version=\"1.0\" standalone=\"no\"?>");
         xmlData.append("<plot>");
         xmlData.append("<title>Chromatogram</title>");
         xmlData.append("<xLabel>Scan Number</xLabel>");
         xmlData.append("<yLabel>Intensity</yLabel>");
         */
        xmlData.append(chromXMLHeader(peptide));

        /*
         xmlData.append("<startRange>").append(peptide.getStartRange()).append("</startRange>");
         xmlData.append("<endRange>").append(peptide.getEndRange()).append("</endRange>");
         xmlData.append("<scanNum>").append(peptide.getScanNum()).append("</scanNum>");
         */
        int startRange = Integer.parseInt(peptide.getStartRange());
        int endRange = Integer.parseInt(peptide.getEndRange());

        StringBuffer sampleData = new StringBuffer();
        StringBuffer refData = new StringBuffer();
        StringBuffer ticks = new StringBuffer();
        sampleData.append("<dataset name=\"sample\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        refData.append("<dataset name=\"reference\" marks=\"none\" connected=\"yes\" stems=\"no\">");

        ticks.append("<xTicks>");

        long[] samArr = new long[l.size()];
        long[] refArr = new long[samArr.length];

        int index = 0;
        int startIndex = 0;
        int endIndex = 0;

        for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
            ChroData data = itr.next();

            ///rpark
            samArr[index] = data.getSampleIntensity();
            refArr[index] = data.getRefIntensity();

            int scanTemp = data.getScanNum();

            if (startRange >= scanTemp) {
                startIndex = index;
            }
            if (endRange >= scanTemp) {
                endIndex = index;
            }

            index++;

            //<tick label="90-100" position="90" />
            sampleData.append("<p x=\"");
            sampleData.append(data.getScanNum());
            sampleData.append("\" y=\"");
            sampleData.append(data.getSampleIntensity());
            sampleData.append("\"/>");

            refData.append("<p x=\"");
            refData.append(data.getScanNum());
            refData.append("\" y=\"");
            refData.append(data.getRefIntensity());
            refData.append("\"/>");

            ticks.append("<tick label=\"");
            ticks.append(data.getScanNum());
            ticks.append("\" position=\"");
            ticks.append(data.getScanNum());
            ticks.append("\" />");

        }

        LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

        double slope = reg.getSlope();
        double intercept = reg.getIntercept();
        endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);

        corrPlot.setData(refArr, samArr, startIndex, endIndex, slope, intercept, reg.getBestShift());
        corrPlot.clear(true);
        corrPlot.repaint();
        correlationPanel.add(corrPlot);

        //proteinRatioDistPanel.add(corrPlot);
        //this.corrCoeffField.setText( CensusHelper.format.format(reg.getCorr()) );
        int dataNum = endIndex - startIndex + 1;
        this.jTextField1.setText(String.valueOf(dataNum));

        if (reg.getCorr() < 0) {
            this.regScoreField.setText("N/A");
            this.rrField.setText("N/A");
            this.areaRatioLogField.setText("N/A");
            this.regressionRatioField.setText("N/A");
            this.areaRatioField.setText("N/A");
            this.pvalueField.setText("N/A");
        } else {
            this.regScoreField.setText(CensusHelper.format.format(reg.getCorr()));
            this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));


            this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
            this.regressionRatioField.setText(CensusHelper.format.format(slope));
            this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));

            double corr = reg.getCorr();

            double tvalue = corr / Math.sqrt((1 - corr * corr) / (dataNum - 2));
            int df = dataNum - 1;
            TDistribution t = new TDistribution(df);

            try {
                double pvalue = 1 - t.cumulativeProbability(tvalue);
                pvalue = pvalue * 2; //two tails
                this.pvalueField.setText(CensusHelper.d5format.format(pvalue));

            } catch (Exception e) {

            }
            //this.snRatioField.setText( CensusHelper.format.format(snRatio) );
        }

        this.shiftField.setText(String.valueOf(reg.getBestShift()));

        ticks.append("</xTicks>");
        sampleData.append("</dataset>");
        refData.append("</dataset>");

        xmlData.append(ticks.toString());
        xmlData.append(sampleData.toString());
        xmlData.append(refData.toString());

        xmlData.append("</plot>");

        this.plot.setScanNum(peptide.getScanNum());
        this.plot.setDtaStartRange(peptide.getDtaStartRange());
        this.plot.setDtaEndRange(peptide.getDtaEndRange());

        drawPlot(xmlData.toString());

    }

    private DataIndepModel calcFragIons(ChroPeptide peptide) {

        List l = peptide.getDataList();

        int startRange = Integer.parseInt(peptide.getStartRange());
        int endRange = Integer.parseInt(peptide.getEndRange());

        int pepLength = ((ChroData) l.get(0)).getResidueLength();

        long[] samArr = new long[l.size()];
        long[] refArr = new long[samArr.length];
        int[] scanNumArr = new int[samArr.length];

        int index = 0;
        int startIndex = 0;
        int endIndex = 0;

        long[][] bsTempArr = new long[pepLength][samArr.length];
        long[][] ysTempArr = new long[pepLength][samArr.length];
        long[][] brTempArr = new long[pepLength][samArr.length];
        long[][] yrTempArr = new long[pepLength][samArr.length];

        for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
            ChroData data = itr.next();

            long bsArr[] = data.getBsIntensity();
            long ysArr[] = data.getYsIntensity();
            long brArr[] = data.getBrIntensity();
            long yrArr[] = data.getYrIntensity();

            for (int i = 0; i < bsArr.length; i++) {
                bsTempArr[i][index] = bsArr[i];
                ysTempArr[i][index] = ysArr[i];
                brTempArr[i][index] = brArr[i];
                yrTempArr[i][index] = yrArr[i];
            }

            scanNumArr[index] = data.getScanNum();

            int scanTemp = data.getScanNum();
            if (startRange >= scanTemp) {
                startIndex = index;
            }
            if (endRange >= scanTemp) {
                endIndex = index;
            }

            index++;

        }

        DataIndepModel model = new DataIndepModel();
        model.startIndex = startIndex;
        model.endIndex = endIndex;

        model.samArr = samArr;
        model.refArr = refArr;
        model.scanNumArr = scanNumArr;

        model.bsTempArr = bsTempArr;
        model.ysTempArr = ysTempArr;
        model.brTempArr = brTempArr;
        model.yrTempArr = yrTempArr;

        return model;

        //
        //return ionList;
    }

    public boolean isLabeled() {
        return labeled;
    }

    public void setLabeled(boolean labeled) {
        this.labeled = labeled;
    }

    public javax.swing.JScrollPane getFragIonScrollPanel() {
        return fragIonScrollPanel;
    }

    public void setFragIonScrollPanel(javax.swing.JScrollPane fragIonScrollPanel) {
        this.fragIonScrollPanel = fragIonScrollPanel;
    }

    public javax.swing.JPanel getCorrelationPanel() {
        return correlationPanel;
    }

    public void setCorrelationPanel(javax.swing.JPanel correlationPanel) {
        this.correlationPanel = correlationPanel;
    }

    public javax.swing.JTabbedPane getPepTabbedPanel() {
        return pepTabbedPanel;
    }

    public void setPepTabbedPanel(javax.swing.JTabbedPane pepTabbedPanel) {
        this.pepTabbedPanel = pepTabbedPanel;
    }

    public javax.swing.JPanel getParamPanel() {
        return paramPanel;
    }

    public void setParamPanel(javax.swing.JPanel paramPanel) {
        this.paramPanel = paramPanel;
    }

    public void setProteinList(ArrayList<ChroProtein> proteinList) {
        this.proteinList = proteinList;
    }

    private class DataIndepModel {

        public int startIndex;
        public int endIndex;

        long[] samArr;
        long[] refArr;
        int[] scanNumArr;

        long[][] bsTempArr;
        long[][] ysTempArr;
        long[][] brTempArr;
        long[][] yrTempArr;

    }

    public void generateInDepFragData(ChroPeptide peptide) {

        if (null == peptide) {
            peptide = currentPeptide;
        }

        if (null == currentPeptide) {
            return;
        }

        fragIonPanel.removeAll();

        StringBuffer xmlData = new StringBuffer();

        xmlData.append(chromXMLHeader(peptide));

        //xmlData.append("<size width=\"590\" height=\"390\"/>");
        StringBuffer sampleData = new StringBuffer();
        StringBuffer refData = new StringBuffer();
        StringBuffer ticks = new StringBuffer();
        sampleData.append("<dataset name=\"sample\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        refData.append("<dataset name=\"reference\" marks=\"none\" connected=\"yes\" stems=\"no\">");

        ticks.append("<xTicks>");

        DataIndepModel model = this.calcFragIons(peptide);

        int startIndex = model.startIndex;
        int endIndex = model.endIndex;
        long[] samArr = model.samArr;
        long[] refArr = model.refArr;
        int[] scanNumArr = model.scanNumArr;

//        FragIonList ionList = CalcUtil.getBestFragIons(bsTempArr, ysTempArr, brTempArr, yrTempArr, startIndex, endIndex, conf.getMaxSpectrumShift());
        FragIonList ionList = CalcUtil.getBestFragIons(model.bsTempArr, model.ysTempArr, model.brTempArr, model.yrTempArr, model.startIndex, model.endIndex, conf.getMaxSpectrumShift());

        //List<FragIon> ionList = CalcUtil.getBestFragIons(bsTempArr, ysTempArr, brTempArr, yrTempArr, startIndex, endIndex, conf.getMaxSpectrumShift());
        int listSize = ionList.size();
        StringBuffer[] ionSb = new StringBuffer[(ionList.getBestIndex() + 1) * 2];
        //StringBuffer[] ionSb = new StringBuffer[listSize*2];
        int tempIndex = 0;

        for (Iterator<FragIon> itr = ionList.iterator(); itr.hasNext();) {
            FragIon ion = itr.next();

            ionSb[tempIndex * 2] = new StringBuffer();
            ionSb[tempIndex * 2 + 1] = new StringBuffer();
            ionSb[tempIndex * 2].append("<dataset name=\"").append(ion.isBion() ? "bSample" : "ySample").append(ion.getIndex()).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");
            ionSb[tempIndex * 2 + 1].append("<dataset name=\"").append(ion.isBion() ? "bRef" : "yRef").append(ion.getIndex()).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");

            long[] tempSArr = ion.getSArr();
            long[] tempRArr = ion.getRArr();

            String ionName = ion.isBion() ? "b" : "y";
            ionName += ion.getIndex();

            FragIonPlot fPlot = new FragIonPlot(ionName, tempSArr, tempRArr, startIndex, (endIndex != 0) ? endIndex : (tempSArr.length - 1), true);
            fragIonPanel.add(fPlot, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 0 + (tempIndex * 60), 270, 60));

            for (int i = 0; i < tempSArr.length; i++) {
                samArr[i] += tempSArr[i];
                refArr[i] += tempRArr[i];

                ionSb[tempIndex * 2].append("<p x=\"");
                ionSb[tempIndex * 2].append(scanNumArr[i]);
                ionSb[tempIndex * 2].append("\" y=\"");
                ionSb[tempIndex * 2].append(tempSArr[i]);
                ionSb[tempIndex * 2].append("\"/>");

                ionSb[tempIndex * 2 + 1].append("<p x=\"");
                ionSb[tempIndex * 2 + 1].append(scanNumArr[i]);
                ionSb[tempIndex * 2 + 1].append("\" y=\"");
                ionSb[tempIndex * 2 + 1].append(tempRArr[i]);
                ionSb[tempIndex * 2 + 1].append("\"/>");
            }

            if (tempIndex == ionList.getBestIndex()) {
                break;
            }

            tempIndex++;
        }

        // display low quality of fragions
        for (int i = tempIndex + 1; i < listSize; i++) {
            FragIon ion = ionList.get(i);

            long[] tempSArr = ion.getSArr();
            long[] tempRArr = ion.getRArr();

            String ionName = ion.isBion() ? "b" : "y";
            ionName += ion.getIndex();

            FragIonPlot fPlot = new FragIonPlot(ionName, tempSArr, tempRArr, startIndex, (endIndex != 0) ? endIndex : (tempSArr.length - 1), false);
            fragIonPanel.add(fPlot, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 0 + (i * 60), 270, 60));
        }

        this.fragIonPanel.invalidate();
        this.fragIonPanel.validate();
        this.fragIonPanel.repaint();

        //this.sigToNoisePanel.removeAll();
        /*
         //SigNoisePlot sigPlot = new SigNoisePlot(this.sigToNoisePanel.getWidth(), this.sigToNoisePanel.getHeight(), ionList, tempIndex);
         sigPlot.setBackground(new Color(255, 255, 255));
         this.sigToNoisePanel.add(sigPlot);
         this.sigToNoisePanel.invalidate();
         this.sigToNoisePanel.validate();
         this.sigToNoisePanel.repaint();
         */
        for (int i = 0; i < scanNumArr.length; i++) {
            //<tick label="90-100" position="90" />
            sampleData.append("<p x=\"");
            sampleData.append(scanNumArr[i]);
            sampleData.append("\" y=\"");
            sampleData.append(samArr[i]);
            sampleData.append("\"/>");

            refData.append("<p x=\"");
            refData.append(scanNumArr[i]);
            refData.append("\" y=\"");
            refData.append(refArr[i]);
            refData.append("\"/>");

            ticks.append("<tick label=\"");
            ticks.append(scanNumArr[i]);
            ticks.append("\" position=\"");
            ticks.append(scanNumArr[i]);
            ticks.append("\" />");
        }

        LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

//
//	           for(int ii=0;ii<samArr.length;ii++)
//		                  System.out.println(samArr[ii] + " " + refArr[ii]);
        double slope = reg.getSlope();
        double intercept = reg.getIntercept();

        corrPlot.setData(refArr, samArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), slope, intercept, reg.getBestShift());
        corrPlot.clear(true);
        corrPlot.repaint();
        correlationPanel.add(corrPlot);
        //proteinRatioDistPanel.add(corrPlot);

        if (reg.getCorr() < 0) {
            this.rrField.setText("N/A");
            this.areaRatioLogField.setText("N/A");
            this.regressionRatioField.setText("N/A");
            this.areaRatioField.setText("N/A");
        } else {
            this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));
            this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
            this.regressionRatioField.setText(CensusHelper.format.format(slope));
            this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));
        }

        this.shiftField.setText(String.valueOf(reg.getBestShift()));

        ticks.append("</xTicks>");
        sampleData.append("</dataset>");
        refData.append("</dataset>");

        xmlData.append(ticks.toString());
        xmlData.append(sampleData.toString());
        xmlData.append(refData.toString());

        PostOptions options = PostOptions.getInstance();
        if (options.isDisplayFragmentIons()) {
            for (int i = 0; i < ionSb.length; i++) {
                ionSb[i].append("</dataset>");
                xmlData.append(ionSb[i].toString());
                //      System.out.print(".");
            }
        }

        xmlData.append("</plot>");

        this.plot.setScanNum(peptide.getScanNum());
        this.plot.setDtaStartRange(peptide.getDtaStartRange());
        this.plot.setDtaEndRange(peptide.getDtaEndRange());

        //System.out.println("indepen" + xmlData.toString());
        drawPlot(xmlData.toString());

    }

    public void generateInDepData(ChroPeptide peptide) {
        if (null == peptide) {
            peptide = currentPeptide;
        }

        if (null == currentPeptide) {
            return;
        }

        StringBuffer xmlData = new StringBuffer();
        xmlData.append(chromXMLHeader(peptide));

        List l = peptide.getDataList();
        int startRange = 0;
        int endRange = 0;

        startRange = Integer.parseInt(peptide.getStartRange());
        endRange = Integer.parseInt(peptide.getEndRange());
        int pepLength = ((ChroData) l.get(0)).getResidueLength();

        long[] samArr = new long[l.size()];
        long[] refArr = new long[samArr.length];

        //xmlData.append("<size width=\"590\" height=\"390\"/>");
        StringBuffer sampleData = new StringBuffer();
        StringBuffer refData = new StringBuffer();
        StringBuffer ticks = new StringBuffer();
        sampleData.append("<dataset name=\"sample\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        refData.append("<dataset name=\"reference\" marks=\"none\" connected=\"yes\" stems=\"no\">");

        StringBuffer[] bufferArr = new StringBuffer[pepLength * 4];
        for (int i = 0; i < pepLength; i++) {
            bufferArr[i] = new StringBuffer();
            bufferArr[i].append("<dataset name=\"bSample").append(i).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        }
        for (int i = pepLength; i < pepLength * 2; i++) {
            bufferArr[i] = new StringBuffer();
            bufferArr[i].append("<dataset name=\"ySample").append(i % pepLength).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        }
        for (int i = pepLength * 2; i < pepLength * 3; i++) {
            bufferArr[i] = new StringBuffer();
            bufferArr[i].append("<dataset name=\"bRef").append(i % pepLength).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        }
        for (int i = pepLength * 3; i < pepLength * 4; i++) {
            bufferArr[i] = new StringBuffer();
            bufferArr[i].append("<dataset name=\"yRef").append(i % pepLength).append("\" marks=\"none\" connected=\"yes\" stems=\"no\">");
        }

        ticks.append("<xTicks>");

        int index = 0;
        int startIndex = 0;
        int endIndex = 0;

        for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
            ChroData data = itr.next();

            long bsArr[] = data.getBsIntensity();
            long ysArr[] = data.getYsIntensity();
            long brArr[] = data.getBrIntensity();
            long yrArr[] = data.getYrIntensity();

            for (int i = 0; i < bufferArr.length; i++) {
                bufferArr[i].append("<p x=\"");
                bufferArr[i].append(data.getScanNum());
                bufferArr[i].append("\" y=\"");

//                if(data.getScanNum() == 4731)
//                    System.out.println(bsArr[i%pepLength] + " " + ysArr[i%pepLength] + " " + brArr[i%pepLength] + " " + yrArr[i%pepLength]);
                switch ((int) (i / pepLength)) {
                    case 0:
                        bufferArr[i].append(bsArr[i % pepLength]);
                        break;
                    case 1:
                        bufferArr[i].append(ysArr[i % pepLength]);
                        break;
                    case 2:
                        bufferArr[i].append(brArr[i % pepLength]);
                        break;
                    case 3:
                        bufferArr[i].append(yrArr[i % pepLength]);
                        break;
                    default:
                        break;
                }

                bufferArr[i].append("\"/>");
            }

            samArr[index] = data.getSampleIntensity();
            refArr[index] = data.getRefIntensity();

            //if( 1==cr.getDataDependency() )
            {
                int scanTemp = data.getScanNum();
                if (startRange >= scanTemp) {
                    startIndex = index;
                }
                if (endRange >= scanTemp) {
                    endIndex = index;
                }
            }

            index++;

            //<tick label="90-100" position="90" />
            sampleData.append("<p x=\"");
            sampleData.append(data.getScanNum());
            sampleData.append("\" y=\"");
            sampleData.append(data.getSampleIntensity());
            sampleData.append("\"/>");

            refData.append("<p x=\"");
            refData.append(data.getScanNum());
            refData.append("\" y=\"");
            refData.append(data.getRefIntensity());
            refData.append("\"/>");

            ticks.append("<tick label=\"");
            ticks.append(data.getScanNum());
            ticks.append("\" position=\"");
            ticks.append(data.getScanNum());
            ticks.append("\" />");

        }

        LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

        double slope = reg.getSlope();
        double intercept = reg.getIntercept();

        //for(int i=startIndex;i<=endIndex;i++)
        //    System.out.println("datapoint\t" + samArr[i] + "\t" + refArr[i]);
        corrPlot.setData(refArr, samArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), slope, intercept, reg.getBestShift());
        corrPlot.clear(true);
        corrPlot.repaint();
        correlationPanel.add(corrPlot);
        //proteinRatioDistPanel.add(corrPlot);

        //this.corrCoeffField.setText( CensusHelper.format.format(reg.getCorr()) );
        if (reg.getCorr() < 0) {
            this.rrField.setText("N/A");
            this.areaRatioLogField.setText("N/A");
            this.regressionRatioField.setText("N/A");
            this.areaRatioField.setText("N/A");
        } else {
            this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));
            this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
            this.regressionRatioField.setText(CensusHelper.format.format(slope));
            this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));
        }

        this.shiftField.setText(String.valueOf(reg.getBestShift()));

        ticks.append("</xTicks>");
        sampleData.append("</dataset>");
        refData.append("</dataset>");

        xmlData.append(ticks.toString());
        xmlData.append(sampleData.toString());
        xmlData.append(refData.toString());

        PostOptions options = PostOptions.getInstance();

        if (options.isDisplayFragmentIons()) {
            for (int i = 0; i < bufferArr.length; i++) {
                bufferArr[i].append("</dataset>");
                xmlData.append(bufferArr[i].toString());
            }
        }

        xmlData.append("</plot>");

        //System.out.println(xmlData.toString());
        this.plot.setScanNum(peptide.getScanNum());
        this.plot.setDtaStartRange(peptide.getDtaStartRange());
        this.plot.setDtaEndRange(peptide.getDtaEndRange());

        drawPlot(xmlData.toString());

    }

    private void proteinSimpleTableKeyPressed(java.awt.event.KeyEvent evt) {

        if (evt.getKeyCode() != KeyEvent.VK_DOWN && evt.getKeyCode() != KeyEvent.VK_UP) {
            return;
        }

        JTable table = (JTable) evt.getSource();

        int selectedRow = 0;

        if (evt.getKeyCode() == KeyEvent.VK_DOWN) {
            if (table.getSelectedRow() + 1 >= table.getRowCount()) {
                return;
            }

            selectedRow = table.getSelectedRow() + 1;
            this.currentProtein = proteinList.get(selectedRow);

        } else if (evt.getKeyCode() == KeyEvent.VK_UP) {
            if (table.getSelectedRow() <= 0) {
                return;
            }

            selectedRow = table.getSelectedRow() - 1;

            this.currentProtein = proteinList.get(selectedRow);
        }

        proteinSimpleTableActionPerformed(table, selectedRow);
    }

    private void proteinSimpleTableMouseClicked(java.awt.event.MouseEvent evt) {

        JTable table = (JTable) evt.getSource();

        currentProtein = proteinList.get(table.getSelectedRow());

        proteinSimpleTableActionPerformed(table, table.getSelectedRow());

    }

    private void proteinSimpleTableActionPerformed(JTable table, int selectedRow) {
        if (null != table) {
            for (int i = selectedRow; i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }
        }

        PostOptions options = PostOptions.getInstance();

        //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.
        if (cr.getExpType() > 0) {
            switch (cr.getExpType()) {
                case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
                    updateiTRAQpeptideInfo(table);
                    break;

                case CensusConstants.MSMS_DATA_INDEPENDENT:

                    if (options.isFilterFragmentIons()) {
                        updateMS2PeptideFilterInfo(null, selectedRow);
                    } else {
                        updatePeptideInfo();
                    }

                    break;

                default:
                    break;
            }
        } else if (cr.getQuantLevel() == 1) {
            updatePeptideInfo();
        } else if (cr.getQuantLevel() == 2) {
            if (!cr.isLabeled() && null != cr.getFileList() && cr.getFileList().size() > 1) //labeling free with multiple samples
            {
                updateLabelFreeMS2PeptideInfo(table);
            } else {
                if (options.isFilterFragmentIons()) {
                    updateMS2PeptideFilterInfo(null, selectedRow);
                } else {
                    updatePeptideInfo();
                }

            }
        }
    }

    //label free ms2 multiple samples
    private void updateLabelFreeMS2PeptideInfo() {
        updateLabelFreeMS2PeptideInfo(null);
    }

    //label free ms2 multiple samples
    private void updateLabelFreeMS2PeptideInfo(JTable table) {
        if (null != table) {
            currentProtein = proteinList.get(table.getSelectedRow());

            for (int i = table.getSelectedRow(); i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }
        }

        List<ChroPeptide> peptideList = currentProtein.getPeptideList();
        List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

        int totalCount = 0;

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
            ChroPeptide peptide = pepItr.next();
            totalCount++;

            List l = peptide.getDataList();

            long[] samArr = new long[l.size()];

            int startRange = Integer.parseInt(peptide.getStartRange());
            int endRange = Integer.parseInt(peptide.getEndRange());

            tempPepList.add(peptide);
        }

        //Remove all existing rows
        int rowCount = peptideTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            peptideTableModel.removeRow(0); //note zero here
        }

        this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(0);

        ChroPeptide peptide;

        for (Iterator<ChroPeptide> itr = currentProtein.getPeptideList().iterator(); itr.hasNext();) {
            peptide = itr.next();
            peptideTableModel.addRow(peptide.getLabelFreePeptideData());

        }

        peptideListTable.changeSelection(0, 0, false, false);//.set.setSelectionBackground(new Color(102, 102, 153));

        PostOptions options = PostOptions.getInstance();

        this.peptideListTable.repaint();
        this.selectQuantType();
        //this.peptidePanel.remove( irisPanel );

        peptidePanel.add(this.quanPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 63, 1394, -1));
        this.peptidePanel.invalidate();
        this.peptidePanel.validate();
        this.peptidePanel.repaint();

        this.proteinLabel.setText(this.currentProtein.getDescription());
    }

    private void updatePeptideInfo() {
        updatePeptideInfo(null);
    }

    private void updatePeptideInfo(JTable table) {
        updatePeptideInfo(table, -1, -1);
    }

    private void updatePeptideInfo(JTable table, int startRange, int endRange) {
        if (null != table) {
            currentProtein = proteinList.get(table.getSelectedRow());

            //System.out.println(currentProtein.getLocus() + " " + table.getSelectedRow());
            for (int i = table.getSelectedRow(); i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }
        }

        List<ChroPeptide> peptideList = currentProtein.getPeptideList();
        List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

        int totalCount = 0;

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
            ChroPeptide peptide = pepItr.next();
            totalCount++;

            List l = peptide.getDataList();

            long[] samArr = new long[l.size()];
            long[] refArr = null;

            if (cr.isLabeled()) {
                refArr = new long[samArr.length];
            }
//	    System.out.println("cur protein1==>>" + peptide.getScoreHt() + " " + currentProtein.getLocus() + " " + peptide.getSequence());

            int index = 0;
            int startIndex = 0;
            int endIndex = samArr.length - 1;

            if (startRange < 0 && endRange < 0) {
                startRange = Integer.parseInt(peptide.getStartRange());
                endRange = Integer.parseInt(peptide.getEndRange());
            }

            LinearRegression reg = null;

            for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext();) {
                ChroData data = dataItr.next();

                samArr[index] = data.getSampleIntensity();

                if (cr.isLabeled()) {
                    refArr[index] = data.getRefIntensity();
                }

                int scanTemp = data.getScanNum();

                if (startRange >= scanTemp) {
                    startIndex = index;
                }
                if (endRange >= scanTemp) {
                    endIndex = index;
                }

              //  System.out.println(scanTemp + " " + startIndex + " " + endIndex + " " + startRange + " " + endRange);
                index++;
            }

            double samIntSum = 0;
            double refIntSum = 0;

            //System.out.println(startIndex + " " + endIndex + " " + samArr.length);
            endIndex = (endIndex != 0) ? endIndex : (samArr.length - 1);
            for (int ii = startIndex; ii <= endIndex; ii++) {
                samIntSum += samArr[ii];

                if (cr.isLabeled()) {
                    refIntSum += refArr[ii];
                }
            }

            if (cr.isLabeled()) {
                reg = new LinearRegression(samArr, refArr, startIndex, endIndex, conf.getMaxSpectrumShift());
                double slope = reg.getSlope();
                double intercept = reg.getIntercept();
                slope = Math.exp( slope > 0 ? (Math.log(slope) ) : 0);

                peptide.setSlope(slope);
                peptide.setCorr(reg.getCorr());
                peptide.setRefIntensity(refIntSum);
            }

            //normalize the ratio in case overall ratios shifted
            peptide.setSamIntensity(samIntSum);

            tempPepList.add(peptide);
        }

        double averageRatio = 0;
        int peptideCount = 0;
        double ratioSum = 0;

        for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
            ChroPeptide each = tempItr.next();

            ratioSum += each.getSlope();
            peptideCount++;
        }

        if (peptideCount > 0) {
            averageRatio = ratioSum / peptideCount;
        }

        double devSum = 0;

        for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
            ChroPeptide each = tempItr.next();

            //if(noFilter || !each.isFilterOut())
            {
                double dev = each.getSlope() - averageRatio;
                devSum += dev * dev;
            }
        }

        //Remove all existing rows
        int rowCount = peptideTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            peptideTableModel.removeRow(0); //note zero here
        }
        this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(0);

        ChroPeptide peptide;

        StringBuffer pepDistSb = new StringBuffer(getPlotHeader());

        pepDistPlot.clearData();

        for (Iterator<ChroPeptide> itr = currentProtein.getPeptideList().iterator(); itr.hasNext();) {
            peptide = itr.next();
            peptideTableModel.addRow(conf.isLabeling() ? peptide.getPeptideData() : peptide.getLabelFreePeptideData());

            if (peptide.getSlope() > 0) {
                this.pepDistPlot.addData(Math.log(peptide.getSlope()), peptide.getCorr() * peptide.getCorr());
                pepDistSb.append("<p x=\"").append(Math.log(peptide.getSlope())).append("\" y=\"").append(peptide.getCorr() * peptide.getCorr()).append("\" />");
            }
        }

        pepDistSb.append("</dataset>");
        pepDistSb.append("</plot>");

        drawPepDist(pepDistSb.toString());

        peptideListTable.changeSelection(0, 0, false, false);//.set.setSelectionBackground(new Color(102, 102, 153));

        PostOptions options = PostOptions.getInstance();

        this.peptideListTable.repaint();
        this.selectQuantType();
        //this.peptidePanel.remove( irisPanel );

        peptidePanel.add(this.quanPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 63, 1394, -1));
        this.peptidePanel.invalidate();
        this.peptidePanel.validate();
        this.peptidePanel.repaint();

        this.proteinLabel.setText(this.currentProtein.getDescription());
    }

    public String getPlotHeader() {
        StringBuffer sb = new StringBuffer();
        sb.append("<?xml version=\"1.0\" standalone=\"no\"?>");
        sb.append("<plot>");
        sb.append("<title>Peptide Ratio Distribution</title>");
        sb.append("<xLabel>Ln(Ratio)</xLabel>");
        sb.append("<yLabel>Determinant Factor(RxR)</yLabel>");
        sb.append("<noGrid/>");
        sb.append("<size width=\"550\" height=\"370\"/>");
        sb.append("<dataset name=\"dmso\" marks=\"none\" connected=\"no\" stems=\"no\">");
        sb.append("<p x=\"-1.0\" y=\"1.0\" />");
        sb.append("<p x=\"0.0\" y=\"1.0\" />");
        sb.append("<p x=\"1.0\" y=\"1.0\" />");

        return sb.toString();
    }

    private void updateMS2PeptideFilterInfo(JTable table, int selectedRow) {
        if (null != table) {
            currentProtein = proteinList.get(selectedRow);

            for (int i = selectedRow; i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }
        }

        List<ChroPeptide> peptideList = currentProtein.getPeptideList();
        List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

        int totalCount = 0;

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
            ChroPeptide peptide = pepItr.next();
            totalCount++;

            DataIndepModel model = this.calcFragIons(peptide);

            int startIndex = model.startIndex;
            int endIndex = model.endIndex;
            long[] samArr = model.samArr;
            long[] refArr = model.refArr;
            int[] scanNumArr = model.scanNumArr;

            FragIonList ionList = CalcUtil.getBestFragIons(model.bsTempArr, model.ysTempArr, model.brTempArr, model.yrTempArr, model.startIndex, model.endIndex, conf.getMaxSpectrumShift());

            int listSize = ionList.size();
            int tempIndex = 0;

            for (Iterator<FragIon> itr = ionList.iterator(); itr.hasNext();) {
                FragIon ion = itr.next();

                long[] tempSArr = ion.getSArr();
                long[] tempRArr = ion.getRArr();

                String ionName = ion.isBion() ? "b" : "y";
                ionName += ion.getIndex();

                for (int i = 0; i < tempSArr.length; i++) {
                    samArr[i] += tempSArr[i];
                    refArr[i] += tempRArr[i];

                }

                if (tempIndex == ionList.getBestIndex()) {
                    break;
                }

                tempIndex++;
            }

            // display low quality of fragions
            for (int i = tempIndex + 1; i < listSize; i++) {
                FragIon ion = ionList.get(i);

                long[] tempSArr = ion.getSArr();
                long[] tempRArr = ion.getRArr();

                String ionName = ion.isBion() ? "b" : "y";
                ionName += ion.getIndex();

            }

            LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

            double slope = reg.getSlope();
            double intercept = reg.getIntercept();

            double samIntSum = 0;
            double refIntSum = 0;

            for (int ii = startIndex; ii <= endIndex; ii++) {
                samIntSum += samArr[ii];
                refIntSum += refArr[ii];
            }

            //normalize the ratio in case overall ratios shifted
            slope = Math.exp( slope > 0 ? (Math.log(slope) ) : 0);

            peptide.setSlope(slope);
            peptide.setCorr(reg.getCorr());
            peptide.setSamIntensity(samIntSum);
            peptide.setRefIntensity(refIntSum);
            //peptide.setSnRatio(snRatio);

            tempPepList.add(peptide);

            //tempPepCount++;
            //}
        }
	//if(!noFilter && tempPepList.size()>3 && pValueSelect)
        //    edu.scripps.pms.stats.GrubbsTest.filter(tempPepList, pValue);

        double averageRatio = 0;
        int peptideCount = 0;
        double ratioSum = 0;

        for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
            ChroPeptide each = tempItr.next();

              //  if(!noFilter && each.isFilterOut())
            //      continue;
            ratioSum += each.getSlope();
            peptideCount++;
            //quantifiedCount++;
        }

        if (peptideCount > 0) {
            averageRatio = ratioSum / peptideCount;
        }

            //proteinSb.append( averageRatio>0?CensusHelper.format.format(averageRatio):"" );
        double devSum = 0;

//            StringBuffer pepSb = new StringBuffer();
        for (Iterator<ChroPeptide> tempItr = tempPepList.iterator(); tempItr.hasNext();) {
            ChroPeptide each = tempItr.next();

//                     if(isUniquePeptide && !each.isUnique())
//                         continue;
            //if(noFilter || !each.isFilterOut())
            {
                double dev = each.getSlope() - averageRatio;
                devSum += dev * dev;

            }
        }

        //Remove all existing rows
        int rowCount = peptideTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            peptideTableModel.removeRow(0); //note zero here
        }

        this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(0);

        ChroPeptide peptide;

        StringBuffer pepDistSb = new StringBuffer();
        pepDistSb.append("<?xml version=\"1.0\" standalone=\"no\"?>");
        pepDistSb.append("<plot>");
        pepDistSb.append("<title>Peptide Ratio Distribution</title>");
        pepDistSb.append("<xLabel>Ln(Ratio)</xLabel>");
        pepDistSb.append("<yLabel>Determinant Factor(RxR)</yLabel>");
        pepDistSb.append("<noGrid/>");
        pepDistSb.append("<size width=\"550\" height=\"370\"/>");

        /*
         pepDistSb.append("<xTicks/>");
         pepDistSb.append("<tick label=\"-3.0\" position=\"-3.0\" />");
         pepDistSb.append("<tick label=\"-2.0\" position=\"-2.0\" />");
         pepDistSb.append("<tick label=\"-1.0\" position=\"-1.0\" />");
         pepDistSb.append("<tick label=\"0.0\" position=\"0.0\" />");
         pepDistSb.append("<tick label=\"1.0\" position=\"1.0\" />");
         pepDistSb.append("<tick label=\"2.0\" position=\"2.0\" />");
         pepDistSb.append("<tick label=\"3.0\" position=\"3.0\" />");
         pepDistSb.append("</xTicks>");

         /*
         pepDistSb.append("<yTicks>");
         pepDistSb.append("<tick label=\"50\" position=\"50\" />");
         pepDistSb.append("<tick label=\"60\" position=\"60\" />");
         pepDistSb.append("<tick label=\"70\" position=\"70\" />");
         pepDistSb.append("</yTicks>");
         */
        pepDistSb.append("<dataset name=\"dmso\" marks=\"none\" connected=\"no\" stems=\"no\">");
        //pepDistSb.append("<p x=\"-3.0\" y=\"0.0\" />");
        //pepDistSb.append("<p x=\"-2.0\" y=\"1.0\" />");
        pepDistSb.append("<p x=\"-1.0\" y=\"1.0\" />");
        pepDistSb.append("<p x=\"0.0\" y=\"1.0\" />");
        pepDistSb.append("<p x=\"1.0\" y=\"1.0\" />");
        //pepDistSb.append("<p x=\"2.0\" y=\"1.0\" />");
        //pepDistSb.append("<p x=\"3.0\" y=\"1.0\" />");

        //sample xml
        //String str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><plot><title>Paths</title><xLabel>Reference Chromatogram</xLabel><yLabel>Target Chromatogram</yLabel><noGrid /><size width=\"750\" height=\"750\" /><xTicks /><tick label=\"1678\" position=\"1678\" /><tick label=\"1311\" position=\"1311\" /><tick label=\"932\" position=\"932\" /><tick label=\"572\" position=\"572\" /><tick label=\"213\" position=\"213\" /><yTicks />" +
        //            "<tick label=\"1766\" position=\"1766\" /><tick label=\"1378\" position=\"1378\" /><tick label=\"1007\" position=\"1007\" /><tick label=\"617\" position=\"617\" /><tick label=\"226\" position=\"226\" /><dataset name=\"/home/rpark/001/nonlabel/IsotopeFree/SetB/rep1/092106_IsotopeFree_B-01_itms.ms1\" sample=\"SetB\" marks=\"none\" connected=\"yes\" stems=\"no\"><p x=\"2059\" y=\"2136\" /><p x=\"2058\" y=\"2135\" /><p x=\"2057\" y=\"2134\" /><p x=\"2056\" y=\"2133\" />" +
        //            "<p x=\"2055\" y=\"2133\" /><p x=\"2054\" y=\"2132\" /></dataset></plot>";
        pepDistPlot.clearData();

        for (Iterator<ChroPeptide> itr = currentProtein.getPeptideList().iterator(); itr.hasNext();) {
            peptide = itr.next();
            peptideTableModel.addRow(conf.isLabeling() ? peptide.getPeptideData() : peptide.getLabelFreePeptideData());

            if (peptide.getSlope() > 0) {

                this.pepDistPlot.addData(Math.log(peptide.getSlope()), peptide.getCorr() * peptide.getCorr());

                pepDistSb.append("<p x=\"").append(Math.log(peptide.getSlope())).append("\" y=\"").append(peptide.getCorr() * peptide.getCorr()).append("\" />");

        //        pepDistSb.append("<p x=\"").append( Math.log(peptide.getSlope()) ).append("\" y=\"").append( peptide.getCorr()*peptide.getCorr() ).append("\"/>");
            }
        }

        pepDistSb.append("</dataset>");
        pepDistSb.append("</plot>");

        drawPepDist(pepDistSb.toString());

        //if(currentProtein.getPeptideList().size()!=tempPepList.size())
        //    System.out.println("x");
        //else
        //    System.out.println("*");
        peptideListTable.changeSelection(0, 0, false, false);//.set.setSelectionBackground(new Color(102, 102, 153));

        PostOptions options = PostOptions.getInstance();

        this.peptideListTable.repaint();

        if (!this.isLabeled() && this.quantLevel == 1) {
            generateNonLabelData(currentPeptide);
        } else if (this.isDataDependent) {
            generateDepData(currentPeptide);
        } else if (options.isFilterFragmentIons()) {
            generateInDepFragData(currentPeptide);
        } else {
            generateInDepData(currentPeptide);
        }

        //this.peptidePanel.remove( irisPanel );
        peptidePanel.add(this.quanPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 63, 1394, -1));
        this.peptidePanel.invalidate();
        this.peptidePanel.validate();
        this.peptidePanel.repaint();

        this.proteinLabel.setText(this.currentProtein.getDescription());

    }

    private void updateiTRAQpeptideInfo() {
        updateiTRAQpeptideInfo(null, -1, -1);
    }

    private void updateiTRAQpeptideInfo(JTable table) {
        updateiTRAQpeptideInfo(table, -1, -1);
    }

    public void updateiTRAQpeptideInfo(JTable table, int startRange, int endRange) {
        if (null != table) {
            currentProtein = proteinList.get(table.getSelectedRow());

            for (int i = table.getSelectedRow(); i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }
        }

        List<ChroPeptide> peptideList = currentProtein.getPeptideList();
        List<ChroPeptide> tempPepList = new Vector<ChroPeptide>();

        int totalCount = 0;

        for (Iterator<ChroPeptide> pepItr = peptideList.iterator(); pepItr.hasNext();) {
            ChroPeptide peptide = pepItr.next();
            totalCount++;

            List<ChroData> l = peptide.getDataList();

            ChroData fdata = l.get(0);

            long[][] dataArr = new long[l.size()][fdata.getIntensityArr().length];

            int index = 0;
            int startIndex = 0;
            int endIndex = dataArr.length - 1;

            if (null != peptide.getStartRange() && null != peptide.getEndRange() && startRange < 0 && endRange < 0) {
                startRange = Integer.parseInt(peptide.getStartRange());
                endRange = Integer.parseInt(peptide.getEndRange());
            }

            for (Iterator<ChroData> dataItr = l.iterator(); dataItr.hasNext();) {
                ChroData data = dataItr.next();

                dataArr[index] = data.getIntensityArr();

                int scanTemp = data.getScanNum();
                if (startRange >= scanTemp) {
                    startIndex = index;
                }
                if (endRange >= scanTemp) {
                    endIndex = index;
                }

                index++;
            }

            switch (cr.getExpType()) {
                case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                    break;

                case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
                    long[] intSumArr = new long[dataArr[0].length];

                    endIndex = (endIndex != 0) ? endIndex : (dataArr.length - 1);
                    for (int ii = startIndex; ii <= endIndex; ii++) {
                        for (int j = 0; j < dataArr[ii].length; j++) {
                            intSumArr[j] += dataArr[ii][j];
                        }
                    }

                    break;

                default:
                    break;
            }

            tempPepList.add(peptide);

        }

        double averageRatio = 0;
        int peptideCount = 0;
        double ratioSum = 0;

        //Remove all existing rows
        int rowCount = peptideTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            peptideTableModel.removeRow(0); //note zero here
        }
        this.currentPeptide = (ChroPeptide) currentProtein.getPeptideList().get(0);

        ChroPeptide peptide;

        StringBuffer pepDistSb = new StringBuffer(getPlotHeader());

        pepDistPlot.clearData();

        for (Iterator<ChroPeptide> itr = currentProtein.getPeptideList().iterator(); itr.hasNext();) {
            peptide = itr.next();
            peptideTableModel.addRow(peptide.getPeptideData());

            if (peptide.getSlope() > 0) {
                this.pepDistPlot.addData(Math.log(peptide.getSlope()), peptide.getCorr() * peptide.getCorr());
                pepDistSb.append("<p x=\"").append(Math.log(peptide.getSlope())).append("\" y=\"").append(peptide.getCorr() * peptide.getCorr()).append("\" />");
            }
        }

        pepDistSb.append("</dataset>");
        pepDistSb.append("</plot>");

        drawPepDist(pepDistSb.toString());

        peptideListTable.changeSelection(0, 0, false, false);//.set.setSelectionBackground(new Color(102, 102, 153));

        PostOptions options = PostOptions.getInstance();

        this.peptideListTable.repaint();
        this.selectQuantType();
        //this.peptidePanel.remove( irisPanel );

        peptidePanel.add(this.quanPanel, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 63, 1394, -1));
        this.peptidePanel.invalidate();
        this.peptidePanel.validate();
        this.peptidePanel.repaint();

        this.proteinLabel.setText(this.currentProtein.getDescription());

    }

    private void proteinTableMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_proteinTableMouseClicked
        // TODO add your handling code here:

        //double click
        if (evt.getClickCount() != 2) {
            return;
        }

        PostOptions options = PostOptions.getInstance();

        tabbedPanel.setSelectedIndex(1);
        tabbedPanel.setEnabledAt(1, true);
        JTable table = (JTable) evt.getSource();

        if (this.isChroFile) {

            //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.
            if (cr.getExpType() > 0) {
                switch (cr.getExpType()) {
                    case CensusConstants.MSMS_DATA_INDEPENDENT:
                        updatePeptideInfo(table);
                        break;
                    case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                    case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
                        updateiTRAQpeptideInfo(table);
                        break;

                    default:
                        break;
                }
            } else if (cr.getQuantLevel() == 1) {
                updatePeptideInfo(table);

            } else if (cr.getQuantLevel() == 2) {
                if (!cr.isLabeled() && null != cr.getFileList() && cr.getFileList().size() > 1) //labeling free with multiple samples
                {
                    updateLabelFreeMS2PeptideInfo(table);
                } else //typical MRM experiment
                {
                    updatePeptideInfo(table);
                }
            }

        } else {
            for (int i = table.getSelectedRow(); i < proteinList.size(); i++) {
                currentProtein = proteinList.get(i);

//                System.out.println( currentProtein.getPeptideList().size() + " " + proteinList.size() + " " + table.getSelectedColumn());
                if (currentProtein.getPeptideList().size() > 0) {
                    break;
                }
            }

            qualPanel.setCurrentProtein(currentProtein);
        }
    }//GEN-LAST:event_proteinTableMouseClicked

    private void exitItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exitItemActionPerformed
        System.exit(0);

        // TODO add your handling code here:
    }//GEN-LAST:event_exitItemActionPerformed

    public void dummyOpenChroFile() {
        Configuration conf = Configuration.getInstance();
        conf.setSimpleIndexGenerator(true);

        this.openChroFile("/home/rpark/001/census_chro.xml");

    }

    private void OpenChroFileActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_OpenChroFileActionPerformed
        // TODO add your handling code here:
        Configuration conf = Configuration.getInstance();
        conf.setSimpleIndexGenerator(true);

        JFileChooser choose = new JFileChooser();
        choose.setMultiSelectionEnabled(false);
        choose.setDialogTitle("Select Chro File");
        choose.addChoosableFileFilter(new SimpleFileNameFilter("xml", "Chro File"));

        if (currentDirectory != null && !"".equals(currentDirectory)) {
            choose.setCurrentDirectory(new File(currentDirectory));
        }

        int returnVal = choose.showOpenDialog(chroPanel);
        this.chroFile = choose.getSelectedFile();

        if (null == chroFile || returnVal == choose.CANCEL_OPTION) {
            return;
        }

        currentDirectory = chroFile.getAbsolutePath();
        currentDirectory = currentDirectory.substring(0, currentDirectory.lastIndexOf(File.separator));

        filePathLabel.setText(chroFile.getAbsolutePath());

        cleanupProteinTableModel();

        this.openChroFile(chroFile.getAbsolutePath());

    }//GEN-LAST:event_OpenChroFileActionPerformed

    private void cleanupProteinTableModel() {
        chromatogramPanel.removeAll();

        int rowCount = proteinTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            proteinTableModel.removeRow(0); //note zero here
        }
        rowCount = peptideTableModel.getRowCount();
        for (int i = 0; i < rowCount; i++) {
            peptideTableModel.removeRow(0); //note zero here
        }
        rowCount = this.proteinSimpleTableModel.getRowCount();

        for (int i = 0; i < rowCount; i++) {
            proteinSimpleTableModel.removeRow(0); //note zero here
        }
    }

    private void resetBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_resetBtnActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_resetBtnActionPerformed

    public ArrayList getProteinList() {
        return proteinList;
    }


    private void extractBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_extractBtnActionPerformed
        // TODO add your handling code here:
        Configuration conf = Configuration.getInstance();
        conf.setSimpleIndexGenerator(false);

        final String isFull = quantModeRadioGrp.getSelection().getActionCommand();

        final ChroProgressDialog chroProgress = new ChroProgressDialog(this, false);
        chroProgress.setLocationRelativeTo(this);
        chroProgress.setResizable(false);
        chroProgress.setVisible(true);

        final long start = System.currentTimeMillis();

        Thread t = new Thread() {
            boolean isSuccessful = false;

            public void run() {
                try {

                    ChroGenerator chro = new ChroGenerator(
                            chroProgress.getProgressBar(),
                            //chroProgress.getProgressText(),
                            null,
                            Integer.parseInt(scanBefore.getText()),
                            Integer.parseInt(scanAfter.getText()),
                            isotopeFileField.getText().trim(),
                            massTolerance
                    );
                    if (isFull.equals("f")) //Full Scan
                    {
                        chro.createFullscanXmlChro();
                    } else //msms scan
                    {
                        chro.createMsmsXmlChro(chroProgress);
                    }

                    isSuccessful = true;

                } catch (IOException e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(chroPanel, "Failed to generate a chro file: " + e, "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);
                    isSuccessful = false;
                } catch (Exception e) {
                    e.printStackTrace();
                    JOptionPane.showMessageDialog(chroPanel, "Failed to generate a chro file: " + e, "Failed to generate a chro file", JOptionPane.ERROR_MESSAGE);

                    isSuccessful = false;
                }

                SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        chroProgress.setVisible(false);
                        chroProgress.hide();


                        /*
                         if(isSuccessful)
                         {
                         JOptionPane.showMessageDialog(
                         chroPanel,
                         "Chro file is successfully created",
                         "Chro file Creation",
                         JOptionPane.PLAIN_MESSAGE);
                         }
                         */
                    }
                });
            }
        };

        try {
            t.start();
        } catch (Exception e) {
            t = null;

            //chroProgress.setVisible(false);
        }

    }//GEN-LAST:event_extractBtnActionPerformed

    private void fileSelectBtnActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_fileSelectBtnActionPerformed
        // TODO add your handling code here:
        JFileChooser choose = new JFileChooser();
        choose.setMultiSelectionEnabled(false);
        choose.setDialogTitle("Select Isotope File");
        choose.addChoosableFileFilter(new SimpleFileNameFilter("ini", "Isotope File"));

        if (currentDirectory != null && !"".equals(currentDirectory)) {
            choose.setCurrentDirectory(new File(currentDirectory));
        }

        int returnVal = choose.showOpenDialog(chroPanel);
        File file = choose.getSelectedFile();

        isotopeFileField.setText(file.getAbsolutePath());
    }//GEN-LAST:event_fileSelectBtnActionPerformed

    public void changePeakArea(int startPeak, int endPeak) {

        int oldStartRange = Integer.parseInt(this.currentPeptide.getStartRange());
        int oldEndRange = Integer.parseInt(this.currentPeptide.getEndRange());

        if (startPeak == oldStartRange && endPeak == oldEndRange) {
            return;
        }

        this.currentPeptide.setStartRange(String.valueOf(startPeak));
        this.currentPeptide.setEndRange(String.valueOf(endPeak));

        PostOptions options = PostOptions.getInstance();

        //we will use exp type only in the future.  No more many if else.. checking quantlevel or labeled check.
        if (cr.getExpType() > 0) {
            switch (cr.getExpType()) {
                case CensusConstants.MSMS_SPECIFIC_SINGLE_MASS: //iTRAQ
                case CensusConstants.MSMS_SPECIFIC_MULTIPLE_MASS: //iTRAQ
                    updateiTRAQpeptideInfo();
                    break;

                case CensusConstants.MSMS_DATA_INDEPENDENT:
                    runDataIndependentPeakArea(options);
                    break;

                default:
                    break;
            }
        } //non label
        else if (!this.isLabeled()) {
            int startRange = startPeak;
            int endRange = endPeak;

            long[] intensitySumArr = new long[cr.getFileList().size()];

            int index = 0;
            int startIndex = 0;
            int endIndex = 0;

            int rowCount = nonlabelTableModel.getRowCount();

            for (int i = 0; i < rowCount; i++) {
                nonlabelTableModel.removeRow(0); //note zero here
            }
            rowCount = nonlabelSummaryTableModel.getRowCount();
            for (int i = 0; i < rowCount; i++) {
                nonlabelSummaryTableModel.removeRow(0); //note zero here
            }
            if (this.quantLevel == 1) {
                List<ChroNonLabelData> l = currentPeptide.getDataList();

                for (Iterator<ChroNonLabelData> dataItr = l.iterator(); dataItr.hasNext();) {
                    ChroNonLabelData eachData = dataItr.next();

                    int[] scanNumArr = eachData.getScanNumArr();
                    long[] intenArr = eachData.getIntensityArr();

                    for (int i = 0; i < scanNumArr.length; i++) {
                        if (scanNumArr[0] >= startRange && scanNumArr[0] <= endRange) {
                            intensitySumArr[i] += intenArr[i];
                        }
                    }

                    int scanTemp = eachData.getScanNum();

                    if (startRange >= scanTemp) {
                        startIndex = index;
                    }
                    if (endRange >= scanTemp) {
                        endIndex = index;
                    }

                    index++;
                }

            } else if (cr.getQuantLevel() == 2) {
                List<ChroNonLabelMSMSData> l = currentPeptide.getDataList();

                for (Iterator<ChroNonLabelMSMSData> dataItr = l.iterator(); dataItr.hasNext();) {
                    ChroNonLabelMSMSData eachData = dataItr.next();

                    int[] scanNumArr = eachData.getScanArr();
                    long[] intenArr = eachData.getTotalIntArr();

                    for (int i = 0; i < scanNumArr.length; i++) {
                        if (scanNumArr[0] >= startRange && scanNumArr[0] <= endRange) {
                            intensitySumArr[i] += intenArr[i];
                        }
                    }

                    int scanTemp = eachData.getScanNum();

                    if (startRange >= scanTemp) {
                        startIndex = index;
                    }
                    if (endRange >= scanTemp) {
                        endIndex = index;
                    }

                    index++;
                }
//                updateLabelFreeMS2PeptideInfo();
            }

            for (int i = 0; i < intensitySumArr.length; i++) {

                String curDir = cr.getFileList().get(i);

                String sampleName = cr.getSampleName(curDir);
                Vector vec = new Vector();
                vec.add(sampleName);
                vec.add(curDir.substring(curDir.lastIndexOf(File.separator) + 1));

                //we don't know if file path is generated from linux or window.
                if (curDir.startsWith("/")) {
                    vec.add(curDir.substring(0, curDir.lastIndexOf("/")));
                } else {
                    vec.add(curDir.substring(0, curDir.lastIndexOf("\\")));
                }

                vec.add(CensusHelper.scientificFormat.format(intensitySumArr[i]));

                this.nonlabelTableModel.addRow(vec);

            }
        } else if (this.isLabeled() && this.quantLevel == 1) //else if(this.isDataDependent)
        {
            List l = currentPeptide.getDataList();

            int startRange = startPeak;
            int endRange = endPeak;

            long[] samArr = new long[l.size()];
            long[] refArr = new long[samArr.length];

            int index = 0;
            int startIndex = 0;
            int endIndex = 0;

            for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
                ChroData data = itr.next();

                samArr[index] = data.getSampleIntensity();
                refArr[index] = data.getRefIntensity();

                int scanTemp = data.getScanNum();

                if (startRange >= scanTemp) {
                    startIndex = index;
                }
                if (endRange >= scanTemp) {
                    endIndex = index;
                }

                index++;
            }

            LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

            double slope = reg.getSlope();
            double intercept = reg.getIntercept();

            corrPlot.setData(refArr, samArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), slope, intercept, reg.getBestShift());
            corrPlot.clear(true);
            corrPlot.repaint();
            correlationPanel.add(corrPlot);
            //proteinRatioDistPanel.add(corrPlot);

            if (reg.getCorr() < 0) {
                this.rrField.setText("N/A");
                this.areaRatioLogField.setText("N/A");
                this.regressionRatioField.setText("N/A");
                this.areaRatioField.setText("N/A");
            } else {
                this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));
                this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
                this.regressionRatioField.setText(CensusHelper.format.format(slope));
                this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));
            }

            this.shiftField.setText(String.valueOf(reg.getBestShift()));
        } else//data independent // if(options.isFilterFragmentIons())
        {
            runDataIndependentPeakArea(options);
        }

    }

    /**
     * @param args the command line arguments
     *
     * public static void main(String args[]) {
     * java.awt.EventQueue.invokeLater(new Runnable() { public void run() { new
     * RelExFileFilter().setVisible(true); } }); }
     */
    /*
     *new javax.swing.table.DefaultTableModel(
     new Object [][] {
     {"fdsa", null, null, null, null, null, null, null},
     {null, "fdsa", null, null, null, null, null, null},
     {null, null, null, null, null, null, null, null},
     {null, null, null, null, null, null, null, null}
     },
     new String [] {
     "Locus", "Sequence Count", "Spectrum Count", "Sequence Coverage", "Length", "MolWt", "pI", "Description"
     }
     ) {
     Class[] types = new Class [] {
     java.lang.String.class, java.lang.String.class, java.lang.String.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class, java.lang.Object.class
     };
     boolean[] canEdit = new boolean [] {
     false, true, true, true, true, true, true, true
     };

     public Class getColumnClass(int columnIndex) {
     return types [columnIndex];
     }

     public boolean isCellEditable(int rowIndex, int columnIndex) {
     return canEdit [columnIndex];
     }
     }*/
    private void runDataIndependentPeakArea(PostOptions options) {
        if (options.isFilterFragmentIons()) {
            fragIonPanel.removeAll();

            List l = currentPeptide.getDataList();

            int startRange = 0;
            int endRange = 0;

            startRange = Integer.parseInt(currentPeptide.getStartRange());
            endRange = Integer.parseInt(currentPeptide.getEndRange());

	    //StringBuffer sampleData = new StringBuffer();
            //StringBuffer refData = new StringBuffer();
            //StringBuffer ticks = new StringBuffer();
            int pepLength = ((ChroData) l.get(0)).getResidueLength();

            long[] samArr = new long[l.size()];
            long[] refArr = new long[samArr.length];
            int[] scanNumArr = new int[samArr.length];

            int index = 0;
            int startIndex = 0;
            int endIndex = 0;

            long[][] bsTempArr = new long[pepLength][samArr.length];
            long[][] ysTempArr = new long[pepLength][samArr.length];
            long[][] brTempArr = new long[pepLength][samArr.length];
            long[][] yrTempArr = new long[pepLength][samArr.length];

            for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
                ChroData data = itr.next();

                long bsArr[] = data.getBsIntensity();
                long ysArr[] = data.getYsIntensity();
                long brArr[] = data.getBrIntensity();
                long yrArr[] = data.getYrIntensity();

                for (int i = 0; i < bsArr.length; i++) {
                    bsTempArr[i][index] = bsArr[i];
                    ysTempArr[i][index] = ysArr[i];
                    brTempArr[i][index] = brArr[i];
                    yrTempArr[i][index] = yrArr[i];
                }

                scanNumArr[index] = data.getScanNum();

                int scanTemp = data.getScanNum();
                if (startRange >= scanTemp) {
                    startIndex = index;
                }
                if (endRange >= scanTemp) {
                    endIndex = index;
                }

                index++;

            }

            FragIonList ionList = CalcUtil.getBestFragIons(bsTempArr, ysTempArr, brTempArr, yrTempArr, startIndex, endIndex, conf.getMaxSpectrumShift());

            int listSize = ionList.size();
            int tempIndex = 0;

            for (Iterator<FragIon> itr = ionList.iterator(); itr.hasNext();) {
                FragIon ion = itr.next();

                long[] tempSArr = ion.getSArr();
                long[] tempRArr = ion.getRArr();

                String ionName = ion.isBion() ? "b" : "y";
                ionName += ion.getIndex();

                FragIonPlot fPlot = new FragIonPlot(ionName, tempSArr, tempRArr, startIndex, (endIndex != 0) ? endIndex : (tempSArr.length - 1), true);
                fragIonPanel.add(fPlot, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 0 + (tempIndex * 60), 270, 60));

                for (int i = 0; i < tempSArr.length; i++) {
                    samArr[i] += tempSArr[i];
                    refArr[i] += tempRArr[i];

                }

                if (tempIndex == ionList.getBestIndex()) {
                    break;
                }

                tempIndex++;
            }

            // display low quality of fragions
            for (int i = tempIndex + 1; i < listSize; i++) {
                FragIon ion = ionList.get(i);

                long[] tempSArr = ion.getSArr();
                long[] tempRArr = ion.getRArr();

                String ionName = ion.isBion() ? "b" : "y";
                ionName += ion.getIndex();

                FragIonPlot fPlot = new FragIonPlot(ionName, tempSArr, tempRArr, startIndex, (endIndex != 0) ? endIndex : (tempSArr.length - 1), false);
                fragIonPanel.add(fPlot, new org.netbeans.lib.awtextra.AbsoluteConstraints(0, 0 + (i * 60), 270, 60));
            }

            this.fragIonPanel.invalidate();
            this.fragIonPanel.validate();
            this.fragIonPanel.repaint();

            /*
             this.sigToNoisePanel.removeAll();

             SigNoisePlot sigPlot = new SigNoisePlot(this.sigToNoisePanel.getWidth(), this.sigToNoisePanel.getHeight(), ionList, tempIndex);
             sigPlot.setBackground(new Color(255, 255, 255));
             this.sigToNoisePanel.add(sigPlot);
             this.sigToNoisePanel.invalidate();
             this.sigToNoisePanel.validate();
             this.sigToNoisePanel.repaint();
             */
            LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

            double slope = reg.getSlope();
            double intercept = reg.getIntercept();

            corrPlot.setData(refArr, samArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), slope, intercept, reg.getBestShift());
            corrPlot.clear(true);
            corrPlot.repaint();
            correlationPanel.add(corrPlot);
            //proteinRatioDistPanel.add(corrPlot);

            if (reg.getCorr() < 0) {
                this.rrField.setText("N/A");
                this.areaRatioLogField.setText("N/A");
                this.regressionRatioField.setText("N/A");
                this.areaRatioField.setText("N/A");
            } else {
                this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));
                this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
                this.regressionRatioField.setText(CensusHelper.format.format(slope));
                this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));
            }

            this.shiftField.setText(String.valueOf(reg.getBestShift()));

            this.plot.setScanNum(currentPeptide.getScanNum());
            this.plot.setDtaStartRange(currentPeptide.getDtaStartRange());
            this.plot.setDtaEndRange(currentPeptide.getDtaEndRange());
        } else {
            List l = currentPeptide.getDataList();
            StringBuffer xmlData = new StringBuffer();

            int startRange = 0;
            int endRange = 0;

            startRange = Integer.parseInt(currentPeptide.getStartRange());
            endRange = Integer.parseInt(currentPeptide.getEndRange());

            int pepLength = ((ChroData) l.get(0)).getResidueLength();

            long[] samArr = new long[l.size()];
            long[] refArr = new long[samArr.length];

            int index = 0;
            int startIndex = 0;
            int endIndex = 0;

            for (Iterator<ChroData> itr = l.iterator(); itr.hasNext();) {
                ChroData data = itr.next();

                long bsArr[] = data.getBsIntensity();
                long ysArr[] = data.getYsIntensity();
                long brArr[] = data.getBrIntensity();
                long yrArr[] = data.getYrIntensity();

                samArr[index] = data.getSampleIntensity();
                refArr[index] = data.getRefIntensity();

                //if( 1==cr.getDataDependency() )
                {
                    int scanTemp = data.getScanNum();
                    if (startRange >= scanTemp) {
                        startIndex = index;
                    }
                    if (endRange >= scanTemp) {
                        endIndex = index;
                    }
                }

                index++;
            }

            LinearRegression reg = new LinearRegression(samArr, refArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), conf.getMaxSpectrumShift());

            double slope = reg.getSlope();
            double intercept = reg.getIntercept();

            corrPlot.setData(refArr, samArr, startIndex, (endIndex != 0) ? endIndex : (samArr.length - 1), slope, intercept, reg.getBestShift());
            corrPlot.clear(true);
            corrPlot.repaint();
            correlationPanel.add(corrPlot);
            //proteinRatioDistPanel.add(corrPlot);

            if (reg.getCorr() < 0) {
                this.rrField.setText("N/A");
                this.areaRatioLogField.setText("N/A");
                this.regressionRatioField.setText("N/A");
                this.areaRatioField.setText("N/A");
            } else {
                this.rrField.setText(CensusHelper.format.format(reg.getCorr() * reg.getCorr()));
                this.areaRatioLogField.setText(CensusHelper.format.format( slope > 0 ? (Math.log(slope) ) : 0));
                this.regressionRatioField.setText(CensusHelper.format.format(slope));
                this.areaRatioField.setText(CensusHelper.format.format(reg.getAreaRatio()));
            }

            this.shiftField.setText(String.valueOf(reg.getBestShift()));
            this.plot.setScanNum(currentPeptide.getScanNum());
            this.plot.setDtaStartRange(currentPeptide.getDtaStartRange());
            this.plot.setDtaEndRange(currentPeptide.getDtaEndRange());

        }

    }

    private DefaultTableModel proteinTableModel
            = new DefaultTableModel(CensusConstants.PROTEIN_COLUMNS, 0) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

    private DefaultTableModel proteinSimpleTableModel
            = new DefaultTableModel(CensusConstants.PROTEIN_SIMPLE_COLUMNS, 0) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

    Configuration conf = Configuration.getInstance();

    private DefaultTableModel peptideTableModel
            = new DefaultTableModel(conf.isLabeling() ? CensusConstants.PEPTIDE_COLUMNS : CensusConstants.PEPTIDE_LABEL_FREE_COLUMNS, 0) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

    public void setCurrentDirectory(File file) {
        this.currentDirectory = file.getAbsolutePath();
    }

    public void setCurrentDirectory(String curDirect) {
        this.currentDirectory = curDirect;
    }

    public String getCurrentDirectory() {
        return this.currentDirectory;
    }


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JMenuItem alignSpectra;
    private javax.swing.JTextField areaRatioField;
    private javax.swing.JTextField areaRatioLogField;
    private javax.swing.JLabel chorNoteLabel1;
    private javax.swing.JLabel chroNoteLabel2;
    private javax.swing.JPanel chroPanel;
    private javax.swing.JPanel chromatogramPanel;
    private javax.swing.JMenuItem confItem;
    private javax.swing.JPanel correlationPanel;
    private javax.swing.JMenuItem exitItem;
    private javax.swing.JButton export;
    private javax.swing.JButton extractBtn;
    private javax.swing.JPanel extractPanel;
    private javax.swing.ButtonGroup extractRadioGrp;
    private javax.swing.JMenu fileMenu;
    private javax.swing.JLabel filePathLabel;
    private javax.swing.JButton fileSelectBtn;
    private javax.swing.JButton filterBtn;
    private javax.swing.JPanel fragIonPanel;
    private javax.swing.JScrollPane fragIonScrollPanel;
    private javax.swing.JRadioButton fullMassScan;
    private javax.swing.JMenu helpMenu;
    private javax.swing.JTextField isotopeFileField;
    private javax.swing.JLabel isotopeLabel;
    private javax.swing.JCheckBox jCheckBox1;
    private javax.swing.JCheckBox jCheckBox2;
    private javax.swing.JCheckBox jCheckBox3;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel4;
    private javax.swing.JRadioButton jRadioButton1;
    private javax.swing.JRadioButton jRadioButton2;
    private javax.swing.JSeparator jSeparator1;
    private javax.swing.JSplitPane jSplitPane;
    private javax.swing.JTextField jTextField1;
    private javax.swing.JTextField jTextField5;
    private javax.swing.JTextField jTextField6;
    private javax.swing.JToolBar jToolBar1;
    private javax.swing.JPanel mainPanel;
    private javax.swing.JLabel measuredRatioLabel;
    private javax.swing.JLabel measuredRatioLabel1;
    private javax.swing.JMenuItem mergeItem;
    private javax.swing.JMenuItem mrmCsv;
    private javax.swing.JRadioButton msmsScan;
    private javax.swing.JButton open;
    private javax.swing.JMenuItem openItem;
    private javax.swing.JMenuItem openSpectra;
    private javax.swing.JMenuItem optionItem;
    private javax.swing.JPanel paramPanel;
    private javax.swing.JTabbedPane pepTabbedPanel;
    private javax.swing.JScrollPane peptideList;
    private javax.swing.JPanel peptideListBox;
    private javax.swing.JPanel peptidePanel;
    private javax.swing.JPanel proteinInfoPanel;
    private javax.swing.JLabel proteinLabel;
    private javax.swing.JScrollPane proteinListPanel;
    private javax.swing.JPanel proteinPanel;
    private javax.swing.JScrollPane proteinSimplePane;
    private javax.swing.JTable proteinTable;
    private javax.swing.JTextField pvalueField;
    private javax.swing.JLabel pvalueLabel;
    private javax.swing.JPanel quanPanel;
    private javax.swing.ButtonGroup quantModeRadioGrp;
    private javax.swing.JLabel regLnLabel;
    private javax.swing.JTextField regScoreField;
    private javax.swing.JLabel regScoreLabel;
    private javax.swing.JTextField regressionRatioField;
    private javax.swing.JMenuItem reportItem;
    private javax.swing.JButton resetBtn;
    private javax.swing.JTextField rrField;
    private javax.swing.JLabel rrLabel;
    private javax.swing.JMenuItem runItem;
    private javax.swing.JMenu runMenu;
    private javax.swing.JMenuItem runNonLabel;
    private javax.swing.JButton save;
    private javax.swing.JTextField scanAfter;
    private javax.swing.JTextField scanBefore;
    private javax.swing.JButton searchBtn;
    private javax.swing.JTextField searchField;
    private javax.swing.JTextField shiftField;
    private javax.swing.JLabel shiftLabel;
    private javax.swing.JTabbedPane tabbedPanel;
    private javax.swing.JMenu toolMenu;
    private javax.swing.JMenuItem versionItem;
    private javax.swing.JLabel welcomeLabel;
    // End of variables declaration//GEN-END:variables
    //private List<Element> proteinList;

    private javax.swing.JTable proteinSimpleTable;
    private javax.swing.JTable peptideListTable;

    public javax.swing.JPanel getProteinPanel() {
        return this.proteinPanel;

    }

    private ArrayList<ChroProtein> proteinList;

    private BaseChroPlot plot;
    private PeptideDistPlot pepDistPlot;

    private CorrelationPlot corrPlot;
    //private CorrelationPlot corrPlot;

    private StringBuffer xmlData;
    private ChroProtein currentProtein;
    private PlotBoxMLParser chroXmlParser;
    private PlotMLParser peptideDistParser;

    private final int INTEGRATION_WINDOW = 15;
    private final double PEAK_THRESHOLD = 0.15;
    private final int PEAK_SCAN_BEFORE = 4;
    private final int PEAK_SCAN_AFTER = 4;
    private ChroXmlReader cr;

    private String currentDirectory = "";
    //private int MAX_SHIFT=0;
//    private DecimalFormat format = new DecimalFormat("0.000");
//    private DecimalFormat twoDigitFormat = new DecimalFormat("0.00");
//    private NumberFormat scientificFormat = new DecimalFormat("0.###E0");

    private float massTolerance = 0.3f;

    //protected PostOptions options = new PostOptions();
    protected ChroPeptide currentPeptide = null;
    protected boolean isDataDependent;

    /**
     * ************** experiment Type ********
     */
    protected int experimentType = 0;

    protected int quantLevel = 1; //default it is for ms1
    protected boolean labeled = true; //default it is labeled

    protected IrisPanel irisPanel = null;
    protected Hashtable ht = null; //indexed Hashtable
    protected File chroFile = null;
    protected File specFile = null;

    protected int searchIndex = 0;

    protected DTASelectFilterReader dtaReader;
    protected boolean isChroFile = true;

    protected QualificationPanel qualPanel = new QualificationPanel();

    protected JPanel nonLabelPanel = new JPanel();
    protected JPanel nonLabelSummaryPanel = new JPanel();

    //protected JTextArea nonlabelText = new JTextArea();
    protected JScrollPane nonlabelScrollPane = new JScrollPane();
    protected JScrollPane nonlabelSummaryScrollPane = new JScrollPane();

    protected DefaultTableModel nonlabelTableModel
            = new DefaultTableModel(CensusConstants.NONLABEL_COLUMNS, 0) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

    protected DefaultTableModel nonlabelSummaryTableModel
            = new DefaultTableModel(CensusConstants.NONLABEL_SUMMARY_COLUMNS, 0) {
                public boolean isCellEditable(int row, int column) {
                    return false;
                }
            };

    protected JTable nonlabelTable = new JTable();
    protected JTable nonlabelSummaryTable = new JTable();
    protected JPanel proteinRatioDistPanel = new JPanel();

    /**
     * ************** Configuration parameters **************
     */
    private java.net.URI fileUri = null;

    public URI getFileUri() {
        return fileUri;
    }

    private static final double MAX_RATIO_LOG2 = 6.64385619;
}
