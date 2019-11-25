package gui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Timer;
import java.util.TimerTask;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Style;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyledDocument;

import logic.DTASelect;

/**
 * ProgramProgressBarPanel
 * Singleton Class
 */
public class ProgramProgressBarPanel extends JFrame {
	/**
	 * Singleton instance
	 */
	private static ProgramProgressBarPanel instance = null;
	
    public final static int ONE_SECOND   = 1000;

    private JLabel infoLabel             = null;
    private JProgressBar progressBar     = null;
    
    private JTextPane programsOutputArea = new JTextPane();
    private StyledDocument doc           = null;
    private Style styleNormal            = null;
    private Style styleError             = null;
    private Style styleWarning           = null;
    
    private JScrollPane scrollPane       = null;
    private String applicationPath       = ((System.getProperty("java.class.path")).split(";"))[0];
	private static String dir            = "images\\";
	private String icon                  = "icon.jpg";
	private int status                   = 0; // status that shows progress in a numerical value (= current percent)
	private int maxStatus                = -1; // -1 = default value, can never be reached by program
	                                           // maximal status (= 100%)
	private Timer timer                  = null;
	private boolean somethingHasChanged  = false;
	
	private ProgramProgressBarPanel() {
		super("DTASelect in Progess ...");
	
		// set icon:
		setIconImage(Toolkit.getDefaultToolkit().getImage(applicationPath + "\\" + dir + icon));
	
		// to make the window closeable:
		this.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				try {
					System.exit(0);
				} catch (Exception e2) {
					e2.printStackTrace();
				}
			}
		});
		
		this.setLayout(new BorderLayout()); // set layout
		
		addComponents();
		
		this.setSize(600, 500);
		this.setResizable(false);
	    this.setLocation((Toolkit.getDefaultToolkit().getScreenSize().width-
                          this.getSize().width) / 2,
                         (Toolkit.getDefaultToolkit().getScreenSize().height-
                          this.getSize().height) / 2);

	    createStyles();
	}
	
	private void createStyles() {
		// Create a style object and then set the style attributes
        styleNormal = doc.addStyle("normal", null);
        StyleConstants.setForeground(styleNormal, Color.BLACK);
        
        styleError = doc.addStyle("error", null);
        StyleConstants.setBold(styleError, true);
        StyleConstants.setForeground(styleError, Color.RED);
        
        styleWarning = doc.addStyle("warning", null);
        StyleConstants.setBold(styleWarning, true);
        StyleConstants.setForeground(styleWarning, Color.ORANGE);
	}

	public static ProgramProgressBarPanel getInstance() {
		if (instance == null) {
			instance = new ProgramProgressBarPanel();
		}
		return instance;
	}
	
	private void addComponents() {
		infoLabel = new JLabel("  Progress ...");
		
		progressBar = new JProgressBar();
		progressBar.setMinimum(0);
		progressBar.setValue(0); // initial value: 0% completed!
		progressBar.setStringPainted(true);
		
		JPanel northPane = new JPanel();
		northPane.setLayout(new GridLayout(0, 1));
		northPane.add(infoLabel);
		northPane.add(progressBar);
		
		add(northPane, BorderLayout.NORTH);
		
		programsOutputArea = new JTextPane(); // textarea for programs output
        doc = (StyledDocument)programsOutputArea.getDocument();
		scrollPane = new JScrollPane(programsOutputArea,
		                                         JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
		                                         JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		programsOutputArea.setEditable(false);
		programsOutputArea.setAutoscrolls(true);
		scrollPane.setAutoscrolls(true);
		add(scrollPane, BorderLayout.CENTER);
		
		// Create time to update ProgramProgessBarPanel
	    int delay = 0;   // delay for 1 sec.
	    int period = 1000;  // repeat every sec.
	    Timer timer = new Timer();
	    
	    timer.scheduleAtFixedRate(new TimerTask() {
	            public void run() {
	            	if (getInstance().isVisible() && somethingHasChanged) // update progressBar and outputArea
	            		doUpdate();
	            }
	        }, delay, period);
	}
	
	public void setVisible(Boolean value) {
		this.setVisible(value);
	}
	
	public void writeOutputToTextArea(String output, int outputCode) {
	    progressBar.paint(progressBar.getGraphics());
	    infoLabel.paint(infoLabel.getGraphics());
	    
	    switch (outputCode) {
	    	case DTASelect.CORRECT_OUTPUT:
		    	// append  normal output in black font color to document
		        try {
					doc.insertString(doc.getLength(), output, styleNormal);
				} catch (BadLocationException e) {
					e.printStackTrace();
				}
	    		break;
	    	case DTASelect.ERROR_OUTPUT:
		    	// write error message in red, to symbolize that an error has occurred	         
		        try {
					doc.insertString(doc.getLength(), output, styleError);
				} catch (BadLocationException e) {
					e.printStackTrace();
				} 
				// and stop application
				showProgramStoppedMessage();
	    		break;
	    	case DTASelect.WARNING_OUTPUT:
		    	// write warning message in orange, to symbolize that an warning is given         
		        try {
					doc.insertString(doc.getLength(), "WARNING: " + output, styleWarning);
				} catch (BadLocationException e) {
					e.printStackTrace();
				} 
	    		break;
	    	default:
	    		try {
					throw new Exception("wrong outputCode passed!");
				} catch (Exception e) {
					e.printStackTrace();
				}
	    }

	    programsOutputArea.setCaretPosition(programsOutputArea.getDocument().getLength()/*-1*/);

	    // necessary to automatically scroll down, when a text is added to the textarea
	    try {
	    	if (programsOutputArea.isVisible())
	    		programsOutputArea.scrollRectToVisible(programsOutputArea.modelToView(programsOutputArea.getDocument().getLength()));
		} catch (BadLocationException e1) {
			e1.printStackTrace();
		}
		
	    programsOutputArea.paint(programsOutputArea.getGraphics());
	    scrollPane.paint(scrollPane.getGraphics());
	    
	    somethingHasChanged = true;
	}
	
	public void showCompletedMessage(String message) {
		// if DTASearch is finished show finished message:
			JOptionPane.showMessageDialog(this, message,
				                          ".: Completed :.", JOptionPane.INFORMATION_MESSAGE);
	}
	
	public void showProgramStoppedMessage() {
		// if DTASearch is finished show finished message:
		JOptionPane.showMessageDialog(this, "Program stopped\nbecause an ERROR (see red output) has occurred!",
			                          ".: Program Stopped :.", JOptionPane.ERROR_MESSAGE);
		
		// after showing Information come to termination:
		System.exit(0);
	}
	
	public void increaseStatus() {
		++status; // increase status for one unit
		if (status == maxStatus) {
			progressBar.setString("100% COMPLETED");
		}
	    progressBar.setValue(status);

		somethingHasChanged = true;
	}
	
	public void setMaxStatus(int maxStatus) {
		this.maxStatus = maxStatus;
		progressBar.setMaximum(maxStatus);
	}
	
	public void setToMaximum() {
		progressBar.setValue(maxStatus);
		somethingHasChanged = true;
	}
	
	private void doUpdate() {
	    progressBar.paint(progressBar.getGraphics());
	    infoLabel.paint(infoLabel.getGraphics());
	    
	    programsOutputArea.setCaretPosition(programsOutputArea.getDocument().getLength()/*-1*/);
	    
	    // necessary to automatically scroll down, when a text is added to the textarea
	    try {
	    	if (programsOutputArea.isVisible())
	    		programsOutputArea.scrollRectToVisible(programsOutputArea.modelToView(programsOutputArea.getDocument().getLength()));
		} catch (BadLocationException e1) {
			e1.printStackTrace();
		}
	    
	    programsOutputArea.paint(programsOutputArea.getGraphics());
	    scrollPane.paint(scrollPane.getGraphics());
	    
	    // set somethingHasChanged variable back since everything is on acutal state
	    somethingHasChanged = false;
	}

	/*==================================================================*/
	/* Getters and Setters                                              */
	/*==================================================================*/
	public JTextPane getProgramsOutputArea() {
		return programsOutputArea;
	}

}