/*
 * TextAreaRenderer.java
 *
 * Created on August 25, 2005, 2:09 PM
 */

package edu.scripps.pms.census;

/**
 *
 * @author rpark
 */

import java.awt.Component;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.TableCellRenderer;

public class TextAreaRenderer extends JTextArea
implements TableCellRenderer {

    public TextAreaRenderer() {
        setLineWrap(true);
        setWrapStyleWord(true);
    }

    public Component getTableCellRendererComponent(JTable jTable,
            Object obj, boolean isSelected, boolean hasFocus, int row,
            int column) {
        setText((String)obj);
        return this;
    }
}