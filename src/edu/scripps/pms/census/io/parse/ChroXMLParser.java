/*
 * ChroXMLParser.java
 *
 * Created on July 24, 2005, 9:16 PM
 */

package edu.scripps.pms.census.io.parse;

import ptolemy.plot.plotml.PlotMLParser;
//import ptolemy.plot.PlotBox;
import ptolemy.plot.Plot;
import com.microstar.xml.XmlException;

/**
 *
 * @author rpark
 */
public class ChroXMLParser extends PlotMLParser {
    
    /** Creates a new instance of ChroXMLParser */
    
    public ChroXMLParser(Plot plot) {
        super(plot);
    }
    
    /** Start an element.
     *  This is called at the beginning of each XML
     *  element.  By the time it is called, all of the attributes
     *  for the element will already have been reported using the
     *  attribute() method.  Unrecognized elements are ignored.
     *  @param elementName The element type name.
     *  @exception XmlException If the element produces an error
     *   in constructing the model.
     */
    public void startElement(String elementName) throws XmlException {
        
        try {
            // NOTE: The elements are alphabetical below...

            if (elementName.equals("barGraph")) {
                System.out.println("bar graph...");
                String widthSpec = (String)_attributes.get("width");
                String offsetSpec = (String)_attributes.get("offset");
                // NOTE: If only one of these is given, then the other
                // is ignored.
                if (widthSpec == null || offsetSpec == null) {
                    ((Plot)_plot).setBars(true);
                } else {
                    double width = (Double.valueOf(widthSpec)).doubleValue();
                    double offset = (Double.valueOf(offsetSpec)).doubleValue();
                    ((Plot)_plot).setBars(width, offset);
                }

            } else if (elementName.equals("dataset")) {
                _currentDataset++;
                _currentPointCount = 0.0;

                String connected = (String)_attributes.get("connected");
                if (connected != null) {
                    if (connected.equals("no")) {
                        ((Plot)_plot).setConnected(false, _currentDataset);
                    } else {
                        ((Plot)_plot).setConnected(true, _currentDataset);
                    }
                }

                String marks = (String)_attributes.get("marks");
                if (marks != null) {
                    ((Plot)_plot).setMarksStyle(marks, _currentDataset);
                }

                String name = (String)_attributes.get("name");
                if (name != null) {
                    ((Plot)_plot).addLegend(_currentDataset, name);
                }

                String stems = (String)_attributes.get("stems");
                if (stems != null) {
                    if (stems.equals("yes")) {
                        ((Plot)_plot).setImpulses(true, _currentDataset);
                    } else {
                        ((Plot)_plot).setImpulses(false, _currentDataset);
                    }
                }

            }else if (elementName.equals("p")) {                                
                _addPoint(true, elementName);
            }else {

                super.startElement(elementName);
            }
        } catch (Exception ex) {
            if (ex instanceof XmlException) {
                throw (XmlException)ex;
            } else {
                // FIXME: Temporary for debugging.
                System.err.println(ex.toString());
                ex.printStackTrace();
                String msg = "XML element \"" + elementName
                    + "\" triggers exception:\n  " + ex.toString();
                throw new XmlException(msg,
                        _currentExternalEntity(),
                        _parser.getLineNumber(),
                        _parser.getColumnNumber());
            }
        }
        
        // NOTE: if super is called, this gets done twice.
        // Any way to avoid it?
        _attributes.clear();
    }

    public void endElement(String elementName) throws Exception {
        super.endElement(elementName);

        if (elementName.equals("dataset")) {
            // Reset the default, in case it was changed for this dataset.
            ((Plot)_plot).setConnected(_connected);
        }
    }

}
