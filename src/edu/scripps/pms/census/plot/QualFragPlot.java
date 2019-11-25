/*
 * BaseChroPlot.java
 *
 * Created on July 25, 2005, 10:55 PM
 */

package edu.scripps.pms.census.plot;

import ptolemy.plot.Plot;
import java.awt.*;
import java.util.*;

/**
 *
 * @author rpark
 */
public class QualFragPlot extends Plot {
    
    private double xcorr;
    private int dtaStartRange;
    private int dtaEndRange;
        
    /** Creates a new instance of BaseChroPlot */
    public QualFragPlot() {
        super();
    }

    public double getXcorr() {
        return xcorr;
    }

    public void setXcorr(double xcorr) {
        this.xcorr = xcorr;
    }
    
    protected synchronized void _drawPlot(
        Graphics graphics, boolean clearfirst, Rectangle drawRect) {
        // Ignore if there is no graphics object to draw on.
        if (graphics == null) return;

        graphics.setPaintMode();

        /* NOTE: The following seems to be unnecessary with Swing...
           if (clearfirst) {
           // NOTE: calling clearRect() here permits the background
           // color to show through, but it messes up printing.
           // Printing results in black-on-black title and axis labels.
           graphics.setColor(_background);
           graphics.drawRect(0, 0, drawRect.width, drawRect.height);
           graphics.setColor(Color.black);
           }
        */
                
        // If an error message has been set, display it and return.
        if (_errorMsg != null) {
            int fheight = _labelFontMetrics.getHeight() + 2;
            int msgy = fheight;
            graphics.setColor(Color.black);
            for (int i = 0; i < _errorMsg.length; i++) {
                graphics.drawString(_errorMsg[i], 10, msgy);
                msgy += fheight;
                System.err.println(_errorMsg[i]);
            }
            return;
        }

        // Make sure we have an x and y range
        if (!_xRangeGiven) {
            if (_xBottom > _xTop) {
                // have nothing to go on.
                _setXRange(0, 0);
            } else {
                _setXRange(_xBottom, _xTop);
            }
        }
        if (!_yRangeGiven) {
            if (_yBottom > _yTop) {
                // have nothing to go on.
                _setYRange(0, 0);
            } else {
                _setYRange(_yBottom, _yTop);
            }
        }

        // Vertical space for title, if appropriate.
        // NOTE: We assume a one-line title.
        int titley = 0;
        int titlefontheight = _titleFontMetrics.getHeight();

        if (_title == null) {
            // NOTE: If the _title is null, then set it to the empty
            // string to solve the problem where the fill button overlaps
            // the legend if there is no title.  The fix here would
            // be to modify the legend printing text so that it takes
            // into account the case where there is no title by offsetting
            // just enough for the button.
            _title = "";
        }
        if (_title != null || _yExp != 0) {
            titley = titlefontheight + _topPadding;
        }

        // Number of vertical tick marks depends on the height of the font
        // for labeling ticks and the height of the window.
        Font previousFont = graphics.getFont();
        graphics.setFont(_labelFont);
        graphics.setColor(_foreground);	// foreground color not set here  --Rob.
        int labelheight = _labelFontMetrics.getHeight();
        int halflabelheight = labelheight/2;
        int maxXTickLength = 0;


        // Draw scaling annotation for x axis.
        // NOTE: 5 pixel padding on bottom.
        int ySPos = drawRect.height - 20;
        int xSPos = drawRect.width - _rightPadding;

        if (_xlog)
            _xExp = (int)Math.floor(_xtickMin);

        if (_xExp != 0 && _xticks == null) {

            String superscript = Integer.toString(_xExp);
            xSPos -= _superscriptFontMetrics.stringWidth(superscript);
            graphics.setFont(_superscriptFont);
            if (!_xlog) {

                //graphics.drawString(superscript, xSPos,
                //                    ySPos - halflabelheight);

                int xlocation = drawRect.width - 100;
                graphics.drawString(superscript, xlocation,
                                  ySPos - halflabelheight);

                xSPos -= _labelFontMetrics.stringWidth("x10");
                xlocation -= _labelFontMetrics.stringWidth("x10");
                graphics.setFont(_labelFont);
                //graphics.drawString("x10", xSPos, ySPos);
                graphics.drawString("x10", xlocation, ySPos);
            }

            //////////Robin come back here to fix padding - original + 5, not 25

            // NOTE: 5 pixel padding on bottom
            _bottomPadding = (3 * labelheight)/2 + 25;
        }

        // NOTE: 5 pixel padding on the bottom.
        if (_xlabel != null && _bottomPadding < labelheight + 5) {
            _bottomPadding = labelheight + 5;
        }

        // Compute the space needed around the plot, starting with vertical.
        // NOTE: padding of 5 pixels below title.
        _uly = titley + 5;
        // NOTE: 3 pixels above bottom labels.
        _lry = drawRect.height-labelheight-_bottomPadding-3;
        int height = _lry-_uly;
        _yscale = height/(_yMax - _yMin);
        _ytickscale = height/(_ytickMax - _ytickMin);

        ////////////////// vertical axis

        // Number of y tick marks.
        // NOTE: subjective spacing factor.
        int ny = 2 + height/(labelheight+10);
        // Compute y increment.
        double yStep = _roundUp((_ytickMax-_ytickMin)/(double)ny);

        // Compute y starting point so it is a multiple of yStep.
        double yStart = yStep*Math.ceil(_ytickMin/yStep);

        // NOTE: Following disables first tick.  Not a good idea?
        // if (yStart == _ytickMin) yStart += yStep;

        // Define the strings that will label the y axis.
        // Meanwhile, find the width of the widest label.
        // The labels are quantized so that they don't have excess resolution.
        int widesty = 0;

        // These do not get used unless ticks are automatic, but the
        // compiler is not smart enough to allow us to reference them
        // in two distinct conditional clauses unless they are
        // allocated outside the clauses.
        String ylabels[] = new String[ny];
        int ylabwidth[] = new int[ny];

        int ind = 0;
        if (_yticks == null) {
            Vector ygrid = null;
            if (_ylog) {
                ygrid = _gridInit(yStart, yStep, true, null);
            }

            // automatic ticks
            // First, figure out how many digits after the decimal point
            // will be used.
            int numfracdigits = _numFracDigits(yStep);

            // NOTE: Test cases kept in case they are needed again.
            // System.out.println("0.1 with 3 digits: " + _formatNum(0.1, 3));
            // System.out.println("0.0995 with 3 digits: " +
            //                    _formatNum(0.0995, 3));
            // System.out.println("0.9995 with 3 digits: " +
            //                    _formatNum(0.9995, 3));
            // System.out.println("1.9995 with 0 digits: " +
            //                    _formatNum(1.9995, 0));
            // System.out.println("1 with 3 digits: " + _formatNum(1, 3));
            // System.out.println("10 with 0 digits: " + _formatNum(10, 0));
            // System.out.println("997 with 3 digits: " + _formatNum(997, 3));
            // System.out.println("0.005 needs: " + _numFracDigits(0.005));
            // System.out.println("1 needs: " + _numFracDigits(1));
            // System.out.println("999 needs: " + _numFracDigits(999));
            // System.out.println("999.0001 needs: "+_numFracDigits(999.0001));
            // System.out.println("0.005 integer digits: " +
            //                    _numIntDigits(0.005));
            // System.out.println("1 integer digits: " + _numIntDigits(1));
            // System.out.println("999 integer digits: " + _numIntDigits(999));
            // System.out.println("-999.0001 integer digits: " +
            //                    _numIntDigits(999.0001));

            double yTmpStart = yStart;
            if (_ylog)
                yTmpStart = _gridStep(ygrid, yStart, yStep, _ylog);

            for (double ypos = yTmpStart; ypos <= _ytickMax;
                 ypos = _gridStep(ygrid, ypos, yStep, _ylog)) {
                // Prevent out of bounds exceptions
                if (ind >= ny) break;
                String yticklabel;
                if (_ylog) {
                    yticklabel = _formatLogNum(ypos, numfracdigits);
                } else {
                    yticklabel = _formatNum(ypos, numfracdigits);
                    //yticklabel = _formatNum(ypos * Math.pow(10,_yExp), numfracdigits);
                }
                ylabels[ind] = yticklabel;
                int lw = _labelFontMetrics.stringWidth(yticklabel);
                ylabwidth[ind++] = lw;
                if (lw > widesty) {widesty = lw;}
            }

        } else {
            // explicitly specified ticks
            Enumeration nl = _yticklabels.elements();
            while (nl.hasMoreElements()) {
                String label = (String) nl.nextElement();
                int lw = _labelFontMetrics.stringWidth(label);
                if (lw > widesty) {widesty = lw;}
            }
        }

        
//        horizontalSpacing(graphics);
        
        // Next we do the horizontal spacing.
        if (_ylabel != null) {
            _ulx = widesty + _labelFontMetrics.stringWidth("W") + _leftPadding;
        } else {
            _ulx = widesty + _leftPadding;
        }
        int legendwidth = _drawLegend(graphics,
                drawRect.width-_rightPadding, _uly);
        _lrx = drawRect.width-legendwidth -_rightPadding;
        int width = _lrx-_ulx;
        _xscale = width/(_xMax - _xMin);

        _xtickscale = width/(_xtickMax - _xtickMin);
        
        
        // Background for the plotting rectangle.
        // Always use a white background because the dataset colors
        // were designed for a white background.
        graphics.setColor(Color.white);
        graphics.fillRect(_ulx, _uly, width, height);

        graphics.setColor(_foreground);
        graphics.drawRect(_ulx, _uly, width, height);

        // NOTE: subjective tick length.
        int tickLength = 5;
        int xCoord1 = _ulx+tickLength;
        int xCoord2 = _lrx-tickLength;

        
        _drawPeak(graphics);
        
        //vertical axis
        if (_yticks == null) {
            // auto-ticks
            Vector ygrid = null;
            double yTmpStart = yStart;
            if (_ylog) {
                ygrid = _gridInit(yStart, yStep, true, null);
                yTmpStart = _gridStep(ygrid, yStart, yStep, _ylog);
                ny = ind;
            }
            ind = 0;
            // Set to false if we don't need the exponent
            boolean needExponent = _ylog;
            for (double ypos = yTmpStart; ypos <= _ytickMax;
                 ypos = _gridStep(ygrid, ypos, yStep, _ylog)) {
                // Prevent out of bounds exceptions
                if (ind >= ny) break;
                int yCoord1 = _lry - (int)((ypos-_ytickMin)*_ytickscale);
                // The lowest label is shifted up slightly to avoid
                // colliding with x labels.
                int offset = 0;
                if (ind > 0 &&  ! _ylog) offset = halflabelheight;
                graphics.drawLine(_ulx, yCoord1, xCoord1, yCoord1);
                graphics.drawLine(_lrx, yCoord1, xCoord2, yCoord1);
                if (_grid && yCoord1 != _uly && yCoord1 != _lry) {
                    graphics.setColor(Color.lightGray);
                    graphics.drawLine(xCoord1, yCoord1, xCoord2, yCoord1);
                    graphics.setColor(_foreground);
                }
                // Check to see if any of the labels printed contain
                // the exponent.  If we don't see an exponent, then print it.
                if (_ylog && ylabels[ind].indexOf('e') != -1 )
                    needExponent = false;

                // NOTE: 4 pixel spacing between axis and labels.
                graphics.drawString(ylabels[ind],
                        _ulx-ylabwidth[ind++]-4, yCoord1+offset);
            }

            if (_ylog) {
                // Draw in grid lines that don't have labels.
                Vector unlabeledgrid  = _gridInit(yStart, yStep, false, ygrid);
                if (unlabeledgrid.size() > 0) {
                    // If the step is greater than 1, clamp it to 1 so that
                    // we draw the unlabeled grid lines for each
                    //integer interval.
                    double tmpStep = (yStep > 1.0)? 1.0 : yStep;

                    for (double ypos = _gridStep(unlabeledgrid , yStart,
                            tmpStep, _ylog);
                         ypos <= _ytickMax;
                         ypos = _gridStep(unlabeledgrid, ypos,
                                 tmpStep, _ylog)) {
                        int yCoord1 = _lry -
                            (int)((ypos-_ytickMin)*_ytickscale);
                        if (_grid && yCoord1 != _uly && yCoord1 != _lry) {
                            graphics.setColor(Color.lightGray);
                            graphics.drawLine(_ulx+1, yCoord1,
                                    _lrx-1, yCoord1);
                            graphics.setColor(_foreground);
                        }
                    }
                }

                if (needExponent) {
                    // We zoomed in, so we need the exponent
                    _yExp = (int)Math.floor(yTmpStart);
                } else {
                    _yExp = 0;
                }
            }

            // Draw scaling annotation for y axis.
            if (_yExp != 0) {
                graphics.drawString("x10", 2, titley*2-20);
                graphics.setFont(_superscriptFont);
                graphics.drawString(Integer.toString(_yExp),
                        _labelFontMetrics.stringWidth("x10") + 2,
                        titley*2-halflabelheight-20);
                graphics.setFont(_labelFont);
            }
        } else {
            // ticks have been explicitly specified
            Enumeration nt = _yticks.elements();
            Enumeration nl = _yticklabels.elements();

            while (nl.hasMoreElements()) {
                String label = (String) nl.nextElement();
                double ypos = ((Double)(nt.nextElement())).doubleValue();
                if (ypos > _yMax || ypos < _yMin) continue;
                int yCoord1 = _lry - (int)((ypos-_yMin)*_yscale);
                int offset = 0;
                if (ypos < _lry - labelheight) offset = halflabelheight;
                graphics.drawLine(_ulx, yCoord1, xCoord1, yCoord1);
                graphics.drawLine(_lrx, yCoord1, xCoord2, yCoord1);
                if (_grid && yCoord1 != _uly && yCoord1 != _lry) {
                    graphics.setColor(Color.lightGray);
                    graphics.drawLine(xCoord1, yCoord1, xCoord2, yCoord1);
                    graphics.setColor(_foreground);
                }
                // NOTE: 3 pixel spacing between axis and labels.
                graphics.drawString(label,
                        _ulx - _labelFontMetrics.stringWidth(label) - 3,
                        yCoord1+offset);
            }
        }

        //////////////////// horizontal axis
        
        int yCoord1 = _uly+tickLength;
        int yCoord2 = _lry-tickLength;
        int charwidth = _labelFontMetrics.stringWidth("8");
        if (_xticks == null) {
            // auto-ticks

            // Number of x tick marks.
            // Need to start with a guess and converge on a solution here.
            int nx = 10;
            double xStep = 0.0;
            int numfracdigits = 0;

            // Limit to 10 iterations
            int count = 0;
            while (count++ <= 10) {
                xStep = _roundUp((_xtickMax-_xtickMin)/(double)nx);
                // Compute the width of a label for this xStep
                numfracdigits = _numFracDigits(xStep);
                // Number of integer digits is the maximum of two endpoints
                int intdigits = _numIntDigits(_xtickMax);
                int inttemp = _numIntDigits(_xtickMin);
                if (intdigits < inttemp) {
                    intdigits = inttemp;
                }
                // Allow two extra digits (decimal point and sign).
                int maxlabelwidth = charwidth *
                    (numfracdigits + 2 + intdigits);
                // Compute new estimate of number of ticks.
                int savenx = nx;
                // NOTE: 10 additional pixels between labels.
                // NOTE: Try to ensure at least two tick marks.
                nx = 2 + width/(maxlabelwidth+10);
                if (nx - savenx <= 1 || savenx - nx <= 1) break;
            }

            xStep = _roundUp((_xtickMax-_xtickMin)/(double)nx);
            numfracdigits = _numFracDigits(xStep);

            // Compute x starting point so it is a multiple of xStep.
            double xStart = xStep*Math.ceil(_xtickMin/xStep);
System.out.println("=================22" + "\t" + xStart + "\t" + _xtickMin);
            // NOTE: Following disables first tick.  Not a good idea?
            // if (xStart == _xMin) xStart += xStep;

            Vector xgrid = null;
            double xTmpStart = xStart;
  //          System.out.println("=================22" + "\t" + _xlog);
            
            // Set to false if we don't need the exponent
            boolean needExponent = _xlog;

            // Label the x axis.  The labels are quantized so that
            // they don't have excess resolution.
            for (double xpos = xTmpStart;
                 xpos <= _xtickMax;
                 xpos = _gridStep(xgrid, xpos, xStep, _xlog)) {

                String xticklabel;

                if (_xlog) {
                    xticklabel = _formatLogNum(xpos, numfracdigits);
                    if (xticklabel.indexOf('e') != -1 )
                        needExponent = false;
                } else {
                    xticklabel = _formatNum(xpos, numfracdigits);
                    //System.out.println("=================" + numfracdigits + "\t" + xpos + "\t" + xTmpStart);
                    
                    //xticklabel = _formatNum(xpos * Math.pow(10,_xExp), numfracdigits);
                }

                xCoord1 = _ulx + (int)((xpos-_xtickMin)*_xtickscale);
                graphics.drawLine(xCoord1, _uly, xCoord1, yCoord1);
                graphics.drawLine(xCoord1, _lry, xCoord1, yCoord2);
                if (_grid && xCoord1 != _ulx && xCoord1 != _lrx) {
                    graphics.setColor(Color.lightGray);
                    graphics.drawLine(xCoord1, yCoord1, xCoord1, yCoord2);
                    graphics.setColor(_foreground);
                }
                int labxpos = xCoord1 -
                    _labelFontMetrics.stringWidth(xticklabel)/2;
                // NOTE: 3 pixel spacing between axis and labels.
                graphics.drawString(xticklabel, labxpos,
                        _lry + 3 + labelheight);
            }

        } else {
            // ticks have been explicitly specified
            Enumeration nt = _xticks.elements();
            Enumeration nl = _xticklabels.elements();

            double preLength = 0.0;
            //int interval = _xticks.size()/10;
            
            int count=0;
            Enumeration tempEnu = _xticks.elements();
            
            while(tempEnu.hasMoreElements())
            {
                double xpos = ((Double)(tempEnu.nextElement())).doubleValue();
                
                if(xpos>_xMin && xpos<_xMax)
                    count++;
            }
            
            int interval = count/10;
            count=0;
            
            
            while (nl.hasMoreElements()) {
                count++;
                String label = (String) nl.nextElement();
                double xpos = ((Double)(nt.nextElement())).doubleValue();
                                
                // If xpos is out of range, ignore.
                if (xpos > _xMax || xpos < _xMin) continue;

                // Find the center position of the label.
                xCoord1 = _ulx + (int)((xpos-_xMin)*_xscale);

                // Find  the start position of x label.
                //int labxpos = xCoord1 - _labelFontMetrics.stringWidth(label)/2;
                int labxpos = xCoord1;
                int labypos = _lry + 3 + labelheight + _labelFontMetrics.stringWidth(label)/2;

                // calculate the length of the label
                preLength = xCoord1
                    + _labelFontMetrics.stringWidth(label)/2 + 10;

                // Draw the label.
                // NOTE: 3 pixel spacing between axis and labels.
                //                    graphics.drawString(label,

                //////////////////XTicks here to be fixed..
                int temp = _labelFontMetrics.stringWidth(label)/3;

                if (maxXTickLength < temp)
                    maxXTickLength = temp;

                Font font = new Font("",Font.PLAIN,9);

                if(interval!=0 && (count%interval != 0))
                        continue;

                Graphics2D g2d = (Graphics2D)graphics;
                g2d.setFont(font);
                // NOTE: Fudge factor so label doesn't touch axis labels.
                g2d.rotate(Math.toRadians(-90), labxpos, labypos);
                g2d.drawString(label, labxpos+5, labypos);
                g2d.rotate(Math.toRadians(90), labxpos, labypos);
                graphics.setFont(previousFont);

                //        labxpos, _lry + 3 + labelheight);

                // Draw the label mark on the axis
                graphics.drawLine(xCoord1, _uly, xCoord1, yCoord1);
                graphics.drawLine(xCoord1, _lry, xCoord1, yCoord2);

                // Draw the grid line
                if (_grid && xCoord1 != _ulx && xCoord1 != _lrx) {
                    graphics.setColor(Color.lightGray);
                    graphics.drawLine(xCoord1, yCoord1, xCoord1, yCoord2);
                    graphics.setColor(_foreground);
                }
            }
        }

        //////////////////// Draw title and axis labels now.

        // Center the title and X label over the plotting region, not
        // the window.
        graphics.setColor(_foreground);

        if (_title != null) {
            graphics.setFont(_titleFont);
            int titlex = _ulx +
                (width - _titleFontMetrics.stringWidth(_title))/2;
            graphics.drawString(_title, titlex, titley + 5);
        }

        graphics.setFont(new Font("",Font.PLAIN,9));
        if (_xlabel != null) {
            int labelx = _ulx +
                (width - _labelFontMetrics.stringWidth(_xlabel))/2;
            //            graphics.drawString(_xlabel, labelx, ySPos);
            graphics.drawString(_xlabel, labelx, ySPos+20);
        }

        int charcenter = 2 + _labelFontMetrics.stringWidth("W")/2;
        if (_ylabel != null) {
            int yl = _ylabel.length();
            if (graphics instanceof Graphics2D) {
                int starty = _uly
                        + (_lry-_uly)/2
                        + _labelFontMetrics.stringWidth(_ylabel)/2
                        - charwidth;
                Graphics2D g2d = (Graphics2D)graphics;
                // NOTE: Fudge factor so label doesn't touch axis labels.
                int startx = charcenter + halflabelheight - 2;
                g2d.rotate(Math.toRadians(-90), startx, starty);
                g2d.drawString(_ylabel, startx, starty);
                g2d.rotate(Math.toRadians(90), startx, starty);
            } else {
                // Not graphics 2D, no support for rotation.
                // Vertical label is fairly complex to draw.
                int starty = _uly
                        + (_lry-_uly)/2 - yl*halflabelheight + labelheight;
                for (int i = 0; i < yl; i++) {
                    String nchar = _ylabel.substring(i, i+1);
                    int cwidth = _labelFontMetrics.stringWidth(nchar);
                    graphics.drawString(nchar, charcenter - cwidth/2, starty);
                    starty += labelheight;
                }
            }
        }

        graphics.setFont(previousFont);
        
        

        // Plot the points in reverse order so that the first colors
        // appear on top.
        for (int dataset = _points.size() - 1; dataset >= 0 ; dataset--) {
            Vector data = (Vector)_points.elementAt(dataset);
            for (int pointnum = 0; pointnum < data.size(); pointnum++) {
                //System.out.println( ((PlotPoint)data.get(pointnum)).x + "\t" + ((PlotPoint)data.get(pointnum)).y);
                _drawPlotPoint(graphics, dataset, pointnum);
            }
        }
        
        //plot Legend
        graphics.setFont(new Font("Helvetica",Font.BOLD,11));
        graphics.setColor(Color.BLUE);
        String sample = "Sample";
        String ref = "Reference";
        String scanNum = "Scan #";
        
        graphics.drawString(sample, _lrx - _labelFontMetrics.stringWidth(sample + ref + scanNum + 10), _uly-2);
        //graphics.drawLine(_lrx - _labelFontMetrics.stringWidth(ref + 10), _uly-6, _lrx - _labelFontMetrics.stringWidth(ref + 10) + 5, _uly-6);
        graphics.setColor(Color.RED);
        graphics.drawString(ref, _lrx - _labelFontMetrics.stringWidth(ref + scanNum + 5), _uly-2);
        graphics.setColor(Color.GREEN);
        graphics.drawString(scanNum, _lrx - _labelFontMetrics.stringWidth(scanNum), _uly-2);
        
        
        
        
    
        _showing = true;        
    }        

    public int getDtaStartRange() {
        return dtaStartRange;
    }

    public void setDtaStartRange(int dtaStartRange) {
        this.dtaStartRange = dtaStartRange;
    }

    public int getDtaEndRange() {
        return dtaEndRange;
    }

    public void setDtaEndRange(int dtaEndRange) {
        this.dtaEndRange = dtaEndRange;
    }
   
}
