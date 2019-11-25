/*
 * CorrelationPlot.java
 *
 * Created on July 28, 2005, 11:00 PM
 */

package edu.scripps.pms.census.plot;

import ptolemy.plot.Plot;
import java.awt.*;
import java.util.*;
import javax.swing.JPanel;
import java.awt.geom.Line2D;
import java.awt.geom.AffineTransform;
import java.awt.event.*;
import edu.scripps.pms.census.util.LinearRegression;

/**
 *
 * @author rpark
 */
public class CorrelationPlot extends Plot {
    
    /** Creates a new instance of CorrelationPlot */
    private long[] arrX;
    private long[] arrY;
    private int startIndex;
    private int endIndex;
    private double slope;
    private double intercept;
    private int shift;
    
    public CorrelationPlot() {        
        super();
    }

    /*
    public CorrelationPlot(long[] arrX, long[] arrY) {
        
        this(arrX, arrY, 0, arrX.length-1);
    }
    
    public CorrelationPlot(long[] arrX, long[] arrY, int startIndex, int endIndex) {
        super();
        
        this.arrX = arrX;
        this.arrY = arrY;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
        
        corrInit();
    }
    */
    public void setData(long[] arrX, long[] arrY, int startIndex, int endIndex, double slope, double intercept, int shift)
    {
        this.arrX = arrX;
        this.arrY = arrY;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
        this.slope = slope;
        this.intercept = intercept;
        this.shift = shift;
        
        corrInit();
    }
    
    public synchronized void addMouseListener(MouseListener l) {  
        //disable zooming
    }
    
    private void corrInit()
    {
        this.setGrid(false);
        
        //_drawPlotPoint(this.getGraphics());
    }
    
    private void _drawPlotPoint(Graphics graphics)
            //int dataset, int index) {
    {
        graphics.drawLine(0, 0, 15, 20);
        
        /*
        // Set the color
        if (_usecolor) {
            int color = dataset % _colors.length;
            graphics.setColor(_colors[color]);
        } else {
            graphics.setColor(_foreground);
        }
*/
        
        //Vector pts = (Vector)_points.elementAt(dataset);
        //PlotPoint pt = (PlotPoint)pts.elementAt(index);
        // Use long here because these numbers can be quite large
        // (when we are zoomed out a lot).
        long ypos = _lry - (long)((10 - _yMin) * _yscale);
        long xpos = _ulx + (long)((5 - _xMin) * _xscale);
        
        // Draw the line to the previous point.
        //long prevx = ((Long)_prevx.elementAt(dataset)).longValue();
        //long prevy = ((Long)_prevy.elementAt(dataset)).longValue();
        // MIN_VALUE is a flag that there has been no previous x or y.
        //if (pt.connected) {
            //_drawLine(graphics, dataset, xpos, ypos, prevx, prevy, true);

            //            rectVec.add(new Rectangle((int)xpos-5, (int)ypos-5, 5, 5));
            //            graphics.drawRect((int)xpos-5, (int)ypos-5, 5, 5);

          //  int[] x = {(int)prevx-1, (int)prevx+1, (int)xpos+1, (int)xpos-1};
          //  int[] y = {(int)prevy-1, (int)prevy+1, (int)ypos+1, (int)ypos-1};

            //            xValues.add( new Double(xpos) );
            //Robin

            //graphics.drawPolygon( x, y, 4);

            //            rectVec.add(new Polygon(x, y, 4));
        //}
        ///comehere

        //        graphics.drawString("" + xpos + " " + ypos, (int)xpos, (int)ypos);

        // Save the current point as the "previous" point for future
        // line drawing.
        //_prevx.setElementAt(new Long(xpos), dataset);
        //_prevy.setElementAt(new Long(ypos), dataset);
/*
        // Draw decorations that may be specified on a per-dataset basis
        Format fmt = (Format)_formats.elementAt(dataset);
        if (fmt.impulsesUseDefault) {
            if (_impulses) _drawImpulse(graphics, xpos, ypos, true);
        } else {
            if (fmt.impulses) _drawImpulse(graphics, xpos, ypos, true);
        }

        // Check to see whether the dataset has a marks directive
        int marks = _marks;
        if (!fmt.marksUseDefault) marks = fmt.marks;
        if (marks != 0) _drawPoint(graphics, dataset, xpos, ypos, true);

        if (_bars) _drawBar(graphics, dataset, xpos, ypos, true);
        if (pt.errorBar)
            _drawErrorBar(graphics, dataset, xpos,
                    _lry - (long)((pt.yLowEB - _yMin) * _yscale),
                    _lry - (long)((pt.yHighEB - _yMin) * _yscale), true);

        // Restore the color, in case the box gets redrawn.
        graphics.setColor(_foreground);
        if (_pointsPersistence > 0 || _xPersistence > 0.0) {
            // Restore paint mode in case axes get redrawn.
            graphics.setPaintMode();
        }
        */
    }
    
    protected synchronized void _drawPlot(
        Graphics graphics, boolean clearfirst, Rectangle drawRect) {

        // Ignore if there is no graphics object to draw on.
        if (graphics == null) return;

        graphics.setPaintMode();

        // Make sure we have an x and y range
        if (!_xRangeGiven) {
            if (_xBottom > _xTop) {
                // have nothing to go on.
                _setXRange(0, 1000);
            } else {
                _setXRange(_xBottom, _xTop);
            }
        }
        if (!_yRangeGiven) {
            if (_yBottom > _yTop) {
                // have nothing to go on.
                _setYRange(0, 1000);
            } else {
                _setYRange(_yBottom, _yTop);
            }
        }
        
        //_setXRange(0, 100);
        //_setYRange(0, 100);
        
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

                xSPos -= _labelFontMetrics.stringWidth("x10");
                xlocation -= _labelFontMetrics.stringWidth("x10");
                graphics.setFont(_labelFont);
                //graphics.drawString("x10", xSPos, ySPos);
                //graphics.drawString("x10", xlocation, ySPos);
            }

            //////////Robin come back here to fix padding - original + 5, not 25

            // NOTE: 5 pixel padding on bottom
            _bottomPadding = (3 * labelheight)/2;
        }

        // NOTE: 5 pixel padding on the bottom.
        if (_xlabel != null && _bottomPadding < labelheight + 5) {
            _bottomPadding = labelheight + 5;
        }

        // Compute the space needed around the plot, starting with vertical.
        // NOTE: padding of 5 pixels below title.
        _uly = 20; //titley + 5;
        // NOTE: 3 pixels above bottom labels.
        _lry = 230;//drawRect.height-labelheight-_bottomPadding;
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

            int numfracdigits = _numFracDigits(yStep);

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

        // Next we do the horizontal spacing.
        /*
        if (_ylabel != null) {
            _ulx = widesty + _labelFontMetrics.stringWidth("W") + _leftPadding;
        } else {
            _ulx = widesty + _leftPadding;
        }
        */
        _ulx=20;
        //int legendwidth = _drawLegend(graphics,
        //        drawRect.width-_rightPadding, _uly);
        _lrx = drawRect.width -10; //_rightPadding;
        _lrx = 230; //drawRect.width -10;
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
       
        //vertical axis
        /*
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
                //graphics.drawString(ylabels[ind],
                //        _ulx-ylabwidth[ind++]-4, yCoord1+offset);
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
        
        } else {
            // ticks have been explicitly specified
            Enumeration nt = _yticks.elements();
            Enumeration nl = _yticklabels.elements();

            
        Vector ntList = new Vector();
        ntList.add(new Double("0.01"));
        ntList.add(new Double("0.02"));
        Vector nlList = new Vector();
        nlList.add("40");
        nlList.add("80");

            Enumeration nt = ntList.elements(); //new Enumeration
            Enumeration nl = nlList.elements();
        */            
            /*
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
                //graphics.drawString(label,
                //        _ulx - _labelFontMetrics.stringWidth(label) - 3,
                //        yCoord1+offset);
            }
        }

*/
        //////////////////// horizontal axis
        int yCoord1 = _uly+tickLength;
        int yCoord2 = _lry-tickLength;
        int charwidth = _labelFontMetrics.stringWidth("8");
        
        /*
        if (_xticks == null) {
            // auto-ticks

            // Number of x tick marks.
            // Need to start with a guess and converge on a solution here.
            int nx = 10;
            double xStep = 0.0;
            int numfracdigits = 0;
            if (_xlog) {
                // X axes log labels will be at most 6 chars: -1E-02
                nx = 2 + width/((charwidth * 6) + 10);
            } else {
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
            }
            xStep = _roundUp((_xtickMax-_xtickMin)/(double)nx);
            numfracdigits = _numFracDigits(xStep);

            // Compute x starting point so it is a multiple of xStep.
            double xStart = xStep*Math.ceil(_xtickMin/xStep);

            // NOTE: Following disables first tick.  Not a good idea?
            // if (xStart == _xMin) xStart += xStep;

            Vector xgrid = null;
            double xTmpStart = xStart;
            if (_xlog) {
                xgrid = _gridInit(xStart, xStep, true, null);
                //xgrid = _gridInit(xStart, xStep);
                xTmpStart = _gridRoundUp(xgrid, xStart);
            }

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

            if (_xlog) {
                // Draw in grid lines that don't have labels.

                // If the step is greater than 1, clamp it to 1 so that
                // we draw the unlabeled grid lines for each
                // integer interval.
                double tmpStep = (xStep > 1.0)? 1.0 : xStep;

                // Recalculate the start using the new step.
                xTmpStart = tmpStep*Math.ceil(_xtickMin/tmpStep);

                Vector unlabeledgrid  = _gridInit(xTmpStart, tmpStep,
                        false, xgrid);
                if (unlabeledgrid.size() > 0 ) {
                    for (double xpos = _gridStep(unlabeledgrid, xTmpStart,
                            tmpStep, _xlog);
                         xpos <= _xtickMax;
                         xpos = _gridStep(unlabeledgrid, xpos,
                                 tmpStep, _xlog)) {
                        xCoord1 = _ulx + (int)((xpos-_xtickMin)*_xtickscale);
                        if (_grid && xCoord1 != _ulx && xCoord1 != _lrx) {
                            graphics.setColor(Color.lightGray);
                            graphics.drawLine(xCoord1, _uly+1,
                                    xCoord1, _lry-1);
                            graphics.setColor(_foreground);
                        }
                    }
                }

            }


        } else {
            // ticks have been explicitly specified
            Enumeration nt = _xticks.elements();
            Enumeration nl = _xticklabels.elements();
            // Code contributed by Jun Wu (jwu@inin.com.au)
            double preLength = 0.0;

            while (nl.hasMoreElements()) {
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

                // If the labels are not overlapped, proceed.
                //if (labxpos > preLength) {
                {
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

                    Font font = new Font("",Font.PLAIN,8);

                    Graphics2D g2d = (Graphics2D)graphics;
                    g2d.setFont(font);
                    // NOTE: Fudge factor so label doesn't touch axis labels.
                    g2d.rotate(Math.toRadians(-90), labxpos, labypos);
                    //g2d.drawString(label, labxpos, labypos);
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
        }
*/
        graphics.setFont(previousFont);

    //    graphics.drawLine(_ulx, _uly, 100, 50);
    //    graphics.drawLine(_lrx, _lry, 100, 50);
    //    graphics.drawLine(_ulx, _lry, 100, 50);

        
        int dia=4;
        double x, y;
        long max=0;
        long min=arrX[startIndex];
        
        long[] shiftedYArr = new long[endIndex-startIndex+1];
        int j=0;
        for(int i=startIndex;i<=endIndex;i++)
        {
            int tempShift = i+shift;
            if(tempShift<0) tempShift=0;
            if(tempShift>=arrX.length) tempShift=arrX.length-1;
            
            shiftedYArr[j++] = arrY[tempShift];
            
            if(max<arrX[i]) max = arrX[i];
                
            if(max<arrY[tempShift]) max = arrY[tempShift];
            
            if(min>arrX[i]) min = arrX[i];
            
            if(min>arrY[tempShift]) min = arrY[tempShift];
        }
        
        
        double scale = (double)(_lrx-_ulx-10)/(max-min); //210 is width
        double adjust =(min-_ulx)*scale;
        //System.out.println("=====+++++++++++" + min + "\t" + max + "\t" + _ulx + "\t" + scale);
        
        j=0;
        double xSum=0;
        double ySum=0;
        double xxSum=0;
        double yySum=0;
        double xySum=0;
        
        for(int i=startIndex;i<=endIndex;i++)
        {            
            x = _ulx + (scale*arrX[i] - adjust); //(int)((width * arrX[i])/max);
            y = _lry - (scale*shiftedYArr[j] - adjust); //(int)((width * arrY[i])/max);
            
            xSum += (scale*arrX[i] - adjust);
            ySum += (scale*shiftedYArr[j] - adjust);
            xxSum += (scale*arrX[i] - adjust)*(scale*arrX[i] - adjust);
            yySum += (scale*shiftedYArr[j] - adjust)*(scale*shiftedYArr[j] - adjust);
            xySum += (scale*arrX[i] - adjust)*(scale*shiftedYArr[j] - adjust);
                
            //System.out.println("[[[[" + "\t" + arrX[i] + "\t" + shiftedYArr[j-1] + "\t" + x +"\t" + y);
             // System.out.println( "===>>" + "\t" + (scale*arrX[i] - adjust) + "\t" + (scale*shiftedYArr[j] - adjust));
            graphics.setColor(Color.BLUE);
            graphics.fillOval( (int)(x-dia/2), (int)(y-dia/2), dia, dia);    
            
            j++;
        }
        
        double tempYAve = ySum/(endIndex-startIndex+1);
        double tempXAve = xSum/(endIndex-startIndex+1);
        double tempSlope = (xySum-(double)xSum*tempYAve)/(xxSum-xSum*tempXAve);
        double tempIntercept = (tempYAve+ -tempSlope*tempXAve);
        
//                intercept = yAve - slope*xAve;
                
        
        
        int startX;
        int startY;
        int endX;
        int endY;
        
        /*
        System.out.println("data start");
        double temp1=0;
        double temp2=0;
        
        /*
        if(startIndex<0) startIndex=0;
        if(startIndex>=arrX.length) startIndex=arrX.length-1;
        if(endIndex<0) endIndex=0;
        if(endIndex>=arrX.length) endIndex=arrX.length-1;
        
        j=0;
        for(int i=this.startIndex;i<=this.endIndex;i++)
        {
            System.out.println(i  + "\t" + this.arrX[i] + "\t" + shiftedYArr[j++]);
            
        }
        System.out.println("intercept===>>" + intercept + "\t" + intercept*scale + "\t" + slope + "\t" + adjust);
        
        */
        
        
        if(intercept>0)
        {
            startX = _ulx;
            startY = _lry - (int)Math.round(tempIntercept);
            
            //x = _ulx +    (scale*arrX[i]) - adjust; //((width * arrX[i])/max);
            //y = _lry - (scale*arrY[i]) + adjust; //((width * arrY[i])/max);
        } 
        else
        {            
            //startX = _ulx + (intercept*scale/slope - adjust);
                                    
            startX = (int)Math.round((-tempIntercept)/tempSlope + _ulx);
            startY = _lry;
        }
        
        
        //Find end point for the slope
        if((width*tempSlope+tempIntercept)>height)
        {
            //endX = _ulx + (int)((_lry-intercept*scale*scale-adjust)/slope);
            //endY = _lry - (int)(_lrx*slope+intercept*scale*scale + adjust);
            endX = _ulx + (int)Math.round( (height-tempIntercept)/tempSlope );            
            endY = _uly; //_lry - (int)(_lrx*slope+intercept*scale*scale + adjust);
        }
        else
        {
            endX = _lrx;
            endY = _lry - (int)Math.round(_lrx*tempSlope+tempIntercept); 
        }
        
//        System.out.println("=======" + startX + "\t" + startY +"\t" + endX + "\t" + endY);
        //System.out.println("intercept*scale===>>" + intercept*scale + "\t" + startX + "\t" + startY);
        
        graphics.setColor(Color.RED);
        graphics.drawLine(startX, startY, endX, endY);            
        
        //plot Legend
        graphics.setFont(new Font("Helvetica",Font.BOLD,11));
        graphics.setColor(Color.BLUE);
        String sample = "Sample Intensity";
        String ref = "Reference";
        String scanNum = "Scan #";
        
        
        
        Graphics2D g = (Graphics2D)graphics;
        g.setFont(new Font("Helvetica",Font.BOLD,11));
        g.setColor(Color.RED);
        g.drawString("Reference Intensity", _ulx + 50, _lry+10);
        AffineTransform at = new AffineTransform();
        at.setToRotation(-Math.PI/2.0);
        g.translate(getWidth()/2, getHeight()/2);
        g.transform(at);//.rotate(-Math.PI/2);
        g.setColor(Color.BLUE);
        g.drawString("Sample Intensity", _ulx-70, _uly-130);

    }
}
