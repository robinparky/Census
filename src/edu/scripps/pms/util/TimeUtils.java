package edu.scripps.pms.util;

import java.text.DateFormat;
import java.util.Calendar;
import java.util.Date;

/**
 * @author  Robin Park
 * @version $Id: TimeUtils.java,v 1.1 2007/03/14 00:38:04 rpark Exp $
 */

public class TimeUtils
{

    // this class can be used to monitor the time used by a process
    private static final DateFormat mediumFormat;
    private static final DateFormat dateMediumFormat;
    private static final DateFormat dateTimeFormat;
    private long startTime; // for timer
    private long endTime;    // for timer
    static
    {
        mediumFormat = DateFormat.getDateInstance(DateFormat.MEDIUM);
        dateMediumFormat = DateFormat.getDateInstance(DateFormat.MEDIUM);
        dateTimeFormat = DateFormat.getDateTimeInstance
	    (DateFormat.MEDIUM, DateFormat.SHORT);
    }

    public void startTiming() {
         startTime = System.currentTimeMillis();
    }
    public void stopTiming() {
         endTime = System.currentTimeMillis();
    }
    // get time used so far without reset the startTime
    public long getTimeUsed() {
        return System.currentTimeMillis() - startTime;
    }
    public long getTimeUsedMillis() {
        return endTime - startTime;
    }
    public long getTimeUsedSeconds() {
        return Math.round((endTime - startTime)/1000.0);
    }
    public long getTimeUsedMinutes() {
        return Math.round((endTime - startTime)/(60.0*1000));
    }
    public long getTimeUsedHours() {
        return Math.round((endTime - startTime)/(60*60.0*1000));
    }
    public static String getYearMonthDayPath ()
    {
        Calendar now = Calendar.getInstance();
        StringBuffer buf = new StringBuffer();

        buf.append(now.get(Calendar.YEAR));
        buf.append('/').append(now.get(Calendar.MONTH) + 1);
        buf.append('/').append(now.get(Calendar.DAY_OF_MONTH));

        return buf.toString();
    }

    public static Date getAheadDate (int howMany)
    {
        Calendar c = Calendar.getInstance();
        c.add(Calendar.DATE, howMany);
        return c.getTime();
    }

    public static Date getOffsetDate (Date date, int offset)
    {
        Calendar c = Calendar.getInstance();
        c.setTime(date);
        c.add(Calendar.DATE, offset);

        return c.getTime();
    }

    public static DateFormat getMediumFormat ()
    {
        return mediumFormat;
    }

    public static DateFormat getDateTimeFormat ()
    {
        return dateTimeFormat;
    }
}
