/*
 * MSIndexFileCreator.java
 *
 * Created on March 25, 2005, 10:36 AM
 */

package edu.scripps.pms.census.hash;

/**
 * @author Robin Park
 * @version $Id: MSIndexFileCreator.java,v 1.10 2014/02/21 22:10:52 rpark Exp $
 */
import java.io.*;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;

public class MSIndexFileCreator {

    private final static String MS1_FILE = "ms1";
    private final static String MS2_FILE = "ms2";
    private final static String MSZM = "mszm";
    public static final String NEW_LINE = System.getProperty("line.separator");
    public static final char TAB = '\t';

    private static final int
            CR  = 13,
            LF  = 10;
    
    private String filePath;
    private static java.text.NumberFormat n = java.text.NumberFormat.getInstance();

    static {
        n.setMaximumFractionDigits(3);
        n.setGroupingUsed(false);
    }

    private static void printError() {
        System.out.println("Usage : java MSIndexFileCreator [-option]");
        System.out.println("-a\t Generate index files for all ms files in the same folder");
        System.out.println("-f\t ms file name");

    }

    public static void main(String args[]) throws IOException {


        createIndexFile(args[0]);
       // createIndexFile("/home/rpark/test_data/berkeley_lori/72/spectra/split/1_399-JB01_02_102816.ms1");

        if(true) return;

        if (args.length < 1 || ("-f".equals(args[0]) && args.length < 2)) {
            printError();

            return;
        }



        System.out.println("indexing files...");

        if ("-a".equals(args[0])) {
            File f = new File(".");
            String[] list = f.list();

            for (int i = 0; i < list.length; i++) {
                if (null != list[i] && (list[i].endsWith(MS1_FILE) || list[i].endsWith(MS2_FILE) || list[i].endsWith(MSZM))) {
                    createIndexFile(list[i]);
                    System.out.println(list[i] + " was indexed.");
                }

            }

            return;
        } else if ("-f".equals(args[0])) {
            createIndexFile(args[1]);
            System.out.println(args[1] + " was indexed.");
        } else {
            printError();
            return;
        }
        
        System.out.println("completed");

    }

    public static void createIndexFileByPath(String path) throws Exception {

//        System.out.println("COME BACK AND REMOVE THIS DEBUGGING LINE");
  //      if(true) return;


        File f = new File(path);
        String[] list = f.list();

        System.out.print("indexing spectral files");
        for (int i = 0; i < list.length; i++) {
            if (null != list[i] && (list[i].endsWith(MS1_FILE) || list[i].endsWith(MS2_FILE) || list[i].endsWith(MSZM))) {

                /*
System.out.println("==================" + path+File.separator+list[i]);
try {
throw new Exception ();
} catch(Exception e) {
e.printStackTrace();
System.exit(0);
} */

                createIndexFile(path+File.separator+list[i]);
                System.out.print(".");
            }

        }
        System.out.println("\ncomplete");
    }

    public MSIndexFileCreator(String filePath) {
        this.filePath = filePath;
    }

    public static void createIndexFile(String fileName) throws IOException {

        try {
            if (fileName.endsWith("ms1") || fileName.endsWith("mszm")) {
                createMS1IndexFileRevised(fileName);
            } else if (fileName.endsWith("ms2") || fileName.endsWith("ms3")) {
                createMS2IndexFileRevised(fileName);
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IOException("Failed to create index file. " + e);
        }
    }

    private static void createMS1IndexFile(String fileName) throws IOException {
        File file = new File(fileName);
        InputStream fisIn = new FileInputStream(fileName);

        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String lastLine;

        FileOutputStream out = new FileOutputStream(fileName + ".index");
        PrintStream p = new PrintStream(out);
        String temp;

        int size = (int) file.length();

        byte[] byteBuffer = new byte[size];

        fisIn.read(byteBuffer, 0, size);
        StringBuffer sb = new StringBuffer();

        long pos = 0;

        for (int i = 0; i < byteBuffer.length; i++) {
            if ((char) byteBuffer[i] == 'S' && ((i == 0) || (char) byteBuffer[i - 1] == '\n')) {
                pos = i;
                int j = i; //skip S char and space
                char ch = (char) byteBuffer[j];
                int tabCount = 0;

                while (ch != '\n') {
                    if (ch == '\t') {
                        tabCount++;
                    }

                    if (tabCount == 2) {
                        sb.append(ch);
                    }

                    ch = (char) byteBuffer[++j];

                }
                
                p.print(Integer.parseInt(sb.toString().trim()));
                p.print("\t");                
                p.println(pos);

                sb.delete(0, sb.length());
            }
           
        }

        fisIn.close();
        p.close();
        out.close();
    }

    private static void createMS2IndexFile(String fileName) throws IOException {

        File file = new File(fileName);
        InputStream fisIn = new FileInputStream(fileName);

        BufferedReader br = new BufferedReader(new FileReader(fileName));

        String lastLine = br.readLine();
        String[] strArr = lastLine.split("\t");

        String extension = fileName.substring(fileName.lastIndexOf(".") + 1);

        boolean isOldVersion = false;

        if (strArr.length > 3) {
            lastLine = strArr[2].substring(strArr[2].lastIndexOf("/") + 1);
            lastLine = lastLine.substring(0, lastLine.indexOf(" "));
            int year = Integer.parseInt(lastLine);

            isOldVersion = (year < 2006) ? true : false;
        }

        FileOutputStream out = new FileOutputStream(fileName + ".index");
        PrintStream p = new PrintStream(out);

        //String temp;
        int size = (int) file.length();

        byte[] byteBuffer = new byte[size];

        fisIn.read(byteBuffer, 0, size);
        StringBuffer sb = new StringBuffer();
        StringBuffer precurSb = new StringBuffer();

        long pos;
        double prevRet = -1;
        double inc = 0.000d;

        for (int i = 0; i < byteBuffer.length; i++) {
            if ((char) byteBuffer[i] == 'S' && (char) byteBuffer[i - 1] == '\n') {
                pos = i;
                int j = i; //skip S char and space
                char ch = (char) byteBuffer[j];
                int tabCount = 0;

                while (ch != '\n') {
                    if (ch == '\t') {
                        tabCount++;
                    }

                    if (tabCount == 2) {
                        sb.append(ch);
                    } else if (tabCount == 3) {
                        precurSb.append(ch);
                    }

                    ch = (char) byteBuffer[++j];

                }

                p.print(Integer.parseInt(sb.toString().trim()));
                p.print("\t");
                p.print(pos);
                p.print("\t");
                p.print(Float.parseFloat(precurSb.toString().trim()));

                if (isOldVersion) {
                    p.println("");
                } else {
                    p.print("\t");
                }

                sb.delete(0, sb.length());
                precurSb.delete(0, precurSb.length());

            }

            if (!isOldVersion && (char) byteBuffer[i] == 'I'
                    && (char) byteBuffer[i + 1] == '\t'
                    && (char) byteBuffer[i + 2] == 'R'
                    && (char) byteBuffer[i + 3] == 'e'
                    && (char) byteBuffer[i + 4] == 't'
                    && (char) byteBuffer[i + 5] == 'T') {

                pos = i;
                int j = i;
                char ch = (char) byteBuffer[j];
                int tabCount = 0;

                while (ch != '\n') {
                    if (ch == '\t') {
                        tabCount++;
                    }

                    if (tabCount == 2) {
                        sb.append(ch);
                    }

                    ch = (char) byteBuffer[++j];

                }

                double curRet = Double.parseDouble(sb.toString().trim());
                if (curRet == prevRet) {
                    inc += 0.001d;
                    curRet += inc;
                } else {
                    inc = 0.000d;
                    prevRet = curRet;
                }

                p.print(n.format(curRet));
                p.print("\t");

                sb.delete(0, sb.length());
            } else if (!isOldVersion && (char) byteBuffer[i] == 'I'
                    && (char) byteBuffer[i + 1] == '\t'
                    && (char) byteBuffer[i + 2] == 'A'
                    && (char) byteBuffer[i + 3] == 'c'
                    && (char) byteBuffer[i + 4] == 't'
                    && (char) byteBuffer[i + 5] == 'i'
                    && (char) byteBuffer[i + 6] == 'v') {

                pos = i;
                int j = i;
                char ch = (char) byteBuffer[j];
                int tabCount = 0;

                while (ch != '\n') {
                    if (ch == '\t') {
                        tabCount++;
                    }

                    if (tabCount == 2) {
                        sb.append(ch);
                    }

                    ch = (char) byteBuffer[++j];

                }

                p.print(sb.toString().trim());
                if ("ms3".equals(extension) || "ms2".equals(extension)) {
                    p.print("\t");
                } else {
                    p.print("\n");
                }

                sb.delete(0, sb.length());
            } else if (!isOldVersion && (char) byteBuffer[i] == 'I'
                    && (char) byteBuffer[i + 1] == '\t'
                    && (char) byteBuffer[i + 2] == 'P'
                    && (char) byteBuffer[i + 3] == 'r'
                    && (char) byteBuffer[i + 4] == 'e'
                    && //(char)byteBuffer[i+5] == 'c' &&
                    //(char)byteBuffer[i+6] == 'u' &&
                    //(char)byteBuffer[i+7] == 'r' &&
                    //(char)byteBuffer[i+8] == 's' &&
                    //(char)byteBuffer[i+9] == 'o' &&
                    //(char)byteBuffer[i+10] == 'r' &&
                    (char) byteBuffer[i + 11] == 'S'
                    && (char) byteBuffer[i + 12] == 'c'
                    && (char) byteBuffer[i + 13] == 'a'
                    && (char) byteBuffer[i + 14] == 'n') {

                pos = i;
                int j = i;
                char ch = (char) byteBuffer[j];
                int tabCount = 0;

                while (ch != '\n') {
                    if (ch == '\t') {
                        tabCount++;
                    }

                    if (tabCount == 2) {
                        sb.append(ch);
                    }

                    ch = (char) byteBuffer[++j];

                }

                p.println(sb.toString().trim());

                sb.delete(0, sb.length());
            }
        }

        p.close();
        out.close();
        fisIn.close();
    }

    private static void createMS1IndexFileRevised(String fileName) throws IOException, ParseException {

       // System.out.println("===================");
        InputStream fisIn = new FileInputStream(fileName);
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        FileOutputStream out = new FileOutputStream(fileName + ".index");
        PrintStream p = new PrintStream(out);
        RandomAccessFile raf = new RandomAccessFile(fileName , "r");

        StringBuilder sb = new StringBuilder();
        String currentLine;
        long pos = 0;
        double prevRet = -1;
        double inc = 0.000d;
        //int lineSepratorLength = NEW_LINE.length();
        boolean hasPrinted = true;
        
        while ((currentLine = br.readLine()) != null) {
            int currentLineLength = currentLine.length();
            if (currentLine.charAt(0) == 'S') {

                if(!hasPrinted)
                {
                    p.print("\n");
                    sb.setLength(0);
                }
                hasPrinted =false;

                int tabCount = 0;                
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                
                p.print(Integer.parseInt(sb.toString().trim()));
                p.print(TAB);
                p.print(pos);
                //p.print(NEW_LINE);
                sb.setLength(0);                
            }
            
            if ( currentLine.charAt(0) == 'I'
                    && currentLine.charAt(1) == '\t'
                    && currentLine.charAt(2) == 'R'
                    && currentLine.charAt(3) == 'e'
                    && currentLine.charAt(4) == 't'
                    && currentLine.charAt(5) == 'T') {

                int tabCount = 0;
                //int currentLineLength = currentLine.length();
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                double curRet = Double.parseDouble(sb.toString().trim());
                if (curRet == prevRet) {
                    inc += 0.001d;
                    curRet += inc;
                } else {
                    inc = 0.000d;
                    prevRet = curRet;
                }
                p.print(TAB);
//                p.print(pos);
                p.print(n.format(curRet));
                
                
                sb.setLength(0);
            }
            
            if ( currentLine.charAt(0) == 'I'
                    && currentLine.charAt(1) == '\t'
                    && currentLine.charAt(2) == 'I'
                    && currentLine.charAt(3) == 'o'
                    && currentLine.charAt(4) == 'n'
                    && currentLine.charAt(5) == 'I') {

                int tabCount = 0;
                //int currentLineLength = currentLine.length();
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                double curRet = Double.parseDouble(sb.toString().trim());
                if (curRet == prevRet) {
                    inc += 0.001d;
                    curRet += inc;
                } else {
                    inc = 0.000d;
                    prevRet = curRet;
                }
                p.print(TAB);
//                p.print(pos);
                p.print(n.format(curRet));
                p.print("\n");


                sb.setLength(0);
                hasPrinted = true;
            }

            raf.seek(pos + currentLineLength);
            if((int)raf.read() == CR && (int)raf.read() == LF){
                pos += currentLineLength + 2;
            }else{
                pos += currentLineLength + 1;
            }
        }

        if(null != raf)
            raf.close();
        if(null != p)
            p.close();
        if(null != out)
            out.close();
        if(null != br)
            br.close();
        if(null != fisIn)
            fisIn.close();
    }

    private static void createMS2IndexFileRevised(String fileName) throws IOException, ParseException {

        InputStream fisIn = new FileInputStream(fileName);
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        FileOutputStream out = new FileOutputStream(fileName + ".index");
        RandomAccessFile raf = new RandomAccessFile(fileName , "r");
        
        PrintStream p = new PrintStream(out);
        String extension = fileName.substring(fileName.lastIndexOf(".") + 1);

        StringBuilder sb = new StringBuilder();
        StringBuilder precurSb = new StringBuilder();
        String currentLine;
        long pos = 0;
        double prevRet = -1;
        double inc = 0.000d;
        //int lineSepratorLength = NEW_LINE.length();
        boolean isOldVersion = false;
        boolean anyLineAdded = false;
        int scan = -1;
        Map<Integer,Integer> scanMap = new HashMap<>();
        while ((currentLine = br.readLine()) != null) {
            int currentLineLength = currentLine.length();
            if (pos == 0) {
                String[] strArr = currentLine.split(String.valueOf(TAB));

                if (strArr.length > 3) {
                    Calendar cal = Calendar.getInstance();
                    SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy");
                    cal.setTime(sdf.parse(strArr[2]));
                    int year = cal.get(Calendar.YEAR);
                    isOldVersion = year < 2006;
                }
            }

            if (currentLine.charAt(0) == 'S') {
                if(anyLineAdded){
                    p.print(NEW_LINE);
                }
                int tabCount = 0;
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    } else if (tabCount == 3) {
                        precurSb.append(ch);
                    }
                }
                scan = Integer.parseInt(sb.toString().trim());
                p.print(scan);
                p.print(TAB);
                p.print(pos);
                p.print(TAB);
                p.print(Float.parseFloat(precurSb.toString().trim()));

                if (isOldVersion) {
                    p.println("");
                } else {
                    p.print(TAB);
                }
                anyLineAdded = true;
                sb.setLength(0);
                precurSb.setLength(0);
            }

            if (!isOldVersion && currentLine.charAt(0) == 'I'
                    && currentLine.charAt(1) == '\t'
                    && currentLine.charAt(2) == 'R'
                    && currentLine.charAt(3) == 'e'
                    && currentLine.charAt(4) == 't'
                    && currentLine.charAt(5) == 'T') {

                int tabCount = 0;
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                double curRet = Double.parseDouble(sb.toString().trim());
                if (curRet == prevRet) {
                    inc += 0.001d;
                    curRet += inc;
                } else {
                    inc = 0.000d;
                    prevRet = curRet;
                }
              //  System.out.println("testing"+n.format(curRet));
                p.print(n.format(curRet));
                p.print(TAB);

                sb.setLength(0);
            } else if (!isOldVersion && currentLine.charAt(0) == 'I'
                    && currentLine.charAt(1) == '\t'
                    && currentLine.charAt(2) == 'A'
                    && currentLine.charAt(3) == 'c'
                    && currentLine.charAt(4) == 't'
                    && currentLine.charAt(5) == 'i'
                    && currentLine.charAt(6) == 'v') {

                int tabCount = 0;
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                p.print(sb.toString().trim());
                p.print(TAB);
                
                sb.setLength(0);
            } else if (!isOldVersion && currentLine.charAt(0) == 'I'
                    && currentLine.charAt(1) == '\t'
                    && currentLine.charAt(2) == 'P'
                    && currentLine.charAt(3) == 'r'
                    && currentLine.charAt(4) == 'e'
                    && //(char)byteBuffer[i+5] == 'c' &&
                    //(char)byteBuffer[i+6] == 'u' &&
                    //(char)byteBuffer[i+7] == 'r' &&
                    //(char)byteBuffer[i+8] == 's' &&
                    //(char)byteBuffer[i+9] == 'o' &&
                    //(char)byteBuffer[i+10] == 'r' &&
                    currentLine.charAt(11) == 'S'
                    && currentLine.charAt(12) == 'c'
                    && currentLine.charAt(13) == 'a'
                    && currentLine.charAt(14) == 'n') {

                int tabCount = 0;
                for (int i = 0; i < currentLineLength; i++) {
                    char ch = currentLine.charAt(i);
                    if (ch == TAB) {
                        tabCount++;
                    }
                    if (tabCount == 2) {
                        sb.append(ch);
                    }
                }
                int precursorScan = Integer.parseInt(sb.toString().trim());
                if(!scanMap.containsKey(precursorScan))
                {
                    scanMap.put(scan,precursorScan);
                }
                else
                {
                    precursorScan = scanMap.get(precursorScan);
                    scanMap.put(scan,precursorScan);
                }
                p.print(precursorScan);
                sb.setLength(0);
            }
            
            raf.seek(pos + currentLineLength);
            if((int)raf.read() == CR && (int)raf.read() == LF){
                pos += currentLineLength + 2;
            }else{
                pos += currentLineLength + 1;
            }
        }

        fisIn.close();
        p.close();
        out.close();
    }

}
