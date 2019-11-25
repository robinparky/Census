package edu.scripps.pms.census.util.io;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: Yates Lab</p>
 *
 * @author Robin Park
 * @version $Id: FileUtil.java,v 1.1 2014/09/09 19:29:52 rpark Exp $
 *
 */
public class FileUtil
{
    private FileUtil()
    {
    }

    public static void copy(String in, String out, boolean isAppending) throws IOException {
        copy(new File(in), new File(out), isAppending);
    }

    public static void copy(File in, File out, boolean isAppending) throws IOException {
        FileInputStream fis  = new FileInputStream(in);
        FileOutputStream fos = new FileOutputStream(out, isAppending);
        byte[] buf = new byte[4096];
        int i = 0;

        while((i=fis.read(buf))!=-1) {
          fos.write(buf, 0, i);
        }

        fis.close();
        fos.close();
    }

    
    public static void makeDir(String targetDir) {
		File dataFile = new File(targetDir);
		if (!dataFile.exists()) {
			dataFile.mkdirs();
		}

	}
    public static void writeJSON(String jsonData, String jsonFileName) {

        try {
            BufferedWriter writer = null;
            File f = new File(jsonFileName);
            if (f.exists()) {
                f.delete();
            }
            writer = new BufferedWriter(new FileWriter(jsonFileName, true));
//			writer.write(JsonWriter.formatJson(jsonData));
            writer.write(jsonData);
            writer.write("\n");
            writer.close();
        } catch (IOException ex) {
            System.out.println("writeJSON :: "+ex.getMessage());
        }catch (Exception e) {
            e.printStackTrace();
        }
    }

    //This method is using new lib from jdk and supposed to be faster then old lib.
    //But before using it, make sure it is working correctly.
    //Currently we are not using this method.
    /*
    public static void copy(FileInputStream source, FileOutputStream dest) throws IOException {
         FileChannel in = null, out = null;
         try
         {
              in = source.getChannel();
              out = dest.getChannel();

              long size = in.size();
              MappedByteBuffer buf = in.map(FileChannel.MapMode.READ_ONLY, 0, size);

              out.write(buf);

         }
         catch(IOException e)
         {
             throw new IOException(e.toString());
         }
         finally {
              if (in != null)          in.close();
              if (out != null)     out.close();
         }
    }
    */
}
