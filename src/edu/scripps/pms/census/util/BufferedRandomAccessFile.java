/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.scripps.pms.census.util;

/**
 *
 * @author rpark
 */


import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;

public class BufferedRandomAccessFile extends RandomAccessFile
{
    boolean reading=true;
    private byte buffer[];
    private int bufferSize = 0;

    private long filePos=0;  
    private long fileLength=0;
    private long bufferStart=0;

    /**
     *
     * @param filename
     * @param mode
     * @param bufsize
     * @throws IOException
     */
    public BufferedRandomAccessFile(String filename, String mode, int bufsize) throws IOException 
    {
        this(new File(filename),mode,bufsize);
    }

    public BufferedRandomAccessFile(File file, String mode, int bufsize) throws IOException 
    {
        super(file, mode);
        fileLength=file.length();
        buffer = new byte[bufsize];
    }

    public final int read() throws IOException 
    {
        if (!reading) switchToReadBuffer();
        while(true)
        {
            if (filePos==fileLength) return -1;
            // read the data
            int readAtIdx=(int) (filePos-bufferStart);
            if (readAtIdx<0 || readAtIdx>=bufferSize)
                updateReadBuffer();
            else
            {
                ++filePos;
                return ((int)buffer[readAtIdx]) & 0xff;
            }
        }
    }

    @Override
    public int read(byte[] b, int off, int len) throws IOException 
    {
        int fileAvailable=(int) (fileLength-filePos);
        if (fileAvailable==0) return -1;
        if (!reading) switchToReadBuffer();
        if (len>fileAvailable) len=fileAvailable;
        int readAtIdx=(int) (filePos-bufferStart);
        if (readAtIdx<0 || readAtIdx>=bufferSize)
        {
            updateReadBuffer();
            readAtIdx=(int) (filePos-bufferStart);
        }
        int availableInBuffer=bufferSize-readAtIdx;
        if (len>availableInBuffer) len=availableInBuffer;
        System.arraycopy(buffer, readAtIdx, b, off, len);
        filePos+=len;
        return len;
    }

    @Override
    public void write(int b) throws IOException 
    {
        if (reading) 
            switchToWriteBuffer();
        while(true)
        {
            if (bufferSize==0)
                bufferStart=filePos;
            int writeAtIdx=(int) (filePos-bufferStart);
            if (writeAtIdx<0 || writeAtIdx>=buffer.length)
                flush();
            else
            {
                buffer[writeAtIdx]=(byte) b;
                if (writeAtIdx==bufferSize) bufferSize++;
                if (++filePos>fileLength) fileLength=filePos;
                return;
            }
        }
    };

    @Override
    public void write(byte[] b, int off, int len) throws IOException 
    {
        if (reading) 
            switchToWriteBuffer();
        int from=off;
        int remaining=len;
        while(remaining>0)
        {
            if (bufferSize==0)
                bufferStart=filePos;
            int writeAtIdx=(int) (filePos-bufferStart);
            if (writeAtIdx<0 || writeAtIdx>=buffer.length)
                flush();
            else
            {
                int todo=buffer.length-writeAtIdx;
                if (todo>remaining) todo=remaining;
                System.arraycopy(b, from, buffer, writeAtIdx, todo);
                writeAtIdx+=todo;
                if (writeAtIdx>bufferSize) bufferSize=writeAtIdx;
                filePos+=todo;
                if (filePos>fileLength) fileLength=filePos;
                remaining-=todo;
                from+=todo;
            }
        }
    }

    private void switchToWriteBuffer()
    {
        bufferSize=0;
        bufferStart=filePos;
        reading=false;
    }

    public void switchToReadBuffer() throws IOException
    {
        flush();
        reading=true;
    }

    private void updateReadBuffer() throws IOException
    {
        super.seek(filePos);
        bufferStart=filePos;
        int n = super.read(buffer, 0, buffer.length);
        if (n < 0) n=0;
        bufferSize = n;
    }

    @Override
    public long getFilePointer() throws IOException 
    {
        return filePos;
    }

    @Override
    public long length() throws IOException 
    {
        return fileLength;
    }

    @Override
    public void seek(long pos) throws IOException
    {
        filePos=pos;
        if (filePos>fileLength) filePos=fileLength;
        if (filePos<0) filePos=0;
    }

    public void flush() throws IOException
    {   
        if (reading) return;
        super.seek(bufferStart);
        super.write(buffer,0,bufferSize);
        bufferSize=0;
    }

    @Override
    public void setLength(long newLength) throws IOException 
    {
        flush();
        super.setLength(newLength);
        fileLength=newLength;
        seek(filePos);
    }

    @Override
    public void close() throws IOException 
    {
        flush();
        super.close();
    }
}
