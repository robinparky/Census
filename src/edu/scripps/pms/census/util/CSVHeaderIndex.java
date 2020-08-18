package edu.scripps.pms.census.util;


public  class CSVHeaderIndex {

    public final String headerName;
    private  int index =-1;

    public CSVHeaderIndex(String headerName) {
        this.headerName = headerName;
    }

    public boolean checkSetIndex(String candidate, int index)
    {
        if(this.index  == -1 && candidate.equals(headerName))
        {
            this.index = index;
            return true;
        }
        else if(this.index !=-1)
        {
            return true;
        }
        return false;
    }

    public void reset()
    {
        index = -1;
    }


    public int getIndex() {
        return index;
    }
}
