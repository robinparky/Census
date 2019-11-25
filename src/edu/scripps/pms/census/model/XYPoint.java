package edu.scripps.pms.census.model;

/**
 * Created by rampuria on 1/31/17.
 */
public class XYPoint {

	private double x;

	private double y;

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}
	public XYPoint()
	{}


	public XYPoint(edu.scripps.pms.census.util.XYPoint xyPoint)
	{
		x = xyPoint.getX();
		y = xyPoint.getY();
	}
}
