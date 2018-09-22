package math;

public class Statistics{
	private int tp;
	private int fp;
	private int fn;

	public Statistics(int tp, int fp, int fn){
		this.tp = tp;
		this.fp = fp;
		this.fn = fn;
	}

	public double sensitivity(){
		if(tp == 0 && fn == 0) return 0;
		return ((double) tp) / (tp + fn);
	}

	public double ppv(){
		if(tp == 0 && fp == 0) return 0;
		return ((double) tp) / (tp + fp);
	}
}
