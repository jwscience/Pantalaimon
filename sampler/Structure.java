package sampler;


public class Structure{
	public static final int hairpinloop = 1;
	public static final int bulgeinterior = 2;
	public static final int stackedpair = 3;
	public static final int multiloop = 4;
	public static final int exteriorloop = 5;

	public static final int hl_no = 1;

	public static final int bp_a = 2;
	public static final int bp_p = 3;
	public static final int bp_special = 4;

	public static final int bi_ba = 5;
	public static final int bi_ab = 6;
	public static final int bi_bab = 7;

	public static final int ml_1 = 8;
	public static final int ml_2 = 9;
	public static final int ml_3 = 10;

	public static final int el_c = 11;
	public static final int el_ca = 12;
	public static final int el_cat = 13;
	public static final int el_a = 14;
	public static final int el_at = 15;

	private int start;
	private int end;
	private int type;
	private int subtype;
	private double probability;

	public Structure(int start, int end, double probability, int type, int subtype){
		this.type = type;
		this.subtype = subtype;
		this.start = start;
		this.end = end;
		this.probability = probability;
	}

	public int getType(){
		return type;
	}

	public int getSubtype(){
		return subtype;
	}

	public int getStart(){
		return start;
	}

	public int getEnd(){
		return end;
	}

	public double getProbability(){
		return probability;
	}
}
