package math;

public class Probability{
	public static final int impossible = -100000; // e^(-1000) = 1 * 10^(-4343)
	public static final int sure = 0;
	private static final int neglect_low = -80000;     // e^(-500) = 4 * 10^(-3475)

	public static double fromNumber(double number){
		if(number <= 0) return impossible;
		return Math.log(number);
	}

	public static double toNumber(double probability){
		if(probability <= neglect_low) return 0;
		return Math.exp(probability);
	}

	public static double fromRatio(double numerator, double denominator){
		if(numerator == 0) return impossible;
		if(numerator >= denominator) return sure;
		return Math.log(numerator) - Math.log(denominator);
	}

	public static double multiply(double p1, double p2){
		if(p1 <= neglect_low) return impossible;
		if(p2 <= neglect_low) return impossible;
		double result = p1 + p2;
		if(result > sure) result = sure;
		return result;
	}

	public static double multiply(double p1, double p2, double p3){
		return multiply(p1, multiply(p2, p3));
	}
	
	public static double multiply(double p1, double p2, double p3, double p4){
		return multiply(p1, multiply(p2, multiply(p3, p4)));
	}

	public static double multiply(double p1, double p2, double p3, double p4, double p5){
		return multiply(p1, multiply(p2, multiply(p3, multiply(p4, p5))));
	}

	public static double divide(double p1, double p2){
		if(p1 <= neglect_low) return impossible;
		if(p1 >= p2) return sure;
		return p1 - p2;
	}

	public static double puredivide(double p1, double p2){
		// The correction factors for {alpha|beta}_hat seem to need this since they
		// are not really probabilities and are (maybe) desired to be greater than 1
		return p1 - p2;
	}

	public static double add(double p1, double p2){
		if(p1 <= neglect_low && p2 <= neglect_low) return impossible;
		if(p1 <= neglect_low) return p2;
		if(p2 <= neglect_low) return p1;

		return p1 + Math.log(1 + Math.exp(p2 - p1));
	}

	public static double power(double p, int n){
		return n * p;
	}

	public static double max(double p1, double p2){
		// log is monotonic
		if(p1 < p2) return p2;
		return p1;
	}

	public static boolean greater(double p1, double p2){
		return p2 > p1;
	}

	public static boolean isPossible(double p1){
		return p1 > neglect_low;
	}
}
