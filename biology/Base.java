package biology;

public class Base{
	public static final int count = 4;

	public static Base adenine;
	public static Base cytosine;
	public static Base guanine;
	public static Base uracil;

	static{
		adenine = new Base('A', 0);
		cytosine = new Base('C', 1);
		guanine = new Base('G', 2);
		uracil = new Base('U', 3);
	}

	public static boolean validpair(Base u, Base v){
		return (u == adenine && v == uracil) ||
			(u == uracil && v == adenine) ||
			(u == cytosine && v == guanine) ||
			(u == guanine && v == cytosine) ||
			(u == guanine && v == uracil) ||
			(u == uracil && v == guanine);
	}

	public static boolean validpair(int u, int v){
		return (u == 0 && v == 3) ||
			(u == 3 && v == 0) ||
			(u == 1 && v == 2) ||
			(u == 2 && v == 1) ||
			(u == 2 && v == 3) ||
			(u == 3 && v == 2);
	}

	private char symbol;
	private int id;

	private Base(char symbol, int id){
		this.symbol = symbol;
		this.id = id;
	}

	public char getSymbol(){
		return symbol;
	}

	public String toString(){
		return Character.toString(symbol);
	}

	public int getId(){
		return id;
	}
}
