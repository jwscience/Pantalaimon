package biology;

public class Secondary{
	public static Secondary undecided;
	public static Secondary unpaired;
	public static Secondary open;
	public static Secondary close;

	static{
		undecided = new Secondary('*', 0);
		unpaired = new Secondary('.', 1);
		open = new Secondary('(', 2);
		close = new Secondary(')', 3);
	}


	private char symbol;
	private int id;

	private Secondary(char symbol, int id){
		this.symbol = symbol;
		this.id = id;
	}

	public char getSymbol(){
		return symbol;
	}

	public int getId(){
		return id;
	}

	public String toString(){
		return Character.toString(getSymbol());
	}
}
