package grammar;

public class Nonterminal extends GrammarSymbol{
	private static int counter = 0;
	private int number;
	private char symbol;

	public Nonterminal(char symbol){
		this.number = counter;
		counter++;
		this.symbol = symbol;
	}

	public int getNumber(){
		return number;
	}

	public String toString(){
		return Character.toString(symbol);
	}
}
