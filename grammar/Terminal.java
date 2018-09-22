package grammar;

public class Terminal extends GrammarSymbol{
	private static int counter = 0;
	private int number;
	private char symbol;

	/**
	 * resets the counting ond enumerating of terminal symbols.
	 *
	 * This can be useful if more than one grammar is defined.
	 */
	public static void reset_counter(){
		counter = 0;
	}

	public Terminal(char c){
		symbol = c;
		number = counter;
		counter++;
	}

	public int getNumber(){
		return number;
	}

	public char getSymbol(){
		return symbol;
	}
	public String toString(){
		return Character.toString(symbol);
	}
}
