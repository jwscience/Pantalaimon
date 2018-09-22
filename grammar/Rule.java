package grammar;

public class Rule{
	private Nonterminal leftside;
	private GrammarSymbol[] rightside;
	private String string;
	private double probability;

	public Rule(Nonterminal leftside, GrammarSymbol[] rightside, double probability){
		this.leftside = leftside;
		this.rightside = rightside;
		this.probability = probability;

		// TODO: use string builder or something
		this.string = leftside.toString()+ " -> ";
		for(GrammarSymbol r : rightside){
			string += r.toString() + " ";
		}
	}

	public Nonterminal getLeftside(){
		return leftside;
	}

	public int length(){
		return rightside.length;
	}

	public GrammarSymbol rightside_char(int n){
		return rightside[n];
	}

	public double getProbability(){
		return probability;
	}

	public String toString(){
		return string;
	}
}
