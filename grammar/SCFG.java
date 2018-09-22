package grammar;

import java.util.Iterator;
import java.util.LinkedList;

public class SCFG{
	private Rule[] rules;
	private LinkedList<Nonterminal> nonterminals;
	private LinkedList<Terminal> terminals;
	private LinkedList<Integer>[] ntlookup;


	public SCFG(int nonterminal_count, Terminal[] terminalsymbols, Rule[] rules){ //, Rule[] ntrules){
		this.rules = new Rule[rules.length];
		System.arraycopy(rules, 0, this.rules, 0, rules.length);

		nonterminals = new LinkedList<>();
		terminals = new LinkedList<>();
		ntlookup = new LinkedList[nonterminal_count];

		for(Terminal t : terminalsymbols){
			terminals.addLast(t);
		}
		for(int i=0; i<rules.length; i++){
			Nonterminal l = rules[i].getLeftside();
			int index = l.getNumber();
			if(ntlookup[index] == null){
				// means we haven't seen this nonterminal before
				nonterminals.addLast(l);
				ntlookup[index] = new LinkedList<>();
			}
			ntlookup[index].addLast(i);
		}
	}

	public Iterator<Integer> getRules(Nonterminal nt){
		return ntlookup[nt.getNumber()].iterator();
	}

	public Rule getRule(int index){
		return rules[index];
	}

	public int getNonterminalCount(){
		return nonterminals.size();
	}

	public Iterator<Nonterminal> getNonterminals(){
		return nonterminals.iterator();
	}

	public Iterator<Terminal> getTerminals(){
		return terminals.iterator();
	}

	public int getTerminalCount(){
		return terminals.size();
	}

	public int getRuleCount(){
		return rules.length;
	}
}
