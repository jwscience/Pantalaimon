package biology;

import math.Statistics;

import java.util.Arrays;
import java.util.Stack;

public class SecondarySequence{
	private Secondary[] sequence;
	private int[] brackets;

	public SecondarySequence(int length){
		sequence = new Secondary[length];
		brackets = new int[length];
		Arrays.fill(brackets, -1);
	}

	public SecondarySequence(String seq){
		int l = seq.length();
		sequence = new Secondary[l];
		brackets = new int[l];
		Stack<Integer> s = new Stack<>();

		for(int i=0; i<l; i++){
			char c = seq.charAt(i);
			if(c == '.' || c == '|'){
				sequence[i] = Secondary.unpaired;
				brackets[i] = -1;
			}else if(c == '('){
				sequence[i] = Secondary.open;
				s.push(i);
			}else if(c == ')'){
				sequence[i] = Secondary.close;
				int j = s.pop();
				brackets[j] = i;
				brackets[i] = j;
			}else if(c == '*'){
				sequence[i] = Secondary.undecided;
				brackets[i] = -1;
			}else{
				sequence[i] = null;
				brackets[i] = -1;
			}
		}
	}

	public SecondarySequence(Secondary[] sec, int[] br){
		sequence = sec;
		brackets = br;
	}

	public Secondary get(int pos){
		return sequence[pos];
	}

	public int getMatchingBracket(int pos){
		return brackets[pos];
	}

	public void setBase(int pos, boolean unpaired){
		if(unpaired){
			sequence[pos] = Secondary.unpaired;
		}else{
			sequence[pos] = Secondary.undecided;
		}
	}

	public void setBasepair(int pos1, int pos2, boolean paired){
		if(paired){
			sequence[pos1] = Secondary.open;
			sequence[pos2] = Secondary.close;
			brackets[pos1] = pos2;
			brackets[pos2] = pos1;
		}else{
			sequence[pos1] = Secondary.undecided;
			sequence[pos2] = Secondary.undecided;
			brackets[pos1] = -1;
			brackets[pos2] = -1;
		}
	}

	public int getLength(){
		return sequence.length;
	}

	public boolean isBasepair(int pos1, int pos2){
		return get(pos1) == Secondary.open && getMatchingBracket(pos1) == pos2;
	}

	public Statistics compare(SecondarySequence other){
		//I am the predicted and the other is the actual solution
		int tp = 0, fp = 0, fn = 0;
		for(int i=0; i<getLength(); i++){
			if(get(i) == Secondary.open){
				int j = getMatchingBracket(i);
				if(other.isBasepair(i, j)){
					tp += 1;
				}else{
					fp += 1;
				}
			}
			if(other.get(i) == Secondary.open){
				int j = other.getMatchingBracket(i);
				if(! isBasepair(i, j)){
					fn += 1;
				}
			}
		}
		return new Statistics(tp, fp, fn);
	}

	public String toString(){
		String result = "";
		for(Secondary c : sequence){
			result += c.toString();
		}
		return result;
	}
}
