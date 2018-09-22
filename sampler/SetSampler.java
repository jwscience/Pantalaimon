package sampler;

import math.Probability;

import java.util.LinkedList;
import java.util.Random;

public class SetSampler{
	private LinkedList<Structure> list;
	private double probsum;
	private Random random;

	public SetSampler(long seed){
		random = new Random(seed);
	}

	public void add(Structure structure){
		list.addLast(structure);
		probsum = Probability.add(probsum, structure.getProbability());
	}

	public int size(){
		return list.size();
	}

	public void debug(){
		double cumulate = Probability.impossible;
		for(Structure structure : list){
			double plus = Probability.divide(structure.getProbability(), probsum);
			cumulate = Probability.add(cumulate, plus);
		}
		System.out.println("probsum is "+ Probability.toNumber(probsum) +", sum of probabilities is "+Probability.toNumber(cumulate));
	}

	public Structure uniformsample(){
		//System.out.println("Sampling from "+size()+" objects with a total probability of "+probsum);
		if(size() == 0) return null;

		double pivot = random.nextDouble();
		double cumulate = Probability.impossible;
		for(Structure structure : list){
			double plus = Probability.divide(structure.getProbability(), probsum);
			cumulate = Probability.add(cumulate, plus);
			if(Probability.toNumber(cumulate) >= pivot){
				return structure;
			}
		}

		assert false;
		return null;
	}

	public void restart(){
		list = new LinkedList<>();
		probsum = Probability.impossible;
	}
}
