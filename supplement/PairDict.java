package supplement;

import java.util.Arrays;

public class PairDict{
	private int length;
	private boolean[][] matrix;
	private int[] starts;
	private int[] ends;

	public PairDict(int length){
		this.length = length;
		matrix = new boolean[length][length];
		for(boolean[] i : matrix){
			Arrays.fill(i, false);
		}
		starts = new int[length];
		Arrays.fill(starts, length);
		ends = new int[length];
		Arrays.fill(ends, -1);
	}

	public void add(int start, int end){
		matrix[start][end] = true;
		if(start < starts[end]) starts[end] = start;
		if(end > ends[start]) ends[start] = end;
	}

	public boolean exists(int start, int end){
		return matrix[start][end];
	}

	public boolean existsEnd(int end){
		return starts[end] != length;
	}

	public boolean existsStart(int start){
		return ends[start] != -1;
	}

	public int minStart(int end){
		return starts[end];
	}

	public int maxEnd(int start){
		return ends[start];
	}

	public void output(){
		for(int i=0; i<matrix.length; i++){
			for(int j=0; j<matrix[i].length; j++){
				if(matrix[i][j]){
					System.out.println("Structure at ("+i+", "+j+")");
				}
			}
		}
	}
}
