package sampler;

public class SampleStatus{
	private int ii;
	private int jj;
	private boolean completed;
	private int k;

	public SampleStatus(int ii, int jj, boolean completed){
		this.ii = ii;
		this.jj = jj;
		this.completed = completed;
		this.k = -1;
	}

	public SampleStatus(int ii, int jj, boolean completed, int k){
		this(ii, jj, completed);
		this.k = k;
	}

	public int getII(){
		return ii;
	}

	public int getJJ(){
		return jj;
	}

	public boolean isCompleted(){
		return completed;
	}

	public int getK(){
		return k;
	}
}
