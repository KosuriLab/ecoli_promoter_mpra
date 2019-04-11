package sketch;

import shared.Tools;
import structures.AbstractBitSet;

public class CompareBuffer {
	
	public CompareBuffer(boolean makeBS){
		if(makeBS){
			cbs=AbstractBitSet.make(0, SketchObject.bitSetBits);
		}else{
			cbs=null;
		}
	}
	
	void set(final int hits_, final int multiHits_, final int unique2_, final int unique3_, final int noHits_,
			final int contamHits_, final int contam2Hits_, final int multiContamHits_,
			final int queryDivisor_, final int refDivisor_, final int querySize_, final int refSize_, final long depthSum_, final double depthSum2_){
		hits=hits_;
		multiHits=multiHits_;
		unique2=unique2_;
		unique3=unique3_;
		noHits=noHits_;
		
		contamHits=contamHits_;
		contam2Hits=contam2Hits_;
		multiContamHits=multiContamHits_;
		
		queryDivisor=queryDivisor_;
		refDivisor=refDivisor_;
		
		querySize=querySize_;
		refSize=refSize_;

		depthSum=depthSum_;
		depthSum2=(float)depthSum2_;
	}
	
	void clear(){
		hits=multiHits=0;
		unique2=unique3=noHits=0;
		contamHits=contam2Hits=multiContamHits=0;
		refDivisor=queryDivisor=0;
		refSize=querySize=0;
		depthSum=0;
		depthSum2=0;
	}
	
	float depth(){
		return depthSum<1 ? 0 : depthSum/Tools.max(1.0f, hits);
	}
	
	float depth2(){
		return depthSum2<=0 ? 0 : depthSum2/Tools.max(1.0f, hits);
	}

	//For WKID
	int minDivisor(){return Tools.max(1, Tools.min(queryDivisor, refDivisor));}
	//For KID
	int maxDivisor(){return Tools.max(1, queryDivisor, refDivisor);}
	int minSize(){return Tools.max(1, Tools.min(querySize, refSize));}
	int maxSize(){return Tools.max(1, querySize, refSize);}

	int uniqueHits(){return hits-multiHits;}
	int uniqueContamHits(){return contamHits-multiContamHits;}
	
	/*--------------------------------------------------------------*/
	
	@Override
	public String toString(){
		return "hits="+hits+", refDivisor="+refDivisor+", queryDivisor="+queryDivisor+", refSize="+refSize+", querySize="+querySize+
				", contamHits="+contamHits+", contam2Hits="+contam2Hits+", multiContamHits="+multiContamHits+", depthSum="+depthSum+", depthSum2="+depthSum2+
				", hits="+hits+", multiHits="+multiHits+", unique2="+unique2+", unique3="+unique3+", noHits="+noHits;
	}
	
	/*--------------------------------------------------------------*/
	
	int hits(){return hits;}
	int multiHits(){return multiHits;}
	int noHits(){return noHits;}
	int unique2(){return unique2;}
	int unique3(){return unique3;}

	int contamHits(){return contamHits;}
	int contam2Hits(){return contam2Hits;}
	int multiContamHits(){return multiContamHits;}
	
	int queryDivisor(){return queryDivisor;}
	int refDivisor(){return refDivisor;}
	
	int querySize(){return querySize;}
	int refSize(){return refSize;}

	long depthSum(){return depthSum;}
	float depthSum2(){return depthSum2;}
	
	/*--------------------------------------------------------------*/
	
	private int hits;
	private int multiHits;
	private int noHits;
	private int unique2;
	private int unique3;

	private int contamHits;
	private int contam2Hits;
	private int multiContamHits;
	
	private int queryDivisor;
	private int refDivisor;
	
	private int querySize;
	private int refSize;

	private long depthSum;
	private float depthSum2;
	
	/*--------------------------------------------------------------*/

	public final AbstractBitSet cbs; //Only for comparisons, not index
	
}
