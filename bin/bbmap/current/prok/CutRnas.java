package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import stream.Read;
import stream.ReadInputStream;

public class CutRnas {
	
	public static void main(String[] args){
		
		Shared.TRIM_READ_COMMENTS=true;
		GffLine.parseAttributes=true;

		boolean invert=false;
		String fna=args[0];
		String gff=args[1];
		String out=args[2];
		String types=args[3];
		if(args.length>4){invert=Tools.parseBoolean(args[4]);}
		
		ArrayList<GffLine> lines=GffLine.loadGffFile(gff, types);
		
		ArrayList<Read> list=ReadInputStream.toReads(fna, FileFormat.FA, -1);
		HashMap<String, Read> map=new HashMap<String, Read>();
		for(Read r : list){map.put(r.id, r);}
		
		ByteStreamWriter bsw=new ByteStreamWriter(out, true, false, false);
		bsw.start();
		
		processStrand(lines, map, 0, bsw, invert);
		for(Read r : list){r.reverseComplement();}
		processStrand(lines, map, 1, bsw, invert);
		
		if(invert){
			for(Read r : list){
				r.reverseComplement();
				bsw.println(r);
			}
		}
		
		bsw.poisonAndWait();
	}
	
	public static void processStrand(ArrayList<GffLine> lines, HashMap<String, Read> map, int strand, ByteStreamWriter bsw, boolean invert){
		for(GffLine gline : lines){
			if(gline.strand==strand){
				Read scaf=map.get(gline.seqid);
				assert(scaf!=null) : "Can't find "+gline.seqid+" in "+map.keySet();
				int start, stop;
				if(strand==0){
					start=gline.start-1;
					stop=gline.stop-1;
				}else{
					start=scaf.length()-gline.stop-1;
					stop=scaf.length()-gline.start-1;
				}
				if(invert){
					byte[] bases=scaf.bases;
					for(int i=start; i<stop; i++){
						if(i>=0 && i<bases.length){
							bases[i]='N';
						}
					}
				}else{
					if(start>=0 && stop<scaf.length()){
						Read r=new Read(Arrays.copyOfRange(scaf.bases, start, stop), null, gline.attributes, 1);
						bsw.println(r);
					}
				}
			}
		}
	}
	
}
