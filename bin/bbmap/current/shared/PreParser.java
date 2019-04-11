package shared;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile1;

/**
 * Pre-parses the command line to handle:
 * Leading hyphens
 * Java flags
 * Output stream redirection
 * Help/version information
 * 
 * @author Brian Bushnell
 * @date November 28, 2017
 *
 */
public class PreParser {
	
	public PreParser(final String[] args, final Class<?> c){
		this(args, System.err, c, false);
	}
	
	public PreParser(final String[] args, final Class<?> c, boolean printVersion){
		this(args, System.err, c, printVersion, true, true);
	}
	
	public PreParser(final String[] args, PrintStream defaultPrintStream, final Class<?> c, boolean printVersion){
		this(args, defaultPrintStream, c, printVersion, true, true);
	}
	
	public PreParser(String[] args0, PrintStream defaultPrintStream, final Class<?> c, boolean printVersion, boolean removeKnown, boolean autoExit){
		original=args0;
		
		Parser.parseHelp(args0, autoExit);
		
		String[] args=Parser.parseConfig(args0);
		this.config=(args!=args0);
		
		final ArrayList<String> unknown=new ArrayList<String>(args.length);
		int removed=0, hyphens=0;

		PrintStream outstream=((defaultPrintStream == null) ? System.err : defaultPrintStream);
		boolean help=false, jflag=false;
		
		for(int i=0; i<args.length; i++){
			String s=args[i];
			boolean remove=false;
			if(Parser.isJavaFlag(s)){
				jflag=true;
				remove=true;
			}else if(s==null){
				remove=true;
			}else{
				
//				while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);} //Another possibility
				
				if(s.length()>0 && s.charAt(0)=='-'){
					int cnt=1;
					while(cnt<s.length() && s.charAt(cnt)=='-'){cnt++;}
					s=s.substring(cnt);
					hyphens+=cnt;
					args[i]=s;//Note that this mutates the original.  Moral:  Don't use hyphens.
				}
				
				final String[] split=s.split("=");
				assert(split.length<3) : "To many '=' signs: "+s;
				final String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b==null){
					//do nothing
				}else if("null".equalsIgnoreCase(b)){//Replace "null" with nothing
					b=null;
					args[i]=s=(a+"=");
				}
				
				if(a.equals("outstream")){
					outstream=parseOutstream(b);
					remove=removeKnown;
				}else if(a.equals("printexecuting")){
					printExecuting=Tools.parseBoolean(b);
					remove=true;
				}else if(a.equals("bufferbf") || a.equals("bufferbf1")){//for testing
					ByteFile1.BUFFERED=Tools.parseBoolean(b);
					remove=true;
				}
			}
			
			if(remove){removed++;}
			else{unknown.add(s);}
		}
		this.args=(removed==0 ? args : unknown.toArray(new String[0]));
		this.outstream=outstream;
		this.help=help;
		this.jflag=jflag;
		this.hyphens=hyphens;
		
		if(!silent){
			if(c!=null && printExecuting){outstream.println("Executing "+c.getName()+" "+Arrays.toString(original));}
			if(printVersion){outstream.println("Version "+Shared.BBMAP_VERSION_STRING/*+" "+Arrays.toString(args)*/);}
			if(c!=null || printVersion){outstream.println();}
		}
	}

	public static int stripHyphens(String[] args){
		int stripped=0;
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			if(Parser.isJavaFlag(arg)){
				//Do nothing
			}else if(arg.length()>0 && arg.charAt(0)=='-'){
				int cnt=1;
				while(cnt<arg.length() && arg.charAt(cnt)=='-'){cnt++;}
				arg=arg.substring(cnt);
				stripped+=cnt;
			}
		}
		return stripped;
	}
	
	private static PrintStream parseOutstream(final String b) {
		
		try {
			if(b==null || b.equalsIgnoreCase("null")){return new PrintStream(new NullOutputStream());}
			else if(b.equalsIgnoreCase("stdout") || b.startsWith("stdout.") || b.equals("System.out") || b.equalsIgnoreCase("sysout")){return System.out;}
			else if(b.equalsIgnoreCase("stderr") || b.startsWith("stderr.") || b.equals("System.err")){return System.err;}
			else{return new PrintStream(b);}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		KillSwitch.kill("Unable to process argument outstream="+b);
		return null; //Unreachable
	}

	/** Writes to nowhere.
	 * Courtesy of https://stackoverflow.com/a/692580 and https://stackoverflow.com/a/691835 */
	public static class NullOutputStream extends OutputStream {
		
		@Override
		public void write(int b) throws IOException {}
		
		@Override
		public void write(byte[] b) throws IOException {}
		
		@Override
		public void write(byte[] b, int off, int len) throws IOException {}
		
	}
	
	public final String[] original;
	public final String[] args;
	public final PrintStream outstream;
	public final boolean help;
	public final boolean config;
	public final boolean jflag;
	public final int hyphens;

	public static boolean printExecuting=true;
	public static boolean silent=false;
	
}
