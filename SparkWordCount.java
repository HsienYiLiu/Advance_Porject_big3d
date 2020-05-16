/* Java imports */
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.lang.Iterable;
import java.util.Iterator;
import java.io.*; //import
/* Spark imports */
import scala.Tuple2;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.PairFlatMapFunction;
public class SparkWordCount {

    public static int line = 0;
    /**
     * args[0]: Input file path on distributed file system
     * args[1]: Output file path on distributed file system
     */
    public static void main(String[] args){
	System.out.println("Hello World from Java");

	String input = args[0];
        
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	System.out.println("Hello World from Java");
	//String start = args[1];
	//String end = args[2];
	String output = args[1];
	
        
	/* essential to run any spark code */
	SparkConf conf = new SparkConf().setAppName("WordCount");
	JavaSparkContext sc = new JavaSparkContext(conf);

	/* load input data to RDD */
	JavaRDD<String> dataRDD = sc.textFile(input);

	JavaPairRDD<String, Integer> counts =
	    dataRDD.flatMapToPair(new PairFlatMapFunction<String, String, Integer>(){
		    public Iterator<Tuple2<String, Integer>> call(String value){
			line++;
			if(line ==1 ){
			    return null;
			}
			String[] words = value.split(",");
			
			List<Tuple2<String, Integer>> retWords =
		     	new ArrayList<Tuple2<String, Integer>>();
			
			
			retWords.add(new Tuple2<String, Integer>(words[0], Integer.parseInt(words[3])));
			 
			
			return retWords.iterator();
		    }
		}).reduceByKey(new Function2<Integer, Integer, Integer>(){
			public Integer call(Integer x, Integer y){
			    return x+y;
			}
		    });
	
	counts.saveAsTextFile(output);
	
    }
}
