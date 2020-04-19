import java.io.IOException;
import org.apache.hadoop.fs.*;
import java.io.*;
import java.util.*;
import java.text.ParseException;
import java.net.URI;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.mapreduce.filecache.DistributedCache;
import java.io.File;
import org.apache.commons.lang.WordUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.*;
import org.apache.hadoop.mapreduce.lib.output.*;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.lang.*;
import java.lang.Object.*;

public class Covid19_3{
    // 4 types declared: Type of input key, type of input value, type of output key, type of output value
    public static class MyMapper extends Mapper<Object, Text, Text, DoubleWritable> {
        private static DoubleWritable one = new DoubleWritable(0);
        public static Map<String, String> hashmap;
        private Text word = new Text();
        public int line = 0;
        public static void readCacheFile(Path cacheFilePath) throws IOException {
            hashmap = new HashMap<String, String>();
            int n = 0;
            BufferedReader reader = new BufferedReader(new FileReader(cacheFilePath.toUri().getPath()));
            String line;
            while ((line = reader.readLine()) != null) {
                String[] words = line.split(",");
                if(words.length > 4){
                    hashmap.put(words[1],words[4]);
                }

            }
            reader.close();
        }
        // The 4 types declared here should match the types that was declared on the top
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            StringTokenizer tok = new StringTokenizer(value.toString(), ";\"\".\t,");
            Configuration conf = context.getConfiguration();
            Path[] localCacheFiles = DistributedCache.getLocalCacheFiles(conf);
            readCacheFile(localCacheFiles[0]);
            //URI[] url = context.getCacheFiles();
            int n = tok.countTokens();
           while (tok.hasMoreTokens()) {
                line++;

                if(line == 1){
                    break;
                }
                if(n == 4){
                    try {
                        SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
                        Date dt1 = format.parse(tok.nextToken());
                        Date dt2 = format.parse("2020-01-01");
                        if(dt1.after(dt2) == false){
                            break;
                        }
                    } catch (ParseException e) {
                        e.printStackTrace();
                    }
                }
                else if( n == 3 ){
                    word = new Text(tok.nextToken());
                    /*if(word.toString().equals("World") && conf.get("t").equals("false")){
                        break;
                    }*/
                }else if( n == 2 ){
                    double a = Double.parseDouble(tok.nextToken());
                    if(hashmap.get(word.toString()) == null){
                        break;
                    }
                    double pop = Double.parseDouble(hashmap.get(word.toString()));

                    a = a/pop*1000000;
                    one = new DoubleWritable(a);

                    context.write(word, one);
                }else{
                    tok.nextToken();
                }
                n--;
            }
        }

    }

    // 4 types declared: Type of input key, type of input value, type of output key, type of output value
    // The input types of reduce should match the output type of map
    public static class MyReducer extends Reducer<Text, DoubleWritable, Text, DoubleWritable> {
        private DoubleWritable total = new DoubleWritable();

        // Notice the that 2nd argument: type of the input value is an Iterable collection of objects
        //  with the same type declared above/as the type of output value from map
        public void reduce(Text key, Iterable<DoubleWritable> values, Context context) throws IOException, InterruptedException {
            double sum = 0;
            for (DoubleWritable tmp: values) {
                sum += tmp.get();
            }
            total.set(sum);
            // This write to the final output
            context.write(key, total);
        }
    }


    public static void main(String[] args)  throws Exception {
        Configuration conf = new Configuration();
        DistributedCache.addCacheFile(new URI(args[1]), conf);
        Job myjob = Job.getInstance(conf, "my word count test");
        //myjob.addCacheFile(new Path(args[1]).toUri());
        myjob.setJarByClass(Covid19_3.class);
        myjob.setMapperClass(MyMapper.class);
        myjob.setReducerClass(MyReducer.class);
        myjob.setOutputKeyClass(Text.class);
        myjob.setOutputValueClass(DoubleWritable.class);
        FileInputFormat.addInputPath(myjob, new Path(args[0]));
        FileOutputFormat.setOutputPath(myjob,  new Path(args[2]));
        System.exit(myjob.waitForCompletion(true) ? 0 : 1);
    }
}
