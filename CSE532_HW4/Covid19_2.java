import java.io.IOException;
import java.util.*;
import java.text.ParseException;
import org.apache.commons.lang.WordUtils;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.*;
import org.apache.hadoop.mapreduce.*;
import org.apache.hadoop.mapreduce.lib.input.*;
import org.apache.hadoop.mapreduce.lib.output.*;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Covid19_2 {
    // 4 types declared: Type of input key, type of input value, type of output key, type of output value
    public static class MyMapper extends Mapper<Object, Text, Text, LongWritable> {
        private static LongWritable one = new LongWritable(1);
        private Text word = new Text();
        public int line = 0;

        // The 4 types declared here should match the types that was declared on the top
        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            StringTokenizer tok = new StringTokenizer(value.toString(), ";\"\".\t,");
            Configuration conf = context.getConfiguration();
            int n = tok.countTokens();
            while (tok.hasMoreTokens()) {
                line++;
                if(line == 1){
                    break;
                }
                if(n == 4){
                    try {
                        SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
                        Date data_dt = format.parse(tok.nextToken());
                        Date dt1 = format.parse(conf.get("start"));
                        Date dt2 = format.parse(conf.get("end"));
                        if(data_dt.after(dt1) == false || data_dt.after(dt2) == true){
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
                }else if( n == 1 ){
                    one = new LongWritable(Long.parseLong(tok.nextToken()));
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
    public static class MyReducer extends Reducer<Text, LongWritable, Text, LongWritable> {
        private LongWritable total = new LongWritable();

        // Notice the that 2nd argument: type of the input value is an Iterable collection of objects
        //  with the same type declared above/as the type of output value from map
        public void reduce(Text key, Iterable<LongWritable> values, Context context) throws IOException, InterruptedException {
            long sum = 0;
            for (LongWritable tmp: values) {
                sum += tmp.get();
            }
            total.set(sum);
            // This write to the final output
            context.write(key, total);
        }
    }


    public static void main(String[] args)  throws Exception {
        Configuration conf = new Configuration();
        SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd");
        Date start = format.parse(args[1]);
        Date end = format.parse(args[2]);
        Date dt1 = format.parse("2019-12-31");
        Date dt2 = format.parse("2020-04-08");
        //start.after(dt1) == false
        if(start.after(dt1) == false || end.after(dt2) == true){
            System.out.println("Check your date, start from 2019-12-31 and end before 2020-04-08");
            System.out.println("System will exit now");
            System.exit(1);
        }else{
            conf.set("start",format.format(start));
            conf.set("end", format.format(end));
        }
        Job myjob = Job.getInstance(conf, "my word count test");
        myjob.setJarByClass(Covid19_2.class);
        myjob.setMapperClass(MyMapper.class);
        myjob.setReducerClass(MyReducer.class);
        myjob.setOutputKeyClass(Text.class);
        myjob.setOutputValueClass(LongWritable.class);
        // Uncomment to set the number of reduce tasks
        // myjob.setNumReduceTasks(2);
        FileInputFormat.addInputPath(myjob, new Path(args[0]));
        FileOutputFormat.setOutputPath(myjob,  new Path(args[3]));
        System.exit(myjob.waitForCompletion(true) ? 0 : 1);
    }
}

