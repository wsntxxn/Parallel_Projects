import java.io.IOException;
import java.util.StringTokenizer;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.FloatWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.FileInputFormat;
import org.apache.hadoop.mapreduce.lib.output.FileOutputFormat;

public class MinMax {
    public static class MinMaxMapper extends Mapper<Object, Text, Text, FloatWritable> {
        //private final static FloatWritable data = new FloatWritable(0);
        private Text keyForReduce = new Text(); 

        private float maxValue = Float.MIN_VALUE;
        private float minValue = Float.MAX_VALUE;

        public void map(Object key, Text value, Context context) throws IOException, InterruptedException {
            String line = value.toString();
            String[] tempdata = line.split(",");
            if (tempdata[2].contains("Temperature"))
                return;

            float data = Float.parseFloat(tempdata[2]);
            if (data > maxValue)
                maxValue = data;
            if (data < minValue)
                minValue = data;
            //data.set(Float.parseFloat(tempdata[2]));
        }

        public void cleanup(Context context) throws IOException, InterruptedException {
            context.write(new Text("Max"), new FloatWritable(maxValue));
            context.write(new Text("Min"), new FloatWritable(minValue));
        }
    }

    public static class MinMaxReducer extends Reducer<Text, FloatWritable, Text, FloatWritable> {
        
        private float maxValue = Float.MIN_VALUE;
        private float minValue = Float.MAX_VALUE;

        public void reduce(Text key, Iterable<FloatWritable> values, Context context) throws IOException, InterruptedException {
            if (key.toString().equals("Max")) {
                for (FloatWritable val: values) {
                    if (val.get() > maxValue)
                        maxValue = val.get();
                }
                context.write(new Text("MaxValue"), new FloatWritable(maxValue));
            }
            else {
                for (FloatWritable val: values) {
                    if (val.get() < minValue)
                        minValue = val.get();
                }
                context.write(new Text("MinValue"), new FloatWritable(minValue));
            }
        }
    }

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration();
        Job job = Job.getInstance(conf, "min max");

        job.setJarByClass(MinMax.class);

        job.setMapperClass(MinMaxMapper.class);
        job.setMapOutputKeyClass(Text.class);
        job.setMapOutputValueClass(FloatWritable.class);

        job.setReducerClass(MinMaxReducer.class);
        job.setOutputKeyClass(Text.class);
        job.setOutputValueClass(FloatWritable.class);

        FileInputFormat.addInputPath(job, new Path(args[0]));
        FileOutputFormat.setOutputPath(job, new Path(args[1]));

        System.exit(job.waitForCompletion(true) ? 0 : 1);
    }
}

