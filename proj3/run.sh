echo ${JAVA_HOME}
export HADOOP_CLASSPATH=${JAVA_HOME}/lib/tools.jar
hadoop com.sun.tools.javac.Main MinMax.java
jar cf minmax.jar MinMax*.class
hadoop fs -mkdir /input
hadoop fs -put data/* /input
hadoop jar minmax.jar MinMax /input /output
