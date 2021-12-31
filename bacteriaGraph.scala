// imports and contexts

import org.apache.spark._
import org.apache.spark.graphx._
import org.apache.spark.rdd.RDD
import scala.collection.mutable.ListBuffer
import scala.collection.mutable.ArrayBuffer
import org.apache.spark.sql.functions.sqrt
import org.apache.spark.sql.functions.col
import org.apache.spark.sql.functions.pow
import org.apache.spark.sql.DataFrame
import java.util.Calendar
import scala.reflect.io.Directory
import java.io._


val sqlContext = new org.apache.spark.sql.SQLContext(sc)
import sqlContext.implicits._


/* 

deatP, reprP : Parameters that rule replication, inputs for binomial function.
               Approx. (deathP)*(strain numerosity) bacteria reproduce or die for each iteration.

deadlyInter, reproductiveInter : strain dictionaries whose keys are interaction subject.
                                 Values represents strenght of interaction (0<=x<=1)

*/

case class Strain(var deathP:Double = 0.0, 
                  var reprP:Double = 0.0,  
                  var deadlyInter: Map[Long, Double] = Map().withDefaultValue(0.0),
                  var reproductiveInter: Map[Long, Double] = Map().withDefaultValue(0.0),
                 )

/*
distribution : Population array for each strain.
                **NB. makeConsistent() method fills the array by default according to number
                      of strains specified in Strains (front-end parameters).
                      Numerosity of each strain is randomly choosen according to range_popSize setting 
                      parameter.
*/
case class Node(var distribution:Array[Long] = Array[Long](),
                var location:Point = new Point(-1,-1),
                var popIn:Long = -1){
    if(distribution.size!=0) popIn = distribution.sum
}

// dist() method : compute euclidean distance between points

case class Point(value_X: Int, value_Y: Int){
    def dist(point_2:Point): Double={
        return math.sqrt(math.pow(point_2.value_X-value_X, 2)+math.pow(point_2.value_Y-value_Y, 2))
    }
}

case class settings(
                var pop:Array[Node] = Array[Node](),
                var strains:Array[Strain] = Array[Strain](),
                var lastSnapshots:Tuple2[DataFrame, DataFrame] = null,
                var range_popSize:Array[Int] = Array(1000, 5*1000),
                var dim_space:Int = 100,
                var maxEdgeDistance:Int = 20,
                var minVerticesDistance:Int = 10,
                var maxTransf:Double = 0.05,
                var transfP:Double = 0.9,
                var nIterations:Int = 40,
                var ceiling:Long = 1000000,
                var halving:Long = 50,
                var expName:String = "Experiment",
                var set_reallyCasual:Boolean = true,
                var in_reallyCasual:Boolean = true,
                var ev_reallyCasual:Boolean = true,
                ){
    var r_set = scala.util.Random
    var r_in = scala.util.Random
    var r_ev = scala.util.Random
    if(!set_reallyCasual) r_set.setSeed(1001L)
    if(!in_reallyCasual) r_in.setSeed(1002L)
    if(!ev_reallyCasual) r_ev.setSeed(1003L)
    
    // "pop" variable is an Array of type Node, which has strain distributions and locations input parameters
    // if we do not specify any of these, getCasualNodes fills randomnly "strain.size"-number of strains.
    //
    // here an example of initialization of pop:
    /*
        var sett = settings()
        sett.pop = Array(Node(distribution=Array(10,20), location=Point(5,5)))
        
        sett.pop(0).distribution -> Output -> Array(10,20)
        sett.pop(0).location -> Output -> Point(10,10)
        sett.pop(0).popIn -> Output -> -1 
        
    */    
    
     
    if(pop.size == 0) getCasualNodes(strains.size);
    def getCasualNodes(n:Int): Unit={
        pop = Array.fill(n)(Node(distribution = Array.fill(n)(0L)))
        for(i<-0 until n) pop(i).distribution(i) = r_set.nextInt(range_popSize(1)-range_popSize(0))+
                                                                range_popSize(0)
    }
    
    
    
    def makeConsistent(): Unit={
        // la variabile nStrains è di default pari alla cardinalità dell'Array "strains", variabile 
        // della classe "settings" e selezionabile dall'utente nell'interfaccia.
        // ** se l'array "strains" non viene riempito dall'utente la size di default è 0
        var nStrains = strains.size
        
        // questo ciclo riempie l'array "pop" nel caso non sia stato riempito dall'utente nell'interfaccia
        for(i<-0 until pop.size){
            if(pop(i).popIn== -1)
                pop(i).popIn = r_set.nextInt(range_popSize(1)-range_popSize(0)) + range_popSize(0)
        }
        if(nStrains==0) throw new Exception("Irregular settings")
        for(i<-0 until pop.size){
            var distribution = Array.fill(nStrains)(0L)
            for(j<-0 until nStrains) if(pop(i).distribution.size>j) distribution(j) = pop(i).distribution(j)
            pop(i).distribution = distribution
        }
        for(i<-0 until pop.size){
            if(pop(i).distribution.sum == 0)
            pop(i).distribution(r_set.nextInt(nStrains)) = pop(i).popIn
        }
    }
}

%%python
from matplotlib import pyplot as plt
import pandas as pd
import pathlib
import os.path as pth
from os import listdir
from pyspark.sql import SQLContext
sqlContext = SQLContext(sc)

def getDataframes(dir):
    folders = ["bacteriaFiles/"+dir+"/"+f for f in listdir("bacteriaFiles/"+dir)]
    Vdf = {}
    Edf = {}
    df = pd.DataFrame()
    VPositionDF = pd.DataFrame()
    EPositionDF = pd.DataFrame()
    for file in folders:
        parts = [s for s in listdir(file) if '.part' not in s and 'SUCCESS' not in s]
        paths = [file+"/"+part for part in parts]
        
        for i in range(len(paths)):
            try:
                df = pd.read_csv(paths[i], header=None)
                break
            except: None
        
        initialIdx = i
        colnames = [i for i in range(len(df.columns))]
        df.columns = colnames
        for path in paths[initialIdx+1:]:
            try:
                df = pd.concat([df, pd.read_csv(path, names=colnames, header=None)])
            except: None
        if 'Vsnap' in file:
            Vdf[file.split('/')[-1]] = df
        elif 'Esnap' in file:
            Edf[file.split('/')[-1]] = df
            

    Vsorted = sorted(Vdf.items())
    Esorted = sorted(Edf.items())
    return [k[1] for k in Vsorted], [k[1] for k in Esorted]

def save(df:org.apache.spark.sql.DataFrame, dir:String, name:String): Boolean={
    try{
        var buildName = name
        if(name.contains("Vsnap") || name.contains("Esnap")){
            buildName = name.substring(0,5)
            val suf = name.substring(5)
            for(i<-0 until 10 - suf.size) buildName += "0"
            buildName += suf
        }
        df.write.csv("bacteriaFiles/"+dir+"/"+buildName)
        return true
    }
    catch{case e: Exception => return false}
}

def load(dir:String, name:String): org.apache.spark.sql.DataFrame={
    try{
        return spark.read.format("csv").option("header", false).load("bacteriaFiles/"+dir+"/"+name)
    }
    catch{case e: Exception => return null}
}

def delete(dir:String, name:String): Boolean={
    val directory = new Directory(new File("bacteriaFiles/"+dir+"/"+name))
    return directory.deleteRecursively()
}

def exists(dir:String, name:String): Boolean={
    val file = new File("bacteriaFiles/"+dir+"/"+name)
    val directory = new Directory(file)
    return directory.exists
}

def getFiles(dir:String, onlyName:Boolean = true): List[String]={
    try{
        val file = new File("bacteriaFiles/"+dir)
        val innerDirs = ListBuffer[String]()
        if(onlyName) file.listFiles.toList.foreach(x => innerDirs += x.toString.split("/")(2))
        else file.listFiles.toList.foreach(x => innerDirs += x.toString.split("/")(0) + "/" +
                                                  x.toString.split("/")(1)+ "/" +x.toString.split("/")(2))
        val output = innerDirs.distinct.toList
        return output.sorted
    } catch{
        case e: Exception => return List[String]()
    }
}

// This method allows us to start a new simulation from the last snapshot stored in "dir" folder

def getLastSnapshots(dir:String): Tuple2[DataFrame, DataFrame]={
    try{
        var files = getFiles(dir, true)
        var lastV = load(dir, files.filter(_ contains("Vsnap")).last)
        var lastE = load(dir, files.filter(_ contains("Esnap")).last)
        return (lastV, lastE)
    }
    catch{
        case e: Exception => return null
    }
}

def deleteAll(dir:String): Unit={
    val directory = new Directory(new File("bacteriaFiles/"+dir))
    directory.deleteRecursively()
}

/* Following properties rule dinamics of each node

buffer: array che viene riempito in run in base al n. di migranti
// edgeInfo: è un dizionario che ha per chiave il nodo destinatario e valore numerosità popolazione in uscita
// nEdges: è il numero di outer links del vertice corrente

*/
case class Properties(var distribution: Array[Long] = Array(), 
                      var buffer: Array[Long] = Array(),
                      var edgeInfo: Map[Long, Long] = Map().withDefaultValue(0L),
                      var nEdges:Int = 0,
                      var point:Point = Point(-1, -1)
                     ){
    
// serialize() method "wraps" current properties in a string. This is mandatory for message passing
    
    if(buffer.size==0) buffer = Array.fill(distribution.size)(0L)
    def serialize(): String = {var output = distribution.mkString(",")+";"+
                                            buffer.mkString(",")+";"+
                                            edgeInfo.mkString(",").replace(" -> ", ":")+";"+
                                            nEdges+";"+
                                            point.value_X+";"+
                                            point.value_Y;
                                //println(output+"-->"+point)
                                return output;
                               }
    
}

// deserialize() rricostruise l'istanza della classe properties a partire dalla stringa output di serialize

def deserialize(serial:String): Properties = {
        try{
        val groups = serial.split(";");
        val distribution = groups(0).split(",").map(x=>x.toLong);
        val buffer = groups(1).split(",").map(x=>x.toLong);
        val groupMap = groups(2).split(",")
        var edgeInfo:Map[Long, Long] = Map().withDefaultValue(0L)
        if(groupMap(0)!=""){
            groupMap.map(x => edgeInfo += (x.split(":")(0).toLong->x.split(":")(1).toLong));
        }
        val nEdges = groups(3).toInt;
        val point = Point(groups(4).toInt, groups(5).toInt);
        return  Properties(distribution, buffer, edgeInfo, nEdges, point)
        }
    catch{
        case e: Exception => {println("Exception: "+serial); return null}
    }
}

case class EdgeProperties(var distance:Double = 0.0, 
                          var traffic:Long = 0L){
    
    def serialize(): String =   distance+";"+
                                traffic
}

def deserializeEdge(serial:String): EdgeProperties = {
    try{
        val groups = serial.split(";");
        val distance = groups(0).toDouble;
        val traffic = groups(1).toLong;
        return  EdgeProperties(distance, traffic)
        }
    catch{
        case e: Exception => {println("Exception: "+serial); return null}
    }
}

// there methods rules graph dynamics by upgrading old graph with a new ones.

// countEdges is the first example of aggregatemessage method usage.
// an Int message = 1 is sent to each source vertex for each edge (for each triplet in specific)
// doing this each vertex know how many edges exit from it (also the self looping)
// each vertex upgrades it's nEdges attribute in property and a new "aware" graph is producted.


def countEdges(g:Graph[Properties,String]):Graph[Properties,String] = {
    val verts = g.aggregateMessages[Int](_.sendToSrc(1), _+_).sortBy(_._1).join(g.vertices);
    verts.foreach(x => x._2._2.nEdges = x._2._1)
    var graph = Graph[Properties,String](verts.map(x => (x._1, x._2._2)), g.edges);
    return graph
    }

def startEdgeCount(bacteriaGraph:Graph[Properties,String]): Graph[Properties,String]={
    println("Edges's number algorithm started!")
    val now = Calendar.getInstance().getTime();
    var initialized_bacteriaGraph = countEdges(bacteriaGraph)
    initialized_bacteriaGraph.vertices.foreach(x =>println("Id: "+x._1+"  ->  "+ x._2.nEdges))
    val ex_time = (Calendar.getInstance().getTime().getTime() - now.getTime())/1000
    println("Execution time (s): "+ex_time+"\n");
    return initialized_bacteriaGraph
}

def distance(x1: Int, y1: Int, x2: Int, y2: Int) = {
    math.sqrt(math.pow(x1 - x2, 2) + math.pow(y1 - y2, 2))
}

// This methods check for minimum distance among nodes. It depends from the "minVerticesDistance" variable
// assigned in front-end. Default value for minimum distance is 10

def checkMinDist(point:Point, verts:ListBuffer[(VertexId, Properties)], minDist:Int): Boolean={
    var dist = verts.map(x => x._2.point.dist(point))
    if(dist.size>0) return dist.toList.min>minDist
    else return true
}

def serializeArray(array:Array[Long]): String={
    var output = ""
    array.foreach(x=> output+=","+x.toString)
    return output.substring(1)
}

def deserializeArray(array:String): Array[Long]={
    return array.split(",").map(x=> x.toLong).toArray
}

def toLong(x: Any): Long = x match{
    case i:Long => i.toLong
    case i:String => i.toLong
    case _ => 0L
}
def toInt(x: Any): Int = x match{
    case i:Int => i.toInt
    case i:String => i.toInt
    case _ => 0
}
def toSstring(x: Any): String = x match{
    case i:String => i
    case _ => ""
}

def getVertex(x:org.apache.spark.sql.Row): (VertexId, Properties)={
    var distribution = deserializeArray(toSstring(x(3)));
    var point = Point(toInt(x(4)), toInt(x(5)));
    return ((toLong(x(0)), Properties(distribution, point = point)))
}

def initializeGraphFromSnapshot(set:settings): Graph[Properties,String]={
    println("Loading graph");
    import scala.util.control.Breaks._
    import sqlContext.implicits._
    
    var snaps = getLastSnapshots(set.expName)
    var Vsnap = snaps._1
    var Esnap = snaps._2
    var vertexRDD: RDD[(VertexId, Properties)] = Vsnap.map(x=> getVertex(x)).rdd
    val edgeRDD:RDD[Edge[String]] = Esnap.rdd.map(row => 
                Edge(toLong(row(0)), toLong(row(1)), toSstring(row(2))+";"+toLong(row(3))))
    //vertexRDD.collect().foreach(println)
    //edgeRDD.collect().foreach(println)
    //var finalGraph = startEdgeCount(Graph(vertexRDD, edgeRDD))
    return Graph(vertexRDD, edgeRDD)
}

/*
This method initialize the graph with user defined settings.
If nodes location points are missing, the latter are asssigned by random choice in a "dim_space" side square

Crosscheck of "i" and "tent" variables allow us to identify badly placed nodes.
In such cases location point is riassigned.
    
    
-breakable- behaviour:
The instruction {break;} works only inside breakable{...} method. 
It allows us to escape all inside breakable{...} when the corresponding if-statement is True
*/

def initializeGraph(set:settings): Graph[Properties,String]={
    println("Initializing graph")
    import scala.util.control.Breaks._
    import sqlContext.implicits._
    //val r_in = scala.util.Random
    //r_in.setSeed(1002L)
    var vertices = new ListBuffer[(VertexId, Properties)]();
    // tolto Long finale per la popolazione            , Long
    var toCreate_Edges = new ListBuffer[(Long, Int, Int)]();
    var i = 0;
    var tent = 0;
    val maxTent = 100*set.pop.size+1000;
    val strainNum = set.pop.map(x => x.distribution.size).max;
    while(vertices.size != set.pop.size ) {
        tent += 1
        if(tent == maxTent) throw new Exception("Can't initialize the positions given the constraints")
        var casual_point = false
        var point = Point(set.pop(i).location.value_X, set.pop(i).location.value_Y)
        if(set.pop(i).location.value_X <0 || set.pop(i).location.value_X >set.dim_space){
            point = Point(set.r_in.nextInt(set.dim_space), set.r_in.nextInt(set.dim_space)); 
            casual_point=true;
        }
        breakable{
            if(!checkMinDist(point, vertices, set.minVerticesDistance) && casual_point) {break;}
            
            var distribution = Array.fill(strainNum)(0L)
            
            if(set.pop(i).distribution.size == 0){
                val index = i%strainNum
                if(set.pop(i).popIn>=0) distribution(index) = set.pop(i).popIn;
                else distribution(index) = set.r_in.nextInt(set.range_popSize(1) - set.range_popSize(0)) + 
                                            set.range_popSize(0)
            }
            else{
                for(j<-0 until distribution.size) distribution(j) = set.pop(i).distribution(j) 
            }
            
            vertices += ((i, Properties(distribution, point=point)))
            // ho tolto distribution.sum che sembra inutile    , distribution.sum
            toCreate_Edges += ((i, point.value_X, point.value_Y))
            
            i+=1
        }
    }
    var now = Calendar.getInstance().getTime();
    var vertexPos = sc.parallelize(toCreate_Edges)
    var edgesRDD = vertexPos.cartesian(vertexPos) // Cross product
                            .map(x => (x._1._1.toLong, x._2._1.toLong, // compute distance
                                       distance(x._1._2, x._1._3, x._2._2, x._2._3)))
                            .filter(_._3 < set.maxEdgeDistance) // filter too long edges
                            // Create Edges, must be String type
                            .map(x => Edge(x._1, x._2, EdgeProperties(x._3, 0L).serialize())) 
                            
    val vertexRDD: RDD[(VertexId, Properties)] = sc.parallelize(vertices)
    
    println("Execution time (s): "+(Calendar.getInstance().getTime().getTime() - now.getTime())/1000);
    return Graph(vertexRDD, edgesRDD)
}

var glob_set = settings()

// Nasty
def addToMap(key: Long, value: Long, map1: Map[Long, Long]) : Map[Long, Long] ={
    var map = map1.withDefaultValue(0L)
    map += (key -> (map1(key) + value))
    return map
}

// Really Nasty
def addMapToMap(map1: Map[Long, Long], map2: Map[Long, Long]) : Map[Long, Long] ={
    var map = map2.withDefaultValue(0L)
    var ot_map = map1.withDefaultValue(0L)
    if(map2.keys.size<map1.keys.size){ map = map1.withDefaultValue(0L); ot_map = map2.withDefaultValue(0L);}
    val keys = ot_map.keys.toList
    keys.foreach(key => map += (key -> (map(key) + ot_map(key))) )
    return map
}

// computeP() method takes the asolute distribution array as input and returns a relative ones

def computeP(act_distr: Array[Long]): Array[Double] = {
    val N = act_distr.sum
    return act_distr.map(x => x.toDouble/N)
}


// binomial() method takes a reproduction rate (or death rate) and the absolute population of a strain
// and returns the expected offsping or deaths.

def binomial(p:Double, n:Long): Long = {
    val threshold = 2000;
    if(n<=0) return 0L;
    var result = 0L;
    if(n>threshold) result = (n*p + glob_set.r_ev.nextGaussian()*math.sqrt(n*p*(1-p))).toLong
    else result = Array.fill(n.toInt)(glob_set.r_ev.nextFloat).map(x => if(x<p) 1L else 0L).reduce(_+_)
    if(result<0) result=0L
    if(result>n) result=n
    return result
}

/*
def binomial(p:Double, n:Long): Long = {
    val threshold = 2000;
    if(n<=0) return 0L;
    var result = 0L;
    if(n>threshold) result = (n*p + glob_set.r_ev.nextGaussian()*math.sqrt(n*p*(1-p))).toLong
    else result = Array.fill(n.toInt)(glob_set.r_ev.nextFloat).map(x => if(x<p) 1L else 0L).reduce(_+_)
    if(result<0) result=0L
    if(result>n) result=n
    return result
}
*/

// Sigmoidal function. Output 0<=1 

def Sigma(x:Double): Double={
    val a = 7
    val b = 0.5
    val s0 = 1/(1+math.exp(a*b))
    val sigma_noS = 1/(1+math.exp(-a*(x-b))) - s0
    var sigma = sigma_noS * 1/(1/(1+math.exp(-a*(1-b))) - s0)
    if(sigma>1) sigma=1
    if(sigma<0) sigma=0
    return sigma
}


// ricomputa la probabilità p (di riproduzione o di morte) in base al flag di interazione
// il flag è calcolato nel metodo replication

def recomputeForInteraction(p:Double, flag:Double): Double = {
    val sum = p + flag
    if (sum<0) return 0.0
    else if (sum>1) return 1.0
    else return sum
}


def newDistribution(distribution: Array[Long],newrates: ArrayBuffer[Tuple2[Double,Double]]):Array[Long] = {
    var newPop = distribution
    for (strain <- 0 until distribution.size){
        val add = binomial(newrates(strain)._2,distribution(strain))
        val sub = binomial(newrates(strain)._1,distribution(strain))
        newPop(strain) += add 
        newPop(strain) -= sub
        if (newPop(strain)<0) newPop(strain) = 0
    }
    return newPop
}


def selectFromDistribution(distribution:Array[Long], P:Array[Double], Nedges:Long): Array[Long]={
    var T_double = Array.fill(distribution.size)(0.0)
    for(i<-0 until distribution.size){
        T_double(i) = binomial(P(i), distribution(i))
    }
    T_double = T_double.map(x=> x/T_double.sum*Nedges)
    for(i<-0 until distribution.size){
        if(T_double(i)>distribution(i)) T_double(i) = distribution(i)
    }
    return T_double.map(x=> x.toLong)
}


def selectToMigrate(distribution: Array[Long], nEdges: Int, edgeAttr:String): Array[Long] = {
    val distance = deserializeEdge(edgeAttr).distance
    val distFact = math.pow(2.72, math.log(0.5)*distance/glob_set.halving);
    val N = distribution.sum
    //val noise = glob_set.r_ev.nextFloat * 0.2 + 0.8
    val Ntransf = binomial(glob_set.transfP, (N * distFact * glob_set.maxTransf).toLong)
    val Tmax = distribution.map(x => x/nEdges)
    var P = computeP(distribution)
    return selectFromDistribution(Tmax, P, Ntransf/nEdges)
}

// questa funzione restituisce un numero tra 0 e 0.51 che è il balancefactor 
// tanto più value (che è la popolazione del vertice si avvicina a limit (che è il ceiling)
// sotto ceiling/2 restituisce 0

def supLim(value:Long, limit:Long): Double={
    if(value < limit/2) return 0.0
    val a = 1/limit.toDouble*10
    val b = limit/2
    return 0.51/(1+math.pow(2.27, -a*(value-b)))
}

def replication(distribution: Array[Long]): Array[Long] = {
    var newPop = distribution
    var nTot = distribution.sum
    var newrates = ArrayBuffer[Tuple2[Double,Double]]()
    
    for(strain<-0 until newPop.size) {
        var dp = 0.0;
        var rp = 0.0;
        var nstrain = distribution(strain)
        // copia i parametri di riproduzione e morte senza interazione per lo strain corrente
        if(glob_set.strains.size > strain) dp = glob_set.strains(strain).deathP
        if(glob_set.strains.size > strain) rp = glob_set.strains(strain).reprP
        
        var dInter= 0.0
        var rInter = 0.0
        var fraction = 0.0
        val balanceFact = supLim(newPop(strain), glob_set.ceiling)
        // ciclo su tutti i possibili strain meno quello fissato dal ciclo for esterno
        for(otherStrain<-0 until newPop.size) if(otherStrain!= strain){
            
            // fraction è un valore compreso tra 0 e 1 e dipende dalla numerosità reciproca
            // degli elementi interagenti. 
            
            // Se ad esempio strain A interagisce con strain B, la fraction è massima quando
            // la numerosità di A è pari alla numerosità di B.
            
            // fraction non dipende dalla numerosità totale (hp di perfetto miscelamento)
            
            var fraction = newPop(otherStrain).toDouble/nTot.toDouble
            
            // questi IF verificano se all'interno del dizionario di interazione è presente "otherStrain"
            // ciò indica che strain interagisce con otherStrain.
            // In questo caso dFlag e rFlag accumulano i valori di interazione del dizionario 
            // moltiplicati per la "fraction" che rende la forza dell'interazione dipendente dalla 
            // popolazione reciproca.
            if(glob_set.strains(strain).deadlyInter.contains(otherStrain))
                dInter += glob_set.strains(strain).deadlyInter(otherStrain) * fraction
            if(glob_set.strains(strain).reproductiveInter.contains(otherStrain))
                rInter += glob_set.strains(strain).reproductiveInter(otherStrain) * fraction
        }
        
        // ciclati tutti gli "otherStrain" viene effettuato l'upgrade definitivo dei tassi di morte e
        // riproduzione dello strain corrente (ciclo esterno)
        
        dp += balanceFact
        rp -= balanceFact
        //println("pre-update strain no "+strain.toString+" dp: "+dp.toString+" rp: "+rp.toString)
    
        dp = recomputeForInteraction(dp, dInter)
        rp = recomputeForInteraction(rp, rInter)
        //println("UPDATED strain no "+strain.toString+" dp: "+dp.toString+" rp: "+rp.toString)
        //println("\n")
        
        // questa variabile raccoglie in un array di 2Tuple le nuove e definitive percentuali 
        // di riproduzione e morte per tutti gli strain.

        newrates += Tuple2(dp,rp)



    }
    
    // in base alle nuovi rates di morte e riproduzione viene effettuato l'upgrade della popolazione 
    // del vertice.

    newPop = newDistribution(distribution,newrates)

    return newPop
}


def show(x: Option[String]) = x match {
    case Some(i) => deserializeEdge(i)
    case None => EdgeProperties()
}

// 
def sendMsg(ec: EdgeContext[String,String,String]): Unit = {
    val src = deserialize(ec.srcAttr)
    val dst = deserialize(ec.dstAttr)
    if(ec.dstId == ec.srcId) {
        val toSendDst = Properties(dst.distribution, edgeInfo=dst.edgeInfo, nEdges=src.nEdges, point=dst.point)
        ec.sendToDst(toSendDst.serialize());
        return
    }
    
    var selected = selectToMigrate(src.distribution, src.nEdges, ec.attr)
    for (i <- 0 until src.distribution.size) 
    {
        src.buffer(i) = -selected(i);
        dst.buffer(i) = selected(i);
    }
    
    var map = Map(ec.dstId -> selected.sum);
    // perché hai mandato tutta questa roba ? Non bastava solo il buffer?
    val toSendDst = Properties(dst.distribution, dst.buffer, map, dst.nEdges, point=dst.point)
    val toSendSrc = Properties(src.distribution, src.buffer, map, src.nEdges, point=src.point)
    ec.sendToDst(toSendDst.serialize());
    ec.sendToSrc(toSendSrc.serialize());
    
}

def mergeMsg(accum: String, msg: String): String = {
    var acc_prop = deserialize(accum)
    var property = deserialize(msg)
    
    property.edgeInfo = addMapToMap(property.edgeInfo, acc_prop.edgeInfo)
    for (i <- 0 until property.buffer.size) property.buffer(i) += acc_prop.buffer(i)
    
    return property.serialize()
}

def compute(msg: String): String = {
    var property = deserialize(msg)
    for (i <- 0 until property.distribution.size) 
    {
        property.distribution(i) += property.buffer(i);
    }
    property.distribution = replication(property.distribution)
    return property.serialize()
}

def computeTraffic(msg: String, id:VertexId): List[Edge[String]] = {
    var property = deserialize(msg)
    var a = ListBuffer[Edge[String]]()
    for(i<-0 until property.edgeInfo.keys.toList.size){
        a += (Edge(property.edgeInfo.keys.toList(i), id, 
                EdgeProperties(0.0, property.edgeInfo(property.edgeInfo.keys.toList(i))).serialize()))
    }
    return a.toList
}

def truncate(array:Array[Long], n:Int, separator:String): String ={
    var trList = ListBuffer[Long]()
    for(i<-0 until math.min(n, array.size)) trList += array(i)
    return trList.toList.mkString(separator)
}
def truncateMap(map:Map[Long, Long], n:Int, separator:String): String ={
    var trList = ListBuffer[String]()
    for(i<-0 until math.min(n, map.keys.toList.size)) 
        trList += map(map.keys.toList(i))+"->"+map.keys.toList(i)
    return trList.toList.mkString(separator)
}

def read(id: VertexId, msg: String): Unit = {
    val property = deserialize(msg)
    val n=20
    println("Id: " + id + "    nLink="+ property.nEdges+"  ->  pop size: "+property.distribution.sum)
    println("    Distrib: " + truncate(property.distribution, n, " ") +"  ")
    println("    Edge Info: " + truncateMap(property.edgeInfo, n, " ") +"  ")
    var incr_str = "+" + truncate(property.buffer, n, " +");
    println("    Increment: " + incr_str.replace("+-", "-").replace("+0", "0") +"  ")
}

def concatenate(rdd: RDD[List[Edge[String]]]): RDD[Edge[String]] ={
    var newEdges = ListBuffer[Edge[String]]()
    rdd.collect.foreach(x => x.foreach(y => newEdges += y))
    return sc.parallelize(newEdges)
}

def mergeEdgeProperties(prop1: EdgeProperties, prop2: EdgeProperties): EdgeProperties={
    return EdgeProperties(math.max(prop1.distance, prop2.distance), prop1.traffic + prop2.traffic)
}

case class Stats(var distribution: Array[Long] = Array()){
    val sum = distribution.sum
    val G = 1 - distribution.map(x => math.pow((x.toDouble/sum),2)).sum
    val G_norm = G/(distribution.size-1)*distribution.size
}


def getSnapShot(vertices:Array[(VertexId, String)], edges:Array[Edge[String]]): Unit ={
    import sqlContext.implicits._
    var VertexList = ListBuffer[(Long, Long, Double, String, Int, Int)]()
    var EdgeList = ListBuffer[(Long, Long, Double, Long)]()
    
    vertices.foreach(x => VertexList += ((x._1, 
                                          Stats(deserialize(x._2).distribution).sum,
                                          Stats(deserialize(x._2).distribution).G_norm,
                                          serializeArray(deserialize(x._2).distribution),
                                          deserialize(x._2).point.value_X,
                                          deserialize(x._2).point.value_Y,
                                         )))
    edges.foreach(x => EdgeList += ((x.srcId, x.dstId, deserializeEdge(x.attr).distance,
                                     deserializeEdge(x.attr).traffic)))
    
    
    val vDF = sc.parallelize(VertexList).toDF()
                .withColumnRenamed("_1", "id")
                .withColumnRenamed("_2", "pop")
                .withColumnRenamed("_3", "heterogeneity")
                .withColumnRenamed("_4", "distribution")
                .withColumnRenamed("_5", "x")
                .withColumnRenamed("_6", "y")
    //vDF.createOrReplaceTempView("Vsnap"+iterations)
    save(vDF, glob_set.expName, "Vsnap"+iterations)
    val eDF = sc.parallelize(EdgeList).toDF()
                .withColumnRenamed("_1", "id1")
                .withColumnRenamed("_2", "id2")
                .withColumnRenamed("_3", "distance")
                .withColumnRenamed("_4", "traffic")
    save(eDF, glob_set.expName, "Esnap"+iterations)
    //eDF.createOrReplaceTempView("Esnap"+iterations)
                //.write.mode("overwrite").saveAsTable("Esnap"+iterations)
}

def clearDataframes(iteration:Int): Boolean={
    try{
        sqlContext.sql("drop table Vsnap"+iterations)
        sqlContext.sql("drop table Esnap"+iterations)
    }
    catch{
        case e: Exception => return false
    }
    return true
}


var iterations = 0;
var maxIterations = 0;
val printAt = Array[Long]() //10, 20, 30, 40, 50, 60, glob_set.nIterations)

def simulation(g:Graph[String,String]):Tuple2[Graph[String,String], Boolean] = {
    val now = Calendar.getInstance().getTime();
    if(iterations >= maxIterations) return (g, false);
    
    print("Iteration %d started  ".format(iterations+1));
    val newVertices = g.aggregateMessages[String](sendMsg, mergeMsg)
                    .map(x => Tuple2(x._1, compute(x._2)));
    var newEdges_map = concatenate(newVertices.map(y => computeTraffic(y._2, y._1)))
                        .map(x=> ((x.srcId, x.dstId), x.attr))
    val prevEdges_map = g.edges.map(x => ((x.srcId, x.dstId), x.attr))
    
    var newEdges = prevEdges_map.leftOuterJoin(newEdges_map)
            .map(x => (x._1, mergeEdgeProperties(deserializeEdge(x._2._1), show(x._2._2)).serialize()))
            .map(x => Edge(x._1._1, x._1._2, x._2))
    
    if(printAt contains iterations+1) {println(); newVertices.collect.foreach(x => read(x._1, x._2))}
    
    
    getSnapShot(g.vertices.collect(), newEdges.collect())
    var newGraph = Graph(newVertices, g.edges);
    newVertices.unpersist()
    newEdges_map.unpersist()
    prevEdges_map.unpersist()
    newEdges.unpersist()
    g.unpersist()
    val ex_time = (Calendar.getInstance().getTime().getTime() - now.getTime()).toFloat/1000
    println(".. and completed in "+ex_time+"s");
    if(ex_time>10) return (newGraph, false)
    if(printAt contains iterations+1) println("    -----")
    
    iterations+=1;
    return (newGraph, true)
    }



def startSimulation(bacteriaGraph:Graph[Properties,String], set:settings): Unit={
    glob_set = set;
    var alreadyDone = getFiles(glob_set.expName).filter(_ contains "Vsnap").size
    if(alreadyDone!=0) iterations = alreadyDone
    else iterations = 0;
    maxIterations = glob_set.nIterations + iterations
    println("|----------------------|\n    Algorithm started!")
    val firstTime = Calendar.getInstance().getTime();
    
    var dynamicGraph = startEdgeCount(bacteriaGraph).mapVertices((x:VertexId, y:Properties) => y.serialize())

    var continue = true;
    while(continue){
        var tuple = simulation(dynamicGraph);
        dynamicGraph = tuple._1
        continue = tuple._2
    }
    
    continue=true
    while(continue){continue = clearDataframes(iterations); iterations+=1}
    
    val ex_time = (Calendar.getInstance().getTime().getTime() - firstTime.getTime()).toFloat/1000
    println("|----------------------|\n")
    println("Total iterations: %d".format(iterations-1));
    println("Execution time: "+ex_time+"s\n");
    println("         -  per iteration: "+ex_time.toFloat/iterations+"s\n");
}



%%python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.animation as animation
import random
import math
import networkx as nx

def deserializeStrainDistr(distribution):
    if type(distribution) == str:
        strList = distribution.split(',')
        output = []
        for el in strList: output.append(int(el))
        return output
    else: return [distribution, ]

class getData:
    def __init__(self, directory):
        
        self.vertexSnaps, self.edgeSnaps = getDataframes(directory)
        totalIterations = len(self.vertexSnaps)
        self.vertexPosition = {}
        self.edgePosition = {}
        self.vertexDynamic = {}
        self.vertexDistribution = {}
        self.cumEdgeDynamic = {}
        self.edgeDynamic = {}

        for i in range(totalIterations):
            for row in self.vertexSnaps[i].iterrows():
                vertex = int(row[1][0])
                popSize = row[1][1]
                etherogen = float(row[1][2])
                if math.isnan(etherogen): etherogen = 0.0
                distribution = row[1][3]
                point = (row[1][4],row[1][5])
                if vertex not in self.vertexPosition:
                    self.vertexPosition[vertex] = point
                if self.vertexDynamic.get(vertex, None) == None: self.vertexDynamic[vertex] = []
                if self.vertexDistribution.get(vertex, None) == None: self.vertexDistribution[vertex] = []
                while(len(self.vertexDynamic[vertex]) < i): 
                    self.vertexDynamic[vertex].append((0, 0))
                    self.vertexDistribution[vertex].append([])
                    
                self.vertexDynamic[vertex].append((popSize, etherogen))
                self.vertexDistribution[vertex].append(deserializeStrainDistr(distribution))
                
                
        for i in range(totalIterations):
            for vertex in self.vertexDynamic:
                while(len(self.vertexDynamic[vertex]) <= i): 
                    self.vertexDynamic[vertex].append((0, 0))
                    self.vertexDistribution[vertex].append([])
                    

        for i in range(totalIterations):
            for row in self.edgeSnaps[i].iterrows():
                edge = (int(row[1][0]),int(row[1][1]))
                value = row[1][3]
                if edge not in self.edgePosition:
                    self.edgePosition[edge] = (self.vertexPosition[edge[0]], self.vertexPosition[edge[1]])

                if self.cumEdgeDynamic.get(edge, None) == None: 
                    self.edgeDynamic[edge] = []
                    self.cumEdgeDynamic[edge] = []
                while(len(self.cumEdgeDynamic[edge]) < i): 
                    self.cumEdgeDynamic[edge].append(0)
                self.cumEdgeDynamic[edge].append(value)
        
        for i in range(totalIterations):
            for edge in self.cumEdgeDynamic:
                while(len(self.cumEdgeDynamic[edge]) <= i): self.cumEdgeDynamic[edge].append(0)
                    
        
        for i in range(totalIterations):
            for edge in self.cumEdgeDynamic:
                if i>1: 
                    self.edgeDynamic[edge].append(self.cumEdgeDynamic[edge][i] - 
                                                      self.cumEdgeDynamic[edge][i-1])
                    if(self.edgeDynamic[edge][-1]<0): self.edgeDynamic[edge][-1]=0
                else: self.edgeDynamic[edge].append(self.cumEdgeDynamic[edge][0])
        
            
        print('Data extracted')
            
    def check(self):
        error = False
        self.nodeNumber = len(self.vertexPosition)
        if(len(self.vertexDynamic) != self.nodeNumber): error = True
        if(len(self.vertexDistribution) != self.nodeNumber): error = True
        self.edgeNumber = len(self.edgePosition)
        if(len(self.edgeDynamic) != self.edgeNumber): error = True
            
        self.iterNumber = len(self.vertexDynamic[list(self.vertexDynamic.keys())[0]])
        for k in self.vertexDynamic: 
            if len(self.vertexDynamic[k]) != self.iterNumber: error = True
        for k in self.vertexDistribution: 
            if len(self.vertexDistribution[k]) != self.iterNumber: error = True
        for k in self.edgeDynamic: 
            if len(self.edgeDynamic[k]) != self.iterNumber: error = True
        
        self.strainNumber = 0
        for k in self.vertexDistribution:
            for i in range(len(self.vertexDistribution[k])):
                lst = self.vertexDistribution[k][i]
                if len(lst) > self.strainNumber: self.strainNumber = len(lst)
        
        for k in self.vertexDistribution:
            for i in range(len(self.vertexDistribution[k])):
                lst = self.vertexDistribution[k][i]
                if len(lst) == 0: lst = [0]*self.strainNumber
                if len(lst) != self.strainNumber: error = True
            
        self.edgeNumber -= self.nodeNumber
        self.edgeNumber = int(self.edgeNumber/2)
        
        print('Errors: '+str(error))
        print('Number of nodes: '+ str(self.nodeNumber))
        print('Number of edges: '+ str(self.edgeNumber))
        print('Total strains: '+ str(self.strainNumber))
        print('Total iterations: '+ str(self.iterNumber))

def brightCol(color):
    col = mcolors.to_rgba(color)
    darkness = 0
    if col[0]<0.2: darkness += 1
    if col[1]<0.2: darkness += 1
    if col[2]<0.2: darkness += 1
    if col[1]!= 0 and col[1] != 0 and col[0]/col[1]>4 and col[1]/col[2]<2: return False
    if darkness < 2: return True
    else: return False
        
def drawGraph(data, nStrains=-1):
    if nStrains < 0: nStrains = data.strainNumber
    if nStrains > data.strainNumber: nStrains = data.strainNumber
    backcolor = 'black'
    labelFacerColor = (1,1,1)
    colorList = [(1,1,1)]
    fig, ax = plt.subplots()
    ax.margins(x=0, y=0)
    fig.set_figheight(11)
    fig.set_figwidth(30)
    rect = fig.patch
    rect.set_facecolor(backcolor)
    G = nx.Graph()
    G.add_nodes_from(data.vertexPosition)
    baseColors = [k for k in mcolors.BASE_COLORS if (k!='r' and k!='black' and brightCol(k))] + \
                [k for k in mcolors.CSS4_COLORS if (k!='red' and k!='black' and brightCol(k))]
    
    colors = [random.choice(baseColors) for i in range(len(data.vertexDynamic))]
    strainNames = {i:'S'+str(i) for i in range(data.strainNumber)}
    strainNames[0]='A'
    strainNames[1]='B'
    strainNames[2]='C'
    strainNames[3]='D'
    strainNames[4]='E'
    strainNames[5]='F'
    strainNames[6]='G'
    strainNames[7]='H'
    strainNames[8]='I'
    strainNames[9]='J'
    strainNames[10]='K'

    for id_ in data.vertexPosition:
        list_ = {id_:data.vertexPosition[id_]}

    def animate(i):
        ax.clear()
        rect = ax.patch
        title = 'iteration: '+str(i)+'/'+str(data.iterNumber)
        title += '\n'+'Number of nodes: '+ str(data.nodeNumber)
        title += '\n'+'Number of edges: '+ str(data.edgeNumber)
        title += '\n'+'Total strains: '+ str(data.strainNumber)
        ax.set_facecolor(backcolor)
        ax.text(-40, 120, s=title, style='italic', size=20, bbox={'facecolor':labelFacerColor,
                                                                'alpha':0.5, 'pad':15})
        
        try:
            nx.draw_networkx_nodes(G, pos={'A':(-20,-20)}, nodelist={'A':0}, node_color='black')
            nx.draw_networkx_nodes(G, pos={'B':(120,-20)}, nodelist={'B':0}, node_color='black')
            nx.draw_networkx_nodes(G, pos={'C':(120,120)}, nodelist={'C':0}, node_color='black')
            nx.draw_networkx_nodes(G, pos={'D':(-20,120)}, nodelist={'D':0}, node_color='black')
            
            #print(data.vertexDynamic)
            history = {id_:sum(data.vertexDynamic[id_][:i][0]) for id_ in data.vertexDynamic}
            tempPop = {id_:data.vertexDynamic[id_][i] for id_ in data.vertexDynamic}
            maxPop = sum([tempPop[k][0] for k in tempPop])
            tempDistr = {id_:data.vertexDistribution[id_][i] for id_ in data.vertexDistribution}
            
            tempTraffic = {id_:data.edgeDynamic[id_][i] for id_ in data.edgeDynamic}
            maxTraf = max([tempTraffic[k] for k in tempTraffic])

            for j in data.vertexPosition:
                list_ = {j:tempPop[j]}
                size = math.sqrt(tempPop[j][0]/maxPop)*40000 + 1
                
                if tempPop[j][0]<1: 
                    size = math.sqrt(1/maxPop)*20000 + 1
                    alpha = 0.4
                    if(history[j]==0): color = 'white'
                    else: color = 'red'
                else: 
                    alpha = 1 - tempPop[j][1]/2
                    color = colors[j]
                    
                nx.draw_networkx_nodes(G, pos=data.vertexPosition, nodelist=list_, 
                                       node_color=color, alpha=alpha, node_size=size)

            for j in data.edgePosition:
                list_ = {j:data.edgePosition[j]}
                if(maxTraf == 0):
                    width = 0
                    alpha = 0
                else: 
                    #maxwidth = int(math.pow(float(max([tempTraffic[k] for k in tempTraffic]))
                    #                                /maxTraf, 4))+2
                    width = int(math.pow(float(tempTraffic[j])/maxTraf, 3)*300+ 1)
                    #width = tempTraffic[j]/400
                    
                    alpha = 0.9*width/maxTraf + (1-0.9)
                    if alpha <0: alpha = 0
                    elif alpha >1: alpha = 1
                    if width>60: width = 60
                nx.draw_networkx_edges(G, pos=data.vertexPosition, edgelist=list_, width=width, 
                                       edge_color='w', alpha=alpha)

            labels = {k:'ID:'+str(k)+'\n'+str(int(tempPop[k][0]))+'\n'+'{:.2f}'.format(tempPop[k][1]) \
                      for k in tempPop if tempPop[k][0]>1}
            for k in labels:
                orderedDistr = np.array(tempDistr[k])/np.array(tempDistr[k]).sum()
                strainValues = {}
                for i in range(len(orderedDistr)):
                    strainValues[i] = int(orderedDistr[i]*1000)/10
                ordered = sorted(strainValues.items(), key=lambda x: x[1], reverse=True)
                for item in ordered[:nStrains]:
                    labels[k] += '\n'+ strainNames[item[0]] +'~'+str(item[1])+'%'
                    
            nx.draw_networkx_labels(G, data.vertexPosition, labels, font_size=20, 
                                    font_color='r', font_weight='bold')
            
            
        except Exception as e: print(e, i)

    anim = animation.FuncAnimation(fig, animate, np.arange(0, data.iterNumber), interval=200, blit=False)
    plt.show()

print('Python functions defined')

import org.apache.spark.graphx.VertexId
import org.apache.spark.graphx


def DropNode(g:Graph[Properties,String], id:Long): Graph[Properties,String]={
    var newVerts = g.vertices.filter(x=> x._1 != id)
    var newEdges = g.edges.filter(x=> x.srcId != id).filter(x=> x.dstId != id)
    return Graph(newVerts, newEdges)
}

def AddNode(g:Graph[Properties,String], node:Node): Graph[Properties,String]={
    var newVerts = g.vertices.map(x => (x._1, x._2))
    var nDistr = newVerts.map(x => x._2.distribution.size).max
    var id = newVerts.map(x => x._1).max + 1
    var distribution = Array.fill(nDistr)(0L)
    var r = scala.util.Random
    if(node.distribution.size==0) distribution(r.nextInt(nDistr)) = node.popIn
    else for(i<-0 until nDistr) if(node.distribution.size>i) distribution(i) += node.distribution(i)
    var vertex = sc.parallelize(ListBuffer((id.toLong, 
                                Properties(
                                    distribution = distribution,
                                    point = node.location,
                                ))))
    newVerts = newVerts ++ vertex
    var newEdges = g.edges ++ sc.parallelize(ListBuffer(Edge(id, id, EdgeProperties().serialize())))
    return Graph(newVerts, newEdges)
}

def AddEdge(g:Graph[Properties,String], id1:Long, id2:Long): Graph[Properties,String]={
    var points = g.vertices.map(x => (x._1, x._2.point)).filter(x => x._1==id1 || x._1==id2)
    if (points.count!=2) {println("Can't create the edge ("+id1+","+id2+")");return g}
    var arr = points.map(x => (x._2.value_X, x._2.value_Y)).collect()
    var dist = Point(arr(0)._1, arr(0)._2).dist(Point(arr(1)._1, arr(1)._2))
    
    var newEdges = g.edges ++ sc.parallelize(ListBuffer(
                                    Edge(id1, id2, EdgeProperties(dist,0L).serialize()),
                                    Edge(id2, id1, EdgeProperties(dist,0L).serialize()),
                                ))
    return Graph(g.vertices, newEdges)
}

def DropEdge(g:Graph[Properties,String], id1:Long, id2:Long): Graph[Properties,String]={
    if(id1==id2) return g;
    var newEdges = g.edges.filter(x => (x.srcId!=id1 && x.dstId!=id2) || (x.srcId!=id2 && x.dstId!=id1))
    return Graph(g.vertices, newEdges)
}

def DropAllEdges(g:Graph[Properties,String], id:Long): Graph[Properties,String]={
    var newEdges = g.edges.filter(x => (x.srcId == x.dstId) || (x.srcId!=id && x.dstId!=id))
    return Graph(g.vertices, newEdges)
}

def AddFullConnection(g:Graph[Properties,String], id_l:Long = -1L): Graph[Properties,String]={
    try{
        var id = id_l
        if(id == -1) id = g.vertices.map(x => x._1).max
        var point = g.vertices.map(x => (x._1, x._2.point)).filter(x => x._1 == id)
        var points = g.vertices.map(x => (x._1, x._2.point)).filter(x => x._1 != id).cartesian(point)
        val map1 = points.map(x => (x._1._1, id, x._1._2.dist(x._2._2)))
        val map2 = map1.map(x => (x._2, x._1, x._3))
        var oldedges = (g.edges.map(x=> (x.srcId, x.dstId)) ++ g.edges.map(x=> (x.dstId, x.srcId))).collect()
        
        var map = map1++map2
        map = map.filter(x => !(oldedges contains (x._1, x._2)))
        
        var newEdges = g.edges ++ map.map(x => Edge(x._1, x._2, EdgeProperties(x._3, 0).serialize()))
        return Graph(g.vertices, newEdges)
    }
    catch{
        case e: Exception => {println("Exception in adding connections"); return g;}
    }
}


def initialize(sett:settings, fromLastSnapshot:Boolean): Graph[Properties,String]={
    sett.makeConsistent()
    if(!fromLastSnapshot) deleteAll(sett.expName)
    sett.lastSnapshots = getLastSnapshots(sett.expName)
    
    if(sett.lastSnapshots == null) return initializeGraph(sett)
    else return initializeGraphFromSnapshot(sett)
}

def Start(sett:settings, bacteriaGraph:Graph[Properties,String], iterations:Int = -1): Graph[Properties,String]={
    if(iterations!= -1) sett.nIterations = iterations
    
    var iterationCycle = Array.fill(sett.nIterations/50 + 1)(50)
    iterationCycle(iterationCycle.size-1) = sett.nIterations%50
    for(i<-0 until iterationCycle.size){
        sett.lastSnapshots = getLastSnapshots(sett.expName)
        sett.nIterations = iterationCycle(i)
        if(i>0) startSimulation(initializeGraphFromSnapshot(sett), sett)
        else startSimulation(bacteriaGraph, sett)
    }
    return initializeGraphFromSnapshot(sett)
}

class BGraph(settIn:settings, fromLastSnapshot:Boolean){
    var internalSettings:settings = settIn
    var graphx: Graph[Properties, String] = initialize(settIn, fromLastSnapshot)
    def start(iterations:Int = -1, sett:settings=null): BGraph= 
        { if(sett!=null) internalSettings = sett; graphx=Start(internalSettings, graphx, iterations); this }
    def addNode(node:Node): BGraph= { graphx=AddNode(graphx, node:Node); this }
    def dropNode(id:Long): BGraph= { graphx=DropNode(graphx, id); this }
    def addEdge(id1:Long, id2:Long): BGraph= { graphx=AddEdge(graphx, id1, id2); this }
    def dropEdge(id1:Long, id2:Long): BGraph= { graphx=DropEdge(graphx, id1, id2); this }
    def addFullConnection(id_l:Long = -1L): BGraph= { graphx=AddFullConnection(graphx, id_l); this }
    def dropAllEdges(id:Long): BGraph= { graphx=DropAllEdges(graphx, id); this }
}

/*
var array1 = ArrayBuffer[Node]()
var array2 = ArrayBuffer[Node]()
var array3 = ArrayBuffer[Node]()
var array4 = ArrayBuffer[Node]()
var array5 = ArrayBuffer[Node]()

for (i <- 0 until 10){
    array1 += Node(distribution=Array(10,5,5))
    array2 += Node(distribution=Array(100,50,50))
    array3 += Node(distribution=Array(1000,500,500))
    array4 += Node(distribution=Array(10000,5000,5000))
    array5 += Node(distribution=Array(100000,50000,50000))
}

var arrays = (array1++array2++array3++array4++array5).toArray
*/

var arr = ArrayBuffer[Node]()

for (i <- 0 until 10){
    arr += Node(distribution=Array(1000,1000,1000,1000,1000))
}

var arrays = arr.toArray

def ProgramDebug(): Unit={
    var mysettings = settings()
    mysettings.expName = "debug_inal"
    deleteAll(mysettings.expName)
    mysettings.strains = Array(
        Strain(deathP = 0.5, reprP = 0.5,reproductiveInter=Map(1L->0.0018)),
        Strain(deathP = 0.5, reprP = 0.5,reproductiveInter=Map(0L->0.0018)),
        Strain(deathP = 0.5, reprP = 0.5,deadlyInter=Map(3L->0.01)),
        Strain(deathP = 0.5, reprP = 0.5,deadlyInter=Map(2L->0.01)),
        Strain(deathP = 0.5, reprP = 0.5))
    
    mysettings.pop = arrays
    
    /*
    mysettings.strains = Array(
        Strain(deathP = 0.5, reprP = 0.5),
        Strain(deathP = 0.5, reprP = 0.5),
        Strain(deathP = 0.0, reprP = 0.5,deadlyInter=Map(3L->1e-4)),
        Strain(deathP = 0.5, reprP = 0.0,reproductiveInter=Map(2L->1e-4)))
                       
    mysettings.pop = Array(
        Node(distribution=Array(6000,3000,0,0)),
        Node(distribution=Array(0,0,500,500)))
    */
                   
    mysettings.range_popSize = Array(100*1000, 500*1000)
    mysettings.dim_space = 100
    mysettings.nIterations = 3000
    mysettings.maxEdgeDistance = 200
    mysettings.minVerticesDistance = 1
    mysettings.maxTransf = 0.5  // put double here!
    mysettings.transfP = 0.1 // put double here!
    mysettings.halving = 40
    mysettings.ceiling = 100000000
    mysettings.set_reallyCasual = false
    mysettings.in_reallyCasual = false
    mysettings.ev_reallyCasual = false
    
    
    var graph = new BGraph(mysettings, fromLastSnapshot=true)
                .start()
    
}



ProgramDebug();

// # %%python
// # data = getData("RWks3")
// # data.check()
// # attributes = [data.vertexPosition, data.edgePosition, 
// #               data.vertexDynamic, data.vertexDistribution, data.edgeDynamic]
// # drawGraph(data)


// var mysettings = settings()
// mysettings.strains = Array(
//     Strain(deathP = 0.5, reprP = 0.5),
//     Strain(deathP = 0.5, reprP = 0.5))

// mysettings.pop = Array(
//     Node(distribution = Array(1000,1000)),
//     Node(distribution = Array(1000,1000)))

//     mysettings.range_popSize = Array(100*1000, 500*1000)
//     mysettings.dim_space = 100
//     mysettings.nIterations = 3000
//     mysettings.maxEdgeDistance = 200
//     mysettings.minVerticesDistance = 1
//     mysettings.maxTransf = 0.9  // put double here!
//     mysettings.transfP = 0.5 // put double here!
//     mysettings.halving = 40
//     mysettings.ceiling = 1000000
//     mysettings.set_reallyCasual = false
//     mysettings.in_reallyCasual = false
//     mysettings.ev_reallyCasual = false

// var graph = new BGraph(mysettings, fromLastSnapshot=true)
