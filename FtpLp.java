import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;

import ilog.concert.*;
import ilog.cplex.*;

class edge {
	int src = 0;
	int dst = 0;
	int cap = 0;
	int delay = 0;
	int t = 0;
	int id;
	double[] flow= new double[30];
	public edge(int sr, int des, int cp, int ti, int idd, int del) {
		src = sr;
		dst = des;
		cap = cp;
		t = ti;
		id = idd;
		delay=del;
		for(int i=0; i<30; i++)
			flow[i]=0;
	}
}

class comm{
	int src=0;
	int des=0;
	int id;
	int req;
	int supsrc;
	int supdes;

	public comm(int sr, int de, int idd, int re, int susr, int sude){
		src = sr;
		des=de;
		id=idd;
		req=re;
		supsrc=susr;
		supdes=sude;
	}
}


class graph {
	int NumCopies = 0;
	int NumVert = 0;
	int NumEdges = 0;
	int NumOrigV;
	int NumOrigE;
	edge[][] tGraph = new edge[20000][20000]; 

	public graph(String filename, int T) throws IOException {
		NumCopies = T;
		String tmp;
		int a, b, c, d;
		String[] res = null;

		FileReader input1 = new FileReader(filename);
		BufferedReader bufRead1 = new BufferedReader(input1);
		tmp = bufRead1.readLine();
		int e = 0;
		res = tmp.split("\\s");
		a = Integer.parseInt(res[0]);
		NumOrigV=a;
		b = Integer.parseInt(res[1]);
		NumOrigE=b;
		if (tmp != null) {
			NumVert = a*NumCopies;
			NumEdges = b* NumCopies;
		}
		tmp = bufRead1.readLine();
		while (tmp != null) {
			if (tmp != null) {

				res = tmp.split("\\s");
				a = Integer.parseInt(res[0]);
				b = Integer.parseInt(res[1]);
				c = Integer.parseInt(res[2]);
				d=Integer.parseInt(res[3]);
				int delay = d;
				for (int i = 0; i < T; i++) {
					if(a!=b && delay+i < T){
						edge edges = new edge( a + (NumOrigV * i), b + (NumOrigV * (i+d)),
								c, i, e, d);
						tGraph[ a + (NumOrigV * i)][b + (NumOrigV * (i+d))] = edges;
					}
					e++;
				}
			}
			tmp = bufRead1.readLine();
		}
		bufRead1.close();
	}
}


class path {
	double flow = 0;
	Vector<edge> seq = new Vector<edge>();
	int source;
	int dest;

	public static void getpaths(comm[] comms, graph G, int nocom, double sol) throws IOException {

		path[][] p1=new path[30][300];
		int i=0;
		int l=0;
		int j=0;
		int k=0;
		double min=1000000;
		double min1=0;
		double flow=0;
		String filename1 = "paths.txt" ;
		FileWriter outFile = new FileWriter(filename1);
		PrintWriter out2 = new PrintWriter(outFile);
		System.out.println("finding paths");
		while(l<nocom)
		{
			//finding paths 
			while(flow< sol*comms[l].req)
			{
				i=comms[l].supsrc;
				p1[l][k]=new path();
				for(j=0;j<G.NumVert*2;j++)
				{   
					if(G.tGraph[i][j]!=null)
					{
						if(G.tGraph[i][j].flow[l] > 0.0)
						{
							p1[l][k].seq.add(G.tGraph[i][j]);
							if(min>G.tGraph[i][j].flow[l])
								min=G.tGraph[i][j].flow[l];
							i=j;
							if(i==comms[l].supdes)
							{
								min1=min;
								min=1000000;
								flow=flow+min1;
								break;
							}
						}
					}
				}
				//changing capacity of paths
				if(p1[l][k]!=null)
				{ 
					//System.out.println("flow=" + min1);
					System.out.println("path "+k +" for comm " + l +" with flow "+ min1 + "="); 
					out2.println("path "+k +" for comm " + l +" with flow "+ min1 + "=");
					for(edge e: p1[l][k].seq)
					{
						e.flow[l]=e.flow[l]-min1;
						System.out.print("("+ e.src+","+ e.dst + ") ->"); 
						if((e.src!=comms[l].supsrc)&& (e.dst!=comms[l].supdes))
						{
						int src=(e.src%G.NumOrigV);
						int dest=(e.dst%G.NumOrigV);
						int t=e.src/G.NumOrigV;
						int t1=e.dst/G.NumOrigV;
						out2.print("( "+ src+", t="+t+ ","+ dest +", t="+ t1 + ") ->");
					}
					   
					}
					k++;
					System.out.println();
					 out2.println();
				}
			}
			k=0;
			l++;
			flow=0;
		}
		System.out.println("done with finding paths");
		out2.close();
	}
}



public class FtpLp {

	public static void main(String[] args) throws IOException{

		String filename = "smallnet5.txt" ;

		String filename2="commdata.txt";
		int timeofgraph=5;
		double x=0 ;
		System.out.println(filename);
		graph G = new graph(filename, timeofgraph);
		System.out.println("enter number of commodities"); 
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		int nocom=Integer.parseInt(in.readLine());

		//writing commodity data to file 
		FileWriter outFile1 = new FileWriter(filename2);
		PrintWriter out1 = new PrintWriter(outFile1);
		out1.println("Comm source sink requirement supersource supersink");

		//Creating random commodity data
		comm[] comms= new comm[nocom]; 
		/*	for (int i=0;i<nocom; i++)
		{
			Random generator = new Random();
			int rand=generator.nextInt();
			int rand1=generator.nextInt();
			int rand2=generator.nextInt();
			rand=rand % (G.NumVert/(2*G.NumCopies));
			rand=Math.abs(rand);
			rand1=(Math.abs(rand1 %(G.NumVert/(2*G.NumCopies)))) + (G.NumVert/(2*G.NumCopies) -1); 
			rand2=(rand2 % 20) + 40;
			comm coms=new comm(rand, rand1, i, rand2, (G.NumVert+rand), (G.NumVert+rand1) );
			comms[i]=coms;
			out1.println(i + "\t" + rand + "\t" + rand1 + "\t"+ rand2+ "\t" + (G.NumVert+rand) + "\t" + (G.NumVert+rand1) );

		} */
	comms[0]= new comm(0, 4, 0, 55, (G.NumVert+0), (G.NumVert+4) );
	comms[1]= new comm(1, 3, 0, 55, (G.NumVert+1), (G.NumVert+3) );
		out1.close();



		//Creating super source and super sink
		int e=G.NumCopies*G.NumVert;
		for(int i=0;i<nocom;i++)
		{
			for(int j=0; j< G.NumCopies;j++)
			{
				edge edge1=new edge((comms[i].des + (G.NumOrigV * j) ), (comms[i].supdes),
						10000, G.NumCopies, e, 0);

				G.tGraph[(comms[i].des + (G.NumOrigV * j)  )][(comms[i].supdes)]=edge1;
				int source=comms[i].des + (G.NumOrigV * j);
				int destination=(comms[i].supdes);
				System.out.println("edges from des =" + source +","+ destination);
				e++;
				edge edges2 = new edge((comms[i].supsrc), (comms[i].src + (G.NumOrigV * j) ),
						10000, G.NumCopies, e, 0);
				//System.out.println(comms[i].src+ ","  + comms[i].des );
				G.tGraph[(comms[i].supsrc)][(comms[i].src + (G.NumOrigV * j) )]=edges2;
				e++;
				source=comms[i].supsrc;
				destination=comms[i].src + (G.NumOrigV * j);
				System.out.println("edges from src =" + source + ","+ destination);
			}
		}


		//lp equations
		try {
			IloCplex cplex = new IloCplex();
			IloMPModeler model=cplex;
			IloNumVar[][][] f = new IloNumVar[2*G.NumVert][2*G.NumVert][nocom];  
			IloNumVar FLOW = model.numVar(.01, Double.MAX_VALUE, IloNumVarType.Float, "flow");
			model.addMaximize(FLOW);
			IloRange[][][] rng=new IloRange[10][2*G.NumVert][2*G.NumVert];


			// variable initialization	
			for (int k=0;k<=G.NumVert*2;k++)
				for(int l=0;l<=G.NumVert*2;l++)
				{
					if(G.tGraph[k][l]!=null)
					{   
						for(int j=0;j<nocom;j++)
						{
							String xname="f" + "_" + G.tGraph[k][l].src + "_" +G.tGraph[k][l].dst + "_" + j ;
							f[k][l][j]=model.numVar(0, Double.MAX_VALUE, IloNumVarType.Float, xname);
						}
					}
				}

			//eq 9

			for (int k=0;k<G.NumVert;k++)
				for(int l=0;l<G.NumVert;l++)
				{
					if(G.tGraph[k][l]!=null)
					{   IloLinearNumExpr lin = cplex.linearNumExpr();
					for(int j=0;j<nocom;j++)
					{
						lin.addTerm(1, f[G.tGraph[k][l].src][G.tGraph[k][l].dst][j]);

					}
					rng[0][k][l]=model.addLe(lin, G.tGraph[k][l].cap);
					}
				}

			//eq 9 for supsrc and supdes


			for(int j=0;j<nocom;j++)
			{
				for(int i=0; i< G.NumCopies;i++)
				{
					IloLinearNumExpr lin = cplex.linearNumExpr();
					lin.addTerm(1, f[G.tGraph[(comms[j].des + (G.NumOrigV * i))][(comms[j].supdes)].src][G.tGraph[(comms[j].des + (G.NumOrigV * i))][(comms[j].supdes)].dst][j]);
					rng[0][(comms[j].des + (G.NumOrigV * i))][(comms[j].supdes)]=model.addLe(lin, G.tGraph[(comms[j].des + (G.NumOrigV * i))][(comms[j].supdes)].cap);    	
					IloLinearNumExpr lin1 = cplex.linearNumExpr();
					lin1.addTerm(1, f[G.tGraph[(comms[j].supsrc)][(comms[j].src + (G.NumOrigV * i) )].src][G.tGraph[(comms[j].supsrc)][(comms[j].src + (G.NumOrigV * i) )].dst][j]);
					rng[0][(comms[j].supsrc)][(comms[j].src + (G.NumOrigV * i) )]=model.addLe(lin1, G.tGraph[(comms[j].supsrc)][(comms[j].src + (G.NumOrigV * i) )].cap);
				}
			}




			//eq 10
			boolean check=false;
			boolean check1=false;


			for(int i=0;i<nocom;i++)
			{
				for(int k=0; k<G.NumVert*2; k++)
				{
					IloLinearNumExpr lin = cplex.linearNumExpr();
					for(int l=0;l<G.NumVert*2;l++)
					{

						if((G.tGraph[k][l]!=null)&&(k!=comms[i].supsrc))
						{
							check=true;
							lin.addTerm(1, f[G.tGraph[k][l].src][G.tGraph[k][l].dst][i]);
							//out.print(" + f"+ "_" + G.tGraph[k][l].src + "_" +G.tGraph[k][l].dst + "_" + i );	
						}
						if((G.tGraph[l][k]!=null)&&(k!=comms[i].supdes) )
						{
							check1=true;	
							lin.addTerm(-1, f[G.tGraph[l][k].src][G.tGraph[l][k].dst][i]);
							//	out.print(" - f"+ "_" + G.tGraph[l][k].src + "_" +G.tGraph[l][k].dst + "_" + i);	
						}

					}
					if( check1|| check)
						rng[1][k][0]=model.addEq(lin, 0);

					check=false;
					check1=false;
				}	
			}

			//eq 11
			for(int i=0;i<nocom;i++)
			{ 
				IloLinearNumExpr lin = cplex.linearNumExpr();
				for(int k=0; k<G.NumVert; k++)
				{

					if(G.tGraph[k][comms[i].supdes]!=null)
					{
						lin.addTerm(1, f[ G.tGraph[k][comms[i].supdes].src ][G.tGraph[k][comms[i].supdes].dst][i]);

					}
				}

				rng[2][i][0]=(IloRange) model.addEq(lin, model.prod(comms[i].req, FLOW));

			}



			cplex.exportModel("lpex1.lp");




			// solve the model and display the solution if one was found
			if ( cplex.solve() ) {
				x     = cplex.getValue(FLOW);
				double[][][] flo = new double [2*G.NumVert][2*G.NumVert][nocom]; 

				cplex.output().println("Solution status = " + cplex.getStatus());
				cplex.output().println("Solution value  = " + cplex.getObjValue());
				cplex.output().println("flow value =" + x);

				for (int k=0;k<=G.NumVert*2;k++)
					for(int l=0;l<=G.NumVert*2;l++)
					{
						if(G.tGraph[k][l]!=null)
						{   
							for(int j=0;j<nocom;j++)
							{
								String xname="f" + "_" + G.tGraph[k][l].src + "_" +G.tGraph[k][l].dst + "_" + j ;
								flo[k][l][j]=0;
								flo[k][l][j]=cplex.getValue(f[k][l][j]);

								G.tGraph[k][l].flow[j]=flo[k][l][j];
								if(G.tGraph[k][l].flow[j]!=0)
									cplex.output().println(" value of " + xname +" =" + G.tGraph[k][l].flow[j]);
							}
						}
					}



			}
			cplex.end();



		} catch (IloException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		path.getpaths(comms, G, nocom, x);

	}
}
